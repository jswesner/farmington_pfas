library(tidyverse)
library(brms)
library(tidybayes)

hg4_taxon = readRDS(file = "models/hg4_taxon.rds")
old_hg4_taxon = readRDS(file = "models/old/hg4_taxon.rds")

merged_d2 = readRDS("data/merged_d2.rds") 

pp_check(hg4_taxon) + scale_x_log10()
# compare to the original hg4_taxon
# pp_check(readRDS(file = "models/old/hg4_taxon.rds")) + scale_x_log10()

post = hg4_taxon$data %>% 
  select(-conc_ppb_s) %>% 
  distinct() %>% 
  left_join(merged_d2 %>% ungroup %>% distinct(type_taxon, max_conc)) %>%
  # mutate(max_conc = hg4_taxon$data2$max_conc_ppb) %>% 
  add_epred_draws(hg4_taxon) %>% 
  group_by(type_taxon, site) %>%
  mutate(.epred = .epred*max_conc) %>% 
  median_qi(.epred)

post %>% 
  ggplot(aes(x = type_taxon, y = .epred)) +
  geom_pointrange(aes(y = .epred + 1e-05, ymin = .lower, ymax = .upper)) +
  facet_wrap(~site) +
  scale_y_log10(limits = c(1e-05, 1e+03)) +
  geom_point(data = merged_d2, aes(y = conc_ppb)) +
  theme(axis.text.x = element_text(angle = 90))



# check influence of megaloptera ------------------------------------------

pfas_orders = tibble(pfas_type = c(
  "PFHxA", "PFHpA" , "PFOA" , "PFNA" ,
  "PFDA", "PFUnA", "PFDoA" ,
  "PFBS" , "PFHxS" , "PFOS" ,
  "6:2FTS" , "8:2FTS")) %>% 
  mutate(order = 1:nrow(.))

# load model
hg4_taxon = readRDS(file = "models/hg4_taxon.rds")
merged_d2 = hg4_taxon$data2$merged_d2

max_conc = unique(merged_d2$max_conc)

pfas_names = read_csv("data/log_kw.csv") %>% 
  mutate(pfas_type = str_remove(pfas_type, " ")) %>% 
  mutate(log_kmw_mean = parse_number(str_sub(log_kmw, 1, 4)),
         log_kmw_sd = parse_number(str_sub(log_kmw, -4, -1)),
         log_kpw_mean = parse_number(log_kpw),
         number_carbons = parse_number(no_carbons)) %>% 
  rename(pfas_category = type) %>% 
  distinct(pfas_category, pfas_type, number_carbons) %>% 
  filter(pfas_type %in% unique(merged_d2$pfas_type)) %>% 
  arrange(pfas_category, number_carbons) %>% 
  mutate(pfas_order = 1:nrow(.))

mod_dat = hg4_taxon$data %>% 
  separate(type_taxon, into = c("type", "taxon"), remove = F) %>% 
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8)) %>% 
  mutate(conc_ppb = conc_ppb_s*max_conc) %>%
  left_join(pfas_names) %>% 
  mutate(pfas_type = reorder(pfas_type, pfas_order))

posts_concentrations = mod_dat %>%
  select(-contains("conc_ppb")) %>%
  distinct() %>%
  # mutate(site = "new") %>% 
  add_epred_draws(seed = 20202, hg4_taxon, re_formula = NULL, dpar = T, ndraws = 500) %>%
  mutate(.epred = .epred*unique(hg4_taxon$data2$max_conc_ppb)) 

posts_withmegaloptera = posts_concentrations %>% 
  group_by(type, site, pfas_type, order, .draw) %>% # average over taxa
  reframe(.epred = mean(.epred)) %>%
  left_join(pfas_names) %>%
  mutate(pfas_type = reorder(pfas_type, pfas_order)) %>%
  group_by(pfas_type, order, type, site) %>%
  median_qi(.epred) %>% 
  left_join(pfas_names) %>% 
  mutate(pfas_type = reorder(pfas_type, pfas_order)) %>% 
  filter(type %in% c("Emergent", 'Larval'))

posts_withoutmegaloptera = posts_concentrations %>% 
  filter(taxon != "Megaloptera") %>% 
  group_by(type, site, pfas_type, order, .draw) %>% # average over taxa
  reframe(.epred = mean(.epred)) %>%
  left_join(pfas_names) %>%
  mutate(pfas_type = reorder(pfas_type, pfas_order)) %>%
  group_by(pfas_type, order, type, site) %>%
  median_qi(.epred) %>% 
  left_join(pfas_names) %>% 
  mutate(pfas_type = reorder(pfas_type, pfas_order)) %>% 
  filter(type %in% c("Emergent", 'Larval')) %>% 
  ungroup %>% 
  select(pfas_type, order, type, site, .epred, .lower, .upper) %>% 
  rename(.epred_without = .epred,
         .lower_without = .lower,
         .upper_without = .upper)

check_megaloptera = posts_withmegaloptera %>% 
  left_join(posts_withoutmegaloptera) %>% 
  filter(type %in% c("Emergent", "Larval")) %>% 
  ggplot(aes(x = .epred, y = .epred_without, color = type)) +
  geom_pointrange(aes(xmin = .lower, xmax = .upper)) +
  geom_linerange(aes(ymin = .lower_without, ymax = .upper_without)) +
  labs(y = "ppb without megaloptera",
       x = "ppb with megaloptera") +
  facet_wrap(~pfas_type) +
  geom_abline() +
  scale_x_log10() +
  scale_y_log10()

ggsave(check_megaloptera, file = "plots/check_megaloptera.jpg", width = 6.5, height = 6.5)





# sum ppb -----------------------------------------------------------------

mod1_taxa = readRDS(file = "models/mod1_taxa.rds")


mod1_taxa_posts = mod1_taxa$data2 %>% 
  distinct(site, type, taxon) %>% 
  add_epred_draws(mod1_taxa) %>% 
  mutate(sum_ppb = (.epred - 0.0001)*unique(mod1_taxa$data2$mean_sum_ppb)) %>% 
  mutate(.epred = case_when(sum_ppb < 0.000101 ~ 0, TRUE ~ sum_ppb - 0.0001)) 

mod1_taxa_posts %>% 
  group_by(type, taxon, .draw) %>% 
  reframe(.epred = mean(.epred)) %>% 
  group_by(type, taxon) %>% 
  median_qi(.epred)






