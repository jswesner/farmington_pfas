library(tidyverse)
library(brms)
library(tidybayes)
library(ggtext)
library(ggthemes)
library(patchwork)
library(cowplot)

hg4_taxon = readRDS(file = "models/hg4_taxon.rds")
mod1 = readRDS(file = "models/mod1.rds")
mod1_taxa = readRDS(file = "models/mod1_taxa.rds")
mod_mass = readRDS(file = "models/mod_mass.rds")
isotopes_n15 = readRDS(file = "models/brm_isotopes_notadjustedformetamorphosis.rds")
tmf_perpfas =  readRDS(file = "models/brm_tmf_iso_notadjustedformetamorphosis.rds")
brm_tmf_iso_sum = readRDS(file = "models/brm_tmf_iso_sum.rds")
brm_isotopes = readRDS(file = "models/brm_isotopes.rds")

# posterior predictive checks ---------------
pp_check_pfas_conc = pp_check(hg4_taxon) + scale_x_log10() + labs(x = "Standardized PFAS concentration (ppb/max(ppb))",
                                                                  subtitle = "a) PFAS ppb")
pp_check_sum_pfas = pp_check(mod1) + scale_x_log10() + labs(x = "Standardized \u221112 PFAS concentration (\u221112 ppb/max(\u221112 ppb))",
                                                            subtitle = "b) \u221112 PFAS ppb")
pp_check_sum_pfas_taxa = pp_check(mod1_taxa) + scale_x_log10()  + labs(x = "Standardized \u221112 PFAS concentration (\u221112 ppb/max(\u22111 2ppb))",
                                                                       subtitle = "c) \u221112 PFAS concentration per taxon")
pp_check_insect_mass = pp_check(mod_mass) + scale_x_log10() + labs(x = "Individual wet mass",
                                                                   subtitle = "d) Mass")
pp_check_isotopes_n15 = pp_check(isotopes_n15) + labs(x = expression(paste(delta^{15}, "N")),
                                                                   subtitle = expression(paste("e) ", delta^{15}, "N (centered)")))
pp_check_tmf_perpfas = pp_check(tmf_perpfas) + labs(x = "PFAS Concentration (log10 tissue ppb)",
                                                      subtitle = "f) TMS per PFAS")
pp_check_tmf_sum = pp_check(brm_tmf_iso_sum) + labs(x = "\u221112 PFAS Concentration (log10 tissue ppb)",
                                                    subtitle = "g) TMS \u221112 PFAS")
pp_check_isotopes_both_N15 = pp_check(brm_isotopes, resp = "d15n") + labs(x = expression(paste(delta^{15}, "N")))
pp_check_isotopes_both_C13 = pp_check(brm_isotopes, resp = "d13c") + labs(x = expression(paste(delta^{13}, "C")))
pp_check_isotopes_both = ((pp_check_isotopes_both_C13 + guides(color = "none")+ labs(subtitle = "h) Multivariate Stable Isotopes")) + pp_check_isotopes_both_N15) 


pp_a = pp_check_pfas_conc$data %>% 
  mutate(is_y_label = case_when(grepl("rep", is_y_label) ~ "*y*<sub>rep",
                                TRUE ~ "*y*")) %>% 
  ggplot(aes(x = value + 0.00001)) +
  geom_density(aes(group = rep_id, color = is_y_label)) +
  scale_x_log10() +
  theme(legend.text = element_markdown(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank())  +
  scale_color_colorblind() +
  labs(x = "Standardized PFAS concentration\n(ppb/max(ppb))",
       subtitle = "a) PFAS ppb",
       color = "")

pp_b = pp_check_sum_pfas$data %>% 
  mutate(is_y_label = case_when(grepl("rep", is_y_label) ~ "*y*<sub>rep",
                                TRUE ~ "*y*")) %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(group = rep_id, color = is_y_label)) +
  scale_x_log10() +
  theme(legend.text = element_markdown(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank()) +
  scale_color_colorblind() +
  labs(x = "Standardized \u221112 PFAS concentration\n(\u221112 ppb/max(\u221112 ppb))",
       subtitle = "b) \u221112 PFAS ppb")


pp_c = pp_check_sum_pfas_taxa$data %>% 
  mutate(is_y_label = case_when(grepl("rep", is_y_label) ~ "*y*<sub>rep",
                         TRUE ~ "*y*")) %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(group = rep_id, color = is_y_label)) +
  scale_x_log10() +
  theme(legend.text = element_markdown(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank()) +
  scale_color_colorblind() +
  labs(x = "Standardized \u221112 PFAS concentration\n(\u221112 ppb/max(\u22111 2ppb))",
       subtitle = "c) \u221112 PFAS concentration per taxon")

pp_d = pp_check_insect_mass$data %>% 
  mutate(is_y_label = case_when(grepl("rep", is_y_label) ~ "*y*<sub>rep",
                                TRUE ~ "*y*")) %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(group = rep_id, color = is_y_label)) +
  scale_x_log10() +
  theme(legend.text = element_markdown(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank())  +
  scale_color_colorblind() +
  labs(x = "Individual wet mass",
       subtitle = "d) Mass")

pp_e = pp_check_isotopes_n15$data  %>% 
  mutate(is_y_label = case_when(grepl("rep", is_y_label) ~ "*y*<sub>rep",
                                TRUE ~ "*y*")) %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(group = rep_id, color = is_y_label)) +
  theme(legend.text = element_markdown(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank())  +
  scale_color_colorblind() +
  labs(x = expression(paste(delta^{15}, "N")),
       subtitle = expression(paste("e) ", delta^{15}, "N (centered)")))
  
pp_f = pp_check_tmf_perpfas$data %>% 
  mutate(is_y_label = case_when(grepl("rep", is_y_label) ~ "*y*<sub>rep",
                                TRUE ~ "*y*")) %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(group = rep_id, color = is_y_label)) +
  theme(legend.text = element_markdown(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank())  +
  scale_color_colorblind() +
  labs(x = "PFAS Concentration (log10 tissue ppb)",
       subtitle = "f) TMS per PFAS")

pp_g = pp_check_tmf_sum$data %>% 
  mutate(is_y_label = case_when(grepl("rep", is_y_label) ~ "*y*<sub>rep",
                                TRUE ~ "*y*")) %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(group = rep_id, color = is_y_label)) +
  theme(legend.text = element_markdown(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank())  +
  scale_color_colorblind() +
  labs(x = "\u221112 PFAS Concentration (log10 tissue ppb)",
       subtitle = "g) TMS \u221112 PFAS")

pp_h1 = pp_check_isotopes_both_N15$data %>% 
  mutate(is_y_label = case_when(grepl("rep", is_y_label) ~ "*y*<sub>rep",
                                TRUE ~ "*y*")) %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(group = rep_id, color = is_y_label)) +
  theme(legend.text = element_markdown(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 8))  +
  scale_color_colorblind() + 
  labs(x = expression(paste(delta^{15}, "N")),
       subtitle = "h) Multivariate Stable Isotopes") +
  guides(color = "none") 

pp_h2 = pp_check_isotopes_both_C13$data %>% 
  mutate(is_y_label = case_when(grepl("rep", is_y_label) ~ "*y*<sub>rep",
                                TRUE ~ "*y*")) %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(group = rep_id, color = is_y_label)) +
  theme(legend.text = element_markdown(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 8))  +
  scale_color_colorblind() + 
  labs(x = expression(paste(delta^{13}, "C")))

pp_h = pp_h1 + pp_h2

pp_checks = plot_grid(pp_a + guides(color = "none") + theme(text = element_text(size = 8)), 
                      pp_b + guides(color = "none") + theme(text = element_text(size = 8)), 
                      pp_c + guides(color = "none") + theme(text = element_text(size = 8)), 
                      pp_d + guides(color = "none") + theme(text = element_text(size = 8)), 
                      pp_e + guides(color = "none") + theme(text = element_text(size = 8)), 
                      pp_f + guides(color = "none") + theme(text = element_text(size = 8)), 
                      pp_g + guides(color = "none") + theme(text = element_text(size = 8)), 
                      pp_h,
                      ncol = 2)

ggsave(pp_checks, file = "plots/ms_plots_tables/pp_checks.jpg",
       width = 6.5, height = 8)


# bayes p-value -----------------------------------------------------------

get_bayes_p <- function(model, response, group = NULL) {
  data <- model$data %>%
    add_predicted_draws(model) %>%
    mutate(diff = .prediction - !!sym(response),
           greater = if_else(diff > 0, 1, 0)) %>%
    ungroup()
  
  if (!is.null(group)) {
    data <- data %>%
      group_by(across(all_of(group))) %>%
      add_tally() %>%
      group_by(across(all_of(c(group, "n"))))
  } else {
    data <- data %>%
      add_tally() %>%
      group_by(n)
  }
  
  data %>%
    reframe(sum = sum(greater)) %>%
    mutate(bayes_p = sum / n)
}


hg4_bayes_p = get_bayes_p(model = hg4_taxon, response = "conc_ppb_s") %>% mutate(model = "PFAS ppb")
mod1_bayes_p = get_bayes_p(model = mod1, response = "sum_ppb_s_01") %>% mutate(model = "Sum PFAS ppb")
mod1_taxa_bayes_p = get_bayes_p(model = mod1_taxa, response = "sum_ppb_s_01") %>% mutate(model = "Sum PFAS ppb per taxon")
mod_mass_bayes_p = get_bayes_p(model = mod_mass, response = "gdw") %>% mutate(model = "Insect mass")
isotopes_n15_bayes_p = get_bayes_p(model = isotopes_n15, response = "mean_centered_15n") %>% mutate(model = "deltaN15")
tmf_perpfas_bayes_p = get_bayes_p(model = tmf_perpfas, response = "log10_median_conc_s") %>% mutate(model = "TMS per PFAS")
tmf_sumpfas_bayes_p = get_bayes_p(model = brm_tmf_iso_sum, response = "log10_median_sum_s") %>% mutate(model = "TMS sum PFAS")
biplot_n15_bayes_p = get_bayes_p(model = brm_isotopes, response = "d15n") %>% mutate(model = "Multivariate Stable Isotopes d15n")
biplot_c15_bayes_p = get_bayes_p(model = brm_isotopes, response = "d13c") %>% mutate(model = "Multivariate Stable Isotopes d13c")



bayes_p = bind_rows(hg4_bayes_p,
          mod1_bayes_p,
          mod1_taxa_bayes_p,
          mod_mass_bayes_p,
          isotopes_n15_bayes_p,
          tmf_perpfas_bayes_p,
          tmf_sumpfas_bayes_p,
          biplot_n15_bayes_p,
          biplot_c15_bayes_p)


write_csv(bayes_p, file = "tables/bayes_p.csv")

hg4_bayes_p_type = get_bayes_p(model = hg4_taxon, response = "conc_ppb_s", group = "type_taxon") %>% mutate(model = "PFAS ppb") %>% rename(group = type_taxon)
mod1_bayes_p_type = get_bayes_p(model = mod1, response = "sum_ppb_s_01", group = "type") %>% mutate(model = "Sum PFAS ppb") %>% rename(group = type)
mod1_taxa_bayes_p_type = get_bayes_p(model = mod1_taxa, response = "sum_ppb_s_01", group = "type") %>% mutate(model = "Sum PFAS ppb per taxon") %>% rename(group = type)
mod_mass_bayes_p_type = get_bayes_p(model = mod_mass, response = "gdw", group = "order") %>% mutate(model = "Insect mass") %>% rename(group = "order")
isotopes_n15_bayes_p_type = get_bayes_p(model = isotopes_n15, response = "mean_centered_15n", group = "taxon") %>% mutate(model = "deltaN15") %>% rename(group = taxon)
tmf_perpfas_bayes_p_type = get_bayes_p(model = tmf_perpfas, response = "log10_median_conc_s", group = "pfas_type") %>% mutate(model = "TMS per PFAS") %>% rename(group = pfas_type)
tmf_sumpfas_bayes_p_type = get_bayes_p(model = brm_tmf_iso_sum, response = "log10_median_sum_s", group = "site") %>% mutate(model = "TMS sum PFAS") %>% rename(group = site)
biplot_n15_bayes_p_type = get_bayes_p(model = brm_isotopes, response = "d15n", group = "sample_type") %>% mutate(model = "Multivariate Stable Isotopes d15n") %>% rename(group = sample_type)
biplot_c15_bayes_p_type = get_bayes_p(model = brm_isotopes, response = "d13c", group = "sample_type") %>% mutate(model = "Multivariate Stable Isotopes d13c") %>% rename(group = sample_type)

bayes_p_type = bind_rows(hg4_bayes_p_type,
                    mod1_bayes_p_type,
                    mod1_taxa_bayes_p_type,
                    mod_mass_bayes_p_type,
                    isotopes_n15_bayes_p_type,
                    tmf_perpfas_bayes_p_type,
                    tmf_sumpfas_bayes_p_type,
                    biplot_n15_bayes_p_type,
                    biplot_c15_bayes_p_type)

write_csv(bayes_p_type, file = "tables/bayes_p_type.csv")


# other -------------------------------------------------------------------


merged_d2 = readRDS("data/merged_d2.rds") 
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

pfas_names = read_csv("data/pfas_names.csv")

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








# isotopes ----------------------------------------------------------------

brm_isotopes = readRDS(file = "models/brm_isotopes.rds")

pp_check(brm_isotopes, resp = "d15n")
pp_check(brm_isotopes, resp = "d13c")
