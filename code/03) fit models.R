library(tidyverse)
library(brms)
# load data
merged_d2 = readRDS("data/merged_d2.rds")

# fit model
# hg4 = readRDS(file = "models/hg4.rds")
# 
# hg4_taxon <- brm(bf(conc_ppb_s ~ type_taxon + (1 + type_taxon|site) + (1 + type_taxon|pfas_type),
#               hu ~ type_taxon + (1 + type_taxon|site) + (1 + type_taxon|pfas_type)),
#            data = merged_d2, 
#            family = hurdle_gamma(link = "log"),
#            prior=c(prior(normal(0,1), class = Intercept), # -5,2
#                    prior(normal(0,0.5), class = b), # was (0,2), then (0,0.1), then (0,0.5)
#                    prior(normal(0,1), class = b, dpar = hu),
#                    prior(normal(-1.5, 1), class = Intercept, dpar= hu)), # was (-2,1), then (-1.5,1)
#            iter = 2000 , chains = 4, cores = 4,
#            seed = 5, 
#            # save_pars = save_pars(all = TRUE),
#            data2 = list(max_conc_ppb = unique(merged_d2$max_conc)))
# 
# saveRDS(hg4_taxon, file = "models/hg4_taxon.rds")

# hg4_taxon = readRDS(file = "models/hg4_taxon.rds")

hg4_taxon = update(hg4_taxon, newdata = merged_d2, iter = 2000 , chains = 4, cores = 4)
hg4_taxon_unfiltered = update(hg4_taxon, 
                              newdata = readRDS("data/merged_d2_unfiltered.rds"), 
                              iter = 500 , chains = 1, cores = 4,
                              data2 = list(merged_d2_unfiltered = readRDS("data/merged_d2_unfiltered.rds")))

# sum pfas ----------------------------------------------------------------
merged_d = readRDS("data/full_data.rds") %>% 
  separate(sample_type, into = c('type', 'taxon'))

d2_sum = merged_d2 %>%
  group_by(sample_id, site, type) %>% 
  reframe(sum_ppb = sum(conc_ppb)) %>% 
  mutate(mean_sum_ppb = mean(sum_ppb, na.rm = T),
         sum_ppb_s = sum_ppb/mean_sum_ppb,
         sum_ppb_s_01 = sum_ppb_s + 0.0001)


mod1 = brm(sum_ppb_s_01 ~ type + (1 + type|site),
           family = Gamma(link = "log"),
           prior = c(prior(normal(0, 1), class = "Intercept"),
                     prior(normal(0, 1), class = "b"),
                     prior(exponential(2), class = "sd")),
           data = d2_sum, cores = 4, 
           data2 = list(raw_data = d2_sum))

saveRDS(mod1, file = "models/mod1.rds")


d2_sum_taxa = merged_d2 %>%
  filter(!is.na(taxon)) %>% 
  group_by(sample_id, site, type, taxon) %>% 
  reframe(sum_ppb = sum(conc_ppb)) %>% 
  mutate(mean_sum_ppb = mean(sum_ppb, na.rm = T),
         sum_ppb_s = sum_ppb/mean_sum_ppb,
         sum_ppb_s_01 = sum_ppb_s + 0.0001)


mod1_taxa = update(mod1, 
                   formula = sum_ppb_s_01 ~ type + (1 + type|site) + (1 + type|taxon),
                   newdata = d2_sum_taxa,
           data2 = list(raw_data = d2_sum_taxa))

saveRDS(mod1_taxa, file = "models/mod1_taxa.rds")


# log_kmw -----------------------------------------------------------------
# fits baf vs log_kmw. This data is created in 05) plots and tables...
posts_taxon_kmw_mean = readRDS(file = "posteriors/posts_taxon_kmw.rds") %>% 
  # group_by(site, name) %>% 
  # mutate(max_value = max(value),
  #        value = value/max_value) %>% 
  # group_by(log_kmw_mean, log_kmw_sd, name, draw) %>% 
  # reframe(value = mean(value, na.rm = T)) %>% 
  group_by(name, site) %>%
  mutate(mean_value = mean(value),
         value = value/mean_value) %>% 
  ungroup() %>% 
  group_by(log_kmw_mean, log_kmw_sd, name, site, mean_value) %>% 
  reframe(baf_mean = mean(value),
          baf_sd = sd(value)) %>% 
  mutate(log_kmw_mean_s = scale(log_kmw_mean),
         baf_site = paste(name, site, sep = "_")) 

posts_taxon_kmw_mean %>% 
  ggplot(aes(x = log_kmw_mean_s, y = baf_mean, color = site)) + 
  geom_pointrange(aes(ymin = baf_mean - baf_sd, ymax = baf_mean + baf_sd)) +
  facet_wrap(~name, scales = "free") +
  scale_y_log10() +
  geom_smooth() +
  NULL

brm_kmw = brm(baf_mean ~ s(log_kmw_mean_s, by = name) + (1|site),
              chains = 4, 
              iter = 2000,
              cores = 4,
              family = Gamma(link = "log"),
              # prior = c(prior(normal(2.7, 2), class = "Intercept")),
              data = posts_taxon_kmw_mean,
              data2 = list(mean_kmw = attributes(posts_taxon_kmw_mean$log_kmw_mean_s)$'scaled:center',
                           sd_kmw = attributes(posts_taxon_kmw_mean$log_kmw_mean_s)$'scaled:scale', 
                           mean_y = posts_taxon_kmw_mean %>% ungroup %>% distinct(mean_value, name, site, baf_site)))

saveRDS(brm_kmw, file = "models/brm_kmw.rds")

pp_check(brm_kmw) + scale_x_log10()



# insect mass -------------------------------------------------------------

# load models and data
insect_mass = read_excel("data/Farmington_InsectMasses.xlsx") %>% clean_names() %>% 
  mutate(units = "grams") %>% 
  mutate(gdw = 0.2*composite_sample_mass_ww) 

mean_insect_mass = insect_mass %>% 
  group_by(site, order, life_stage) %>% 
  reframe(mean_gdw = mean(0.2*composite_sample_mass_ww, na.rm = T), # 0.2 converts wet to dry mass
          sd_gdw = sd(0.2*composite_sample_mass_ww, na.rm = T)) %>% 
  rename(taxon = order) %>% 
  mutate(type = case_when(life_stage == "Adult" ~ "Emergent", T ~ life_stage),
         type_taxon = paste0(type, "_", taxon)) %>% 
  mutate(site = case_when(site == "Pequabuck Brook" ~ "Pequabuck River", TRUE ~ site))

mean_insect_mass_pertaxon = insect_mass %>% 
  group_by(order, life_stage) %>% 
  reframe(mean_gdw = median(0.2*composite_sample_mass_ww, na.rm = T), # 0.2 converts wet to dry mass
          sd_gdw = sd(0.2*composite_sample_mass_ww, na.rm = T)) %>% 
  rename(taxon = order) %>% 
  mutate(type = case_when(life_stage == "Adult" ~ "Emergent", T ~ life_stage),
         type_taxon = paste0(type, "_", taxon)) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>% 
  mutate(mean_sd = paste(mean_gdw, " \u00b1 ", sd_gdw))

write_csv(mean_insect_mass_pertaxon, file = "tables/mean_insect_mass_pertaxon.csv")

hg4_taxon = readRDS(file = "models/hg4_taxon.rds")
merged_d2 = hg4_taxon$data2$merged_d2

max_conc = unique(merged_d2$max_conc)

# plot insect masses
mean_insect_mass %>% 
  ggplot(aes(x = taxon, color = life_stage,
             y = mean_gdw)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_gdw - sd_gdw, ymax = mean_gdw + sd_gdw)) +
  facet_wrap(~site)

mod_mass = brm(gdw ~ order + (1 + order|site) + (1 + order|life_stage),
               data = insect_mass,
               family = Gamma(link = "log"),
               prior = c(prior(normal(0, 1), class = b),
                         prior(exponential(2), class = sd)))
saveRDS(mod_mass, file = "models/mod_mass.rds")

mean_insect_mass = insect_mass %>% 
  distinct(order, site, life_stage) %>% 
  add_epred_draws(mod_mass, re_formula = NULL) %>% 
  mutate(type = case_when(life_stage == "Adult" ~ "Emergent", T ~ life_stage)) %>% 
  mutate(site = case_when(site == "Pequabuck Brook" ~ "Pequabuck River", TRUE ~ site)) %>% 
  mutate(taxon = order,
         type_taxon = paste0(type, "_", taxon)) %>% 
  group_by(site, type, type_taxon) %>% 
  reframe(mean_gdw = median(.epred),
          sd_gdw = sd(.epred))

saveRDS(mean_insect_mass, file = "posteriors/mean_insect_mass.rds")
