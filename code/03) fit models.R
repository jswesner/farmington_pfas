library(tidyverse)
library(brms)
library(readxl)


# pfas ppb ----------------------------------------------------------------


# load data
merged_d2 = readRDS("data/merged_d2.rds")

# fit model
# hg4 = readRDS(file = "models/hg4.rds")
# 
# hg4_taxon <- brm(bf(conc_ppb_s ~ type_taxon + (1 + type_taxon|site) + (1 + type_taxon|pfas_type),
#               hu ~ type_taxon + (1 + type_taxon|site) + (1 + type_taxon|pfas_type)),
#            data = merged_d2,
#            family = hurdle_gamma(link = "log"),
#            prior=c(prior(normal(-2,1), class = Intercept), # -5,2
#                    prior(normal(0,0.5), class = b), # was (0,2), then (0,0.1), then (0,0.5)
#                    prior(normal(0,1), class = b, dpar = hu),
#                    prior(normal(-1.5, 1), class = Intercept, dpar= hu)), # was (-2,1), then (-1.5,1)
#            iter = 2000 , chains = 4, cores = 4,
#            seed = 5,
#            # save_pars = save_pars(all = TRUE),
#            data2 = merged_d2))
# 
# saveRDS(hg4_taxon, file = "models/hg4_taxon.rds")

# hg4_taxon = readRDS(file = "models/hg4_taxon.rds")

hg4_taxon = update(hg4_taxon, newdata = merged_d2, iter = 2000 , chains = 4, cores = 4,
                   data2 = merged_d2)

hg4_taxon_unfiltered = update(hg4_taxon, 
                              newdata = readRDS("data/merged_d2_unfiltered.rds"), 
                              iter = 500 , chains = 1, cores = 4,
                              data2 = list(merged_d2_unfiltered = readRDS("data/merged_d2_unfiltered.rds")))

# sum pfas ----------------------------------------------------------------
merged_d2_sum = readRDS(file = "data/merged_d2_sum.rds")


mod1 = brm(sum_ppb_s_01 ~ type + (1 + type|site),
           family = Gamma(link = "log"),
           prior = c(prior(normal(0, 1), class = "Intercept"),
                     prior(normal(0, 1), class = "b"),
                     prior(exponential(2), class = "sd")),
           data = d2_sum, cores = 4, 
           data2 = merged_d2_sum)

saveRDS(mod1, file = "models/mod1.rds")


merged_d2_sum_taxa = readRDS(file = "data/merged_d2_sum_taxa.rds")

mod1_taxa = update(mod1, 
                   formula = sum_ppb_s_01 ~ type + (1 + type|site) + (1 + type|taxon),
                   newdata = d2_sum_taxa,
           data2 = merged_d2_sum_taxa)

saveRDS(mod1_taxa, file = "models/mod1_taxa.rds")


# insect mass -------------------------------------------------------------

# load models and data
insect_mass = readRDS(file = "data/insect_mass.rds")

# mod_mass = brm(gdw ~ order + (1 + order|site) + (1 + order|life_stage),
#                data = insect_mass,
#                family = Gamma(link = "log"),
#                prior = c(prior(normal(0, 1), class = b),
#                          prior(exponential(2), class = sd)))
# saveRDS(mod_mass, file = "models/mod_mass.rds")

mod_mass = readRDS(file = "models/mod_mass.rds")

insect_mass_posts = insect_mass %>% 
  distinct(order, site, life_stage) %>% 
  add_epred_draws(mod_mass, re_formula = NULL) %>% 
  mutate(type = case_when(life_stage == "Adult" ~ "Emergent", T ~ life_stage)) %>% 
  mutate(site = case_when(site == "Pequabuck Brook" ~ "Pequabuck River", TRUE ~ site)) %>% 
  mutate(taxon = order,
         type_taxon = paste0(type, "_", taxon)) %>% 
  rename(mean_gww = .epred)

saveRDS(insect_mass_posts, file = "posteriors/insect_mass_posts.rds")

mean_insect_mass = insect_mass %>% 
  distinct(order, site, life_stage) %>% 
  add_epred_draws(mod_mass, re_formula = NULL) %>% 
  mutate(type = case_when(life_stage == "Adult" ~ "Emergent", T ~ life_stage)) %>% 
  mutate(site = case_when(site == "Pequabuck Brook" ~ "Pequabuck River", TRUE ~ site)) %>% 
  mutate(taxon = order,
         type_taxon = paste0(type, "_", taxon)) %>% 
  group_by(site, type, type_taxon) %>% 
  reframe(mean_gww = median(.epred),
          sd_gww = sd(.epred),
          lower = quantile(.epred, probs = 0.025),
          upper = quantile(.epred, probs = 0.975))

saveRDS(mean_insect_mass, file = "posteriors/mean_insect_mass.rds")

insect_mass_averaged_over_sites = insect_mass %>% 
  distinct(order, site, life_stage) %>% 
  add_epred_draws(mod_mass, re_formula = NULL) %>% 
  mutate(type = case_when(life_stage == "Adult" ~ "Emergent", T ~ life_stage)) %>% 
  mutate(site = case_when(site == "Pequabuck Brook" ~ "Pequabuck River", TRUE ~ site)) %>% 
  mutate(taxon = order,
         type_taxon = paste0(type, "_", taxon)) %>% 
  # group_by(type, type_taxon, .epred) %>% 
  # reframe(.epred = mean(.epred)) %>% 
  group_by(type, type_taxon) %>% 
  reframe(mean_gww = median(.epred),
          sd_gww = sd(.epred),
          lower = quantile(.epred, probs = 0.025),
          upper = quantile(.epred, probs = 0.975))

saveRDS(insect_mass_averaged_over_sites, file = "posteriors/insect_mass_averaged_over_sites.rds")

# tmf -----------------------------------------------------
mod1_taxa = readRDS(file = "models/mod1_taxa.rds") # has just the insect taxa
mod1 = readRDS(file = "models/mod1.rds") # has biofilm and spiders 

iso_posts_notadjustedformetamorphosis = readRDS(file = "posteriors/iso_posts_notadjustedformetamorphosis.rds") 

pfas_sum_trophic_insects = mod1_taxa$data2 %>% 
  distinct(type, taxon, site, mean_sum_ppb) %>% 
  add_epred_draws(mod1_taxa) %>% 
  mutate(.epred = (.epred - 0.0001)*mean_sum_ppb) %>% 
  group_by(type, taxon, site) %>% 
  reframe(mean_sum = mean(.epred),
          median_sum = median(.epred),
          sd_sum = sd(.epred),
          lower_sum = quantile(.epred, probs = 0.025),
          upper_sum = quantile(.epred, probs = 0.975)) 

pfas_sum_trophic_biofilm_spiders = mod1$data2 %>% 
  distinct(type, site, mean_sum_ppb) %>% 
  filter(type %in% c("Biofilm", "Tetragnathidae")) %>% 
  add_epred_draws(mod1) %>% 
  mutate(.epred = (.epred - 0.0001)*mean_sum_ppb) %>% 
  group_by(type, site) %>% 
  reframe(mean_sum = mean(.epred),
          median_sum = median(.epred),
          sd_sum = sd(.epred),
          lower_sum = quantile(.epred, probs = 0.025),
          upper_sum = quantile(.epred, probs = 0.975)) %>%
  rename(taxon = type) %>% 
  mutate(taxon = case_when(taxon == "Tetragnathidae" ~ "Spider", TRUE ~ taxon)) 

pfas_sum_trophic = bind_rows(pfas_sum_trophic_insects, pfas_sum_trophic_biofilm_spiders) %>% 
  left_join(iso_post_summaries_notadjustedformetamorphosis, relationship = "many-to-many") %>% 
  mutate(log10_mean_sum = log10(mean_sum),
         log10_sd_sum = sd(log10_mean_sum),
         log10_mean_sum_s = scale(log10_mean_sum),
         log10_median_sum = log10(median_sum),
         log10_median_sum_s = scale(log10_median_sum)) 

pfas_sum_trophic %>% 
  ggplot(aes(x = mean_n15, y = log10_median_sum_s)) + 
  geom_point()

# fit model

# brm_tmf_iso_sum = brm(log10_median_sum_s ~ mean_n15 + (1 + mean_n15|site),
#                       data = pfas_sum_trophic,
#                       family = gaussian(),
#                       prior = c(prior(normal(0, 1), class = Intercept),
#                                 prior(normal(0, 1), class = b),
#                                 prior(exponential(2), class = sd)),
#                       iter = 2000, chains = 4)
# 
# saveRDS(brm_tmf_iso_sum, file = "models/brm_tmf_iso_sum.rds")
