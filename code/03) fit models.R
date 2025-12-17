library(tidyverse)
library(brms)
library(readxl)

#~~~~This code fits eight models and takes ~20 hours to run. It also requires and installation of rstan. 
#~~~~The model codes are silenced with "#". The fitted models can instead be loaded directly below:
hg4_taxon = readRDS(file = "models/hg4_taxon.rds")    # pfas ppb
mod1 = readDRS(file = "models/mod1.rds")              # sum pfas
mod1_taxa = readRDS(file = "models/mod1_taxa.rds")    # sum pfas per taxon
mod_mass = readRDS(file = "models/mod_mass.rds")      # insect mass
brm_isotopes_notadjustedformetamorphosis = readRDS(file = "models/brm_isotopes_notadjustedformetamorphosis.rds")  # isotopes for n15 only
brm_tmf_iso_notadjustedformetamorphosis = readRDS(file = "models/brm_tmf_iso_notadjustedformetamorphosis.rds")    # tmf per pfas
brm_tmf_iso_sum = readRDS("models/brm_tmf_iso_sum.rds")   # tmf for sum pfas
brm_isotopes = readRDS(file = "models/brm_isotopes.rds")  # isotopes c13 and n15 multivariate (bi-plot)

# pfas ppb ----------------------------------------------------------------
# load data
merged_d2 = readRDS("data/merged_d2.rds")

# fit model
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

hg4_taxon = readRDS(file = "models/hg4_taxon.rds")


# sum pfas ----------------------------------------------------------------
merged_d2_sum = readRDS(file = "data/merged_d2_sum.rds")


# mod1 = brm(sum_ppb_s_01 ~ type + (1 + type|site),
#            family = Gamma(link = "log"),
#            prior = c(prior(normal(0, 1), class = "Intercept"),
#                      prior(normal(0, 1), class = "b"),
#                      prior(exponential(2), class = "sd")),
#            data = d2_sum, cores = 4, 
#            data2 = merged_d2_sum)
# 
# saveRDS(mod1, file = "models/mod1.rds")


merged_d2_sum_taxa = readRDS(file = "data/merged_d2_sum_taxa.rds")

# mod1_taxa = update(mod1, 
#                    formula = sum_ppb_s_01 ~ type + (1 + type|site) + (1 + type|taxon),
#                    newdata = merged_d2_sum_taxa,
#            data2 = merged_d2_sum_taxa)
# 
# saveRDS(mod1_taxa, file = "models/mod1_taxa.rds")


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

# isotopes N15 only -------------------------------------------------------

# load data
isotopes = readRDS(file = "data/isotopes.rds")

# brm_isotopes_notadjustedformetamorphosis = brm(mean_centered_15n ~ taxon + (1 + taxon|site),
#                                                family = gaussian(),
#                                                prior = c(prior(normal(0, 1), class = Intercept),
#                                                          prior(normal(0, 1), class = b),
#                                                          prior(exponential(2), class = sd)),
#                                                data = isotopes)
# 
# saveRDS(brm_isotopes_notadjustedformetamorphosis, file = "models/brm_isotopes_notadjustedformetamorphosis.rds")
# plot(conditional_effects(brm_isotopes_notadjustedformetamorphosis), points = T)

# posteriors 
brm_isotopes_notadjustedformetamorphosis = readRDS(file = "models/brm_isotopes_notadjustedformetamorphosis.rds")

baseline_n15 = readRDS(file = "posteriors/iso_posts_notadjustedformetamorphosis.rds") %>% 
  filter(taxon == "Biofilm") %>% 
  mutate(baseline_n15_raw = .epred + center_15n) %>% 
  ungroup %>% 
  select(site, .draw, baseline_n15_raw) 

iso_posts_notadjustedformetamorphosis = isotopes %>%
  distinct(taxon, site, mean_15n) %>% 
  rename(center_15n = mean_15n) %>% 
  add_epred_draws(brm_isotopes_notadjustedformetamorphosis, re_formula = ~ NULL) %>% 
  group_by(taxon) %>% 
  mutate(median = median(.epred)) %>% 
  left_join(baseline_n15) %>% 
  mutate(raw_n15 = .epred + center_15n,
         trophic_level = 1 + (raw_n15 - baseline_n15_raw)/3.4)

saveRDS(iso_posts_notadjustedformetamorphosis, file = "posteriors/iso_posts_notadjustedformetamorphosis.rds")

# tmf per pfas ------------------------------------------------
iso_posts_notadjustedformetamorphosis = readRDS(file = "posteriors/iso_posts_notadjustedformetamorphosis.rds") 

pfas_posts = merged_d2 %>% 
  distinct(site, type_taxon, taxon, type, max_conc, pfas_type) %>% 
  filter(type != "Sediment") %>% 
  filter(type != "Detritus") %>%
  filter(type != "Water") %>% 
  filter(type != "Seston") %>% 
  mutate(taxon = case_when(is.na(taxon) ~ type,TRUE ~ taxon)) %>% 
  add_epred_draws(hg4_taxon, re_formula = NULL) %>% 
  mutate(.epred = .epred*max_conc) 

iso_post_summaries_notadjustedformetamorphosis = iso_posts_notadjustedformetamorphosis %>% 
  group_by(taxon, site, center_15n) %>% 
  reframe(mean_n15 = mean(.epred),
          median_n15 = median(.epred),
          sd_n15 = sd(.epred),
          lower_n15 = quantile(.epred, probs = 0.025),
          upper_n15 = quantile(.epred, probs = 0.975)) %>% 
  filter(site != "Russell Brook") # No spiders and only 2 dots for insects

saveRDS(iso_post_summaries_notadjustedformetamorphosis, file = "posteriors/iso_post_summaries_notadjustedformetamorphosis.rds")

pfas_conc_isotopes_notadjustedformetamorphosis = pfas_posts %>% 
  filter(site != "Russell Brook") %>% # only 2 data points. Can fit regression here and model below was having trouble so I removed.
  group_by(site, taxon, type, pfas_type)  %>% 
  reframe(mean_conc = mean(.epred),
          median_conc = median(.epred),
          sd_conc = sd(.epred),
          lower_conc = quantile(.epred, probs = 0.025),
          upper_conc = quantile(.epred, probs = 0.975)) %>% 
  mutate(taxon = case_when(taxon == "Tetragnathidae" ~ "Spider", TRUE ~ taxon)) %>%
  left_join(iso_post_summaries_notadjustedformetamorphosis, relationship = "many-to-many") %>% 
  mutate(log10_mean_conc = log10(mean_conc),
         log10_sd_conc = sd(log10_mean_conc),
         log10_mean_conc_s = scale(log10_mean_conc),
         log10_median_conc = log10(median_conc),
         log10_median_conc_s = scale(log10_median_conc)) 

saveRDS(pfas_conc_isotopes_notadjustedformetamorphosis, "data/pfas_conc_isotopes_notadjustedformetamorphosis.rds")

pfas_conc_isotopes_notadjustedformetamorphosis %>%
  # filter(site == "Burr Pond Brook") %>%
  # filter(pfas_type == "PFOS") %>% 
  ggplot(aes(x = mean_n15, color = pfas_type, y = log10_median_conc)) + 
  geom_point() +
  # geom_text(aes(label = taxon)) +
  facet_grid2(site ~ pfas_type, scales = "free_y") +
  scale_color_custom() +
  # scale_y_log10() +
  # scale_x_log10() +
  guides(fill = "none",
         color = "none") +
  # geom_smooth(method = lm, se = T) +
  NULL

# fit model

# brm_tmf_iso_notadjustedformetamorphosis = brm(log10_median_conc_s ~ mean_n15 + (1 + mean_n15|site) + (1 + mean_n15|pfas_type),
#                   data = pfas_conc_isotopes_notadjustedformetamorphosis,
#                   data2 = pfas_conc_isotopes_notadjustedformmetamorphosis,
#                   family = gaussian(),
#                   prior = c(prior(normal(0, 1), class = Intercept),
#                              prior(normal(0, 1), class = b),
#                              prior(exponential(2), class = sd)),
#                   iter = 2000, chains = 4)
# 
# saveRDS(brm_tmf_iso_notadjustedformetamorphosis, file = "models/brm_tmf_iso_notadjustedformetamorphosis.rds")

brm_tmf_iso_notadjustedformetamorphosis = readRDS(file = "models/brm_tmf_iso_notadjustedformetamorphosis.rds")
# tmf sum -----------------------------------------------------

iso_post_summaries_notadjustedformetamorphosis = readRDS(file = "posteriors/iso_post_summaries_notadjustedformetamorphosis.rds")

mod1_taxa = readRDS(file = "models/mod1_taxa.rds") # has just the insect taxa
mod1 = readRDS(file = "models/mod1.rds") # has biofilm and spiders 

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
         log10_median_sum_s = scale(log10_median_sum)) %>% 
  filter(site != "Russell Brook") # No spiders and only 2 dots for insects

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

brm_tmf_iso_sum = readRDS("models/brm_tmf_iso_sum.rds")

# isotopes C13 N15 multivariate ----------------------------------------------------------------
isotopes = readRDS("data/isotopes.rds")

# filter for just spiders and emergers. Remove Russell Brook b/c it only has 1 emerger value and no spiders measured
isotopes_filtered = isotopes %>% 
  filter(sample_type %in% c("Spider", "Emergent Diptera", "Emergent Ephemeroptera", "Emergent Odonata", "Emergent Plecoptera",
                            "Emergent Trichoptera")) %>% 
  filter(site != "Russell Brook") %>% 
  mutate(d13c = mean_centered_13c,
         d15n = mean_centered_15n)

# fit_model
# brm_isotopes = brm(bf(mvbind(d13c, d15n) ~ sample_type + (1 + sample_type|site),
#               set_rescor(TRUE)),
#            data = isotopes_filtered,
#            prior = c(prior(normal(0,1), class = Intercept, resp = c("d13c", "d15n")),
#                      prior(normal(0,1), class = b, resp = c("d13c", "d15n")),
#                      prior(exponential(2), class = sigma, resp = c("d13c", "d15n"))),
#            chains = 4,
#            cores = 4)
# 
# saveRDS(brm_isotopes, file = "models/brm_isotopes.rds")
brm_isotopes = readRDS(file = "models/brm_isotopes.rds")

# filter for just larvae. Remove Russell Brook b/c it only has 1 emerger value and no spiders measured
isotopes_notfiltered = isotopes %>%
  filter(site != "Russell Brook") %>% 
  mutate(d13c = mean_centered_13c,
         d15n = mean_centered_15n)

# fit_model
# brm_isotopes_notfiltered = update(brm_isotopes, newdata = isotopes_notfiltered,data2 = isotopes_notfiltered)
# # 
# saveRDS(brm_isotopes_notfiltered, file = "models/brm_isotopes_notfiltered.rds")
brm_isotopes_notfiltered = readRDS(file = "models/brm_isotopes_notfiltered.rds")

