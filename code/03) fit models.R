library(tidyverse)
library(brms)
# load data
merged_d2 = readRDS("data/merged_d2.rds") %>% 
  mutate(type_taxon = sample_type) %>% 
  separate(sample_type, into = c('type', 'taxon')) %>% 
  mutate(max_conc = max(conc_ppb, na.rm = T)) %>% 
  mutate(conc_ppb_s = conc_ppb/max(conc_ppb))

# fit model
hg4 = readRDS(file = "models/hg4.rds")

hg4_taxon <- brm(bf(conc_ppb_s ~ type_taxon + (1 + type_taxon|site) + (1 + type_taxon|pfas_type),
              hu ~ type_taxon + (1 + type_taxon|site) + (1 + type_taxon|pfas_type)),
           data = merged_d2, 
           family = hurdle_gamma(link = "log"),
           prior=c(prior(normal(0,1), class = Intercept), # -5,2
                   prior(normal(0,0.5), class = b), # was (0,2), then (0,0.1), then (0,0.5)
                   prior(normal(0,1), class = b, dpar = hu),
                   prior(normal(-1.5, 1), class = Intercept, dpar= hu)), # was (-2,1), then (-1.5,1)
           iter = 2000 , chains = 4, cores = 4,
           seed = 5, 
           # save_pars = save_pars(all = TRUE),
           data2 = list(max_conc_ppb = unique(merged_d2$max_conc)))

saveRDS(hg4_taxon, file = "models/hg4_taxon.rds")

hg4_taxon = readRDS(file = "models/hg4_taxon.rds")

# sum pfas ----------------------------------------------------------------
merged_d = readRDS("data/full_data.rds") %>% 
  separate(sample_type, into = c('type', 'taxon'))

d2_sum = merged_d %>% 
  group_by(site, type) %>% 
  reframe(sum_ppb = sum(conc_ppb)) %>%
  mutate(sum_ppb_1 = sum_ppb + 0.01) %>% 
  group_by(type) %>%
  mutate(mean_sum_ppb_1 = mean(sum_ppb_1, na.rm = T),
         mean_sum_ppb = mean(sum_ppb, na.rm = T),
         sum_ppb_s_1 = sum_ppb_1/mean_sum_ppb_1,
         sum_ppb_s = sum_ppb/mean_sum_ppb)


d2_sum %>% 
  ggplot(aes(x = type, y = sum_ppb_s_1)) + 
  geom_point() 

mod1 = brm(sum_ppb_s_1 ~ type + (1 + type|site),
           family = Gamma(link = "log"),
           prior = c(prior(normal(1, 0.5), class = "Intercept"),
                     prior(normal(0, 1), class = "b"),
                     prior(exponential(2), class = "sd")),
           data = d2_sum, cores = 4)

saveRDS(mod1, file = "models/mod1.rds")


d2_max = merged_d2 %>% 
  mutate(conc_ppb_s = conc_ppb/max(conc_ppb))

mod2 = brm(sum_ppb_s_1 ~ type + (1 + type|site),
           family = Gamma(link = "log"),
           prior = c(prior(normal(1, 0.5), class = "Intercept"),
                     prior(normal(0, 1), class = "b"),
                     prior(exponential(2), class = "sd")),
           data = d2_sum, cores = 4)