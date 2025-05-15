library(brms)
library(tidyverse)
library(tidybayes)
library(janitor)

# load data
log_kmw = read_csv("data/log_kw.csv") %>% 
  mutate(pfas_type = str_remove(pfas_type, " ")) %>% 
  mutate(log_kmw_mean = parse_number(str_sub(log_kmw, 1, 4)),
         log_kmw_sd = parse_number(str_sub(log_kmw, -4, -1)),
         log_kpw_mean = parse_number(log_kpw))

merged_d <- read_csv("data/PFAS_merged_d.csv") %>% 
  mutate(Site = as.factor(Site),
         PFAS_type = as.factor(PFAS_type)) %>% 
  clean_names() %>% 
  left_join(log_kmw %>% select(pfas_type, log_kmw_mean, log_kmw_sd, log_kpw))

saveRDS(merged_d, file = "data/full_data.rds")

# determine which PFAS compounds were detected in >5 % of samples
xxx <- merged_d %>% 
  select(-pfas_type) %>% 
  select(-sample_type) %>% 
  select(-conc_ppb) %>% 
  select(-size_class) %>%
  group_by(sample_id)

xxx <- xxx[!duplicated(xxx), ] ## 258 samples total
xxx <- merged_d %>% 
  group_by(pfas_type) %>%  
  count(conc_ppb > 0)
xxx %>% 
  filter(`conc_ppb > 0` == "TRUE") %>% 
  summarize(percent = (n/258)*100)

# filter for these compounds
merged_d2 <- merged_d %>% 
  filter(pfas_type == '6:2FTS' | pfas_type == 'PFBS' | pfas_type == 'PFDA' | 
           pfas_type == 'PFDoA' | pfas_type == 'PFHpA' | pfas_type == 'PFHxA' | 
           pfas_type == 'PFHxS' | pfas_type == 'PFNA' | pfas_type == 'PFOA' | 
           pfas_type == 'PFOS' | pfas_type == 'PFUnA' |pfas_type == '8:2FTS') %>%
  droplevels() %>% 
  mutate(log_kmw_mean_s = scale(log_kmw_mean)) %>% 
  mutate(type_taxon = sample_type) %>% 
  separate(sample_type, into = c('type', 'taxon')) %>% 
  mutate(max_conc = max(conc_ppb, na.rm = T)) %>% 
  mutate(conc_ppb_s = conc_ppb/max(conc_ppb)) %>% 
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8))

saveRDS(merged_d2, file = "data/merged_d2.rds")

merged_d2_unfiltered <- merged_d %>% 
  # filter(pfas_type == '6:2FTS' | pfas_type == 'PFBS' | pfas_type == 'PFDA' | 
  #          pfas_type == 'PFDoA' | pfas_type == 'PFHpA' | pfas_type == 'PFHxA' | 
  #          pfas_type == 'PFHxS' | pfas_type == 'PFNA' | pfas_type == 'PFOA' | 
  #          pfas_type == 'PFOS' | pfas_type == 'PFUnA' |pfas_type == '8:2FTS') %>%
  droplevels() %>% 
  mutate(log_kmw_mean_s = scale(log_kmw_mean)) %>% 
  mutate(type_taxon = sample_type) %>% 
  separate(sample_type, into = c('type', 'taxon')) %>% 
  mutate(max_conc = max(conc_ppb, na.rm = T)) %>% 
  mutate(conc_ppb_s = conc_ppb/max(conc_ppb)) %>% 
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8))

saveRDS(merged_d2_unfiltered, file = "data/merged_d2_unfiltered.rds")
