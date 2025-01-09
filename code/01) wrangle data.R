library(brms)
library(tidyverse)
library(tidybayes)
library(here)
library(janitor)


# load data
merged_d <- read_csv(here("data/PFAS_merged_d.csv")) %>% 
  mutate(Site = as.factor(Site),
         PFAS_type = as.factor(PFAS_type)) %>% 
  clean_names()

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
  filter(pfas_type == '6:2FTS' | pfas_type == 'PFBS' | pfas_type == 'PFDA' | pfas_type == 'PFDoA' | pfas_type == 'PFHpA' | pfas_type == 'PFHxA' | pfas_type == 'PFHxS' | pfas_type == 'PFNA' | pfas_type == 'PFOA' | pfas_type == 'PFOS' | pfas_type == 'PFUnA') %>%
  droplevels()

saveRDS(merged_d2, file = "data/merged_d2.rds")
