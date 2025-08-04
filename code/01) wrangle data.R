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

merged_d2 = read_csv("data/Farmington_PFAS_FWstudy_Datarelease.csv") %>% 
  clean_names() %>% 
  rename(sample_id = field_id,
         pfas_type = compound) %>% 
  left_join(log_kmw %>% select(pfas_type, log_kmw_mean, log_kmw_sd, log_kpw)) %>% 
  mutate(nondetects = case_when(conc == "ND" ~ "ND", # not detected
                                conc == "NM" ~ "NM", # not measured
                                conc == 'na' ~ 'na')) %>%  # not available
  mutate(conc_ppb = parse_number(conc)) %>% 
  mutate(conc_ppb = case_when(nondetects == "ND" ~ 0, TRUE ~ conc_ppb)) %>% 
  filter(pfas_type == 'X6.2FTS' | pfas_type == 'PFBS' | pfas_type == 'PFDA' | pfas_type == 'PFDoA' | pfas_type == 'PFHpA' | 
           pfas_type == 'PFHxA' | pfas_type == 'PFHxS' | pfas_type == 'PFNA' | pfas_type == 'PFOA' | pfas_type == 'PFOS' | 
           pfas_type == 'PFUnA' | pfas_type == 'X8.2FTS') %>%
  droplevels() %>% 
  mutate(pfas_type = case_when(pfas_type == "X6.2FTS" ~ "6:2FTS",
                               pfas_type == "X8.2FTS" ~ "8:2FTS",
                               TRUE ~ pfas_type)) %>% 
  filter(!grepl("lank", sample_id)) %>% 
  mutate(log_kmw_mean_s = scale(log_kmw_mean)) %>% 
  mutate(type_taxon = sample_type) %>% 
  separate(sample_type, into = c('type', 'taxon')) %>% 
  mutate(max_conc_ppb = max(conc_ppb, na.rm = T),
         max_conc = max_conc_ppb) %>% 
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



# sum pfas ----------------------------------------------------------------

merged_d2_sum = read_csv("data/Farmington_PFAS_FWstudy_Datarelease.csv") %>% 
  clean_names() %>% 
  rename(sample_id = field_id,
         pfas_type = compound) %>% 
  left_join(log_kmw %>% select(pfas_type, log_kmw_mean, log_kmw_sd, log_kpw)) %>% 
  mutate(nondetects = case_when(conc == "ND" ~ "ND", # not detected
                                conc == "NM" ~ "NM", # not measured
                                conc == 'na' ~ 'na')) %>%  # not available
  mutate(conc_ppb = parse_number(conc)) %>% 
  mutate(conc_ppb = case_when(nondetects == "ND" ~ 0, TRUE ~ conc_ppb)) %>% 
  # filter(pfas_type == 'X6.2FTS' | pfas_type == 'PFBS' | pfas_type == 'PFDA' | pfas_type == 'PFDoA' | pfas_type == 'PFHpA' | 
  #          pfas_type == 'PFHxA' | pfas_type == 'PFHxS' | pfas_type == 'PFNA' | pfas_type == 'PFOA' | pfas_type == 'PFOS' | 
  #          pfas_type == 'PFUnA' | pfas_type == 'X8.2FTS') %>%
  droplevels() %>% 
  mutate(pfas_type = case_when(pfas_type == "X6.2FTS" ~ "6:2FTS",
                               pfas_type == "X8.2FTS" ~ "8:2FTS",
                               TRUE ~ pfas_type)) %>% 
  filter(!grepl("lank", sample_id)) %>% 
  mutate(log_kmw_mean_s = scale(log_kmw_mean)) %>% 
  mutate(type_taxon = sample_type) %>% 
  separate(sample_type, into = c('type', 'taxon')) %>% 
  mutate(max_conc_ppb = max(conc_ppb, na.rm = T)) %>% 
  mutate(conc_ppb_s = conc_ppb/max(conc_ppb),
         max_conc = max_conc_ppb) %>% 
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8)) %>%
  group_by(sample_id, site, type) %>% 
  reframe(sum_ppb = sum(conc_ppb)) %>% 
  mutate(mean_sum_ppb = mean(sum_ppb, na.rm = T),
         sum_ppb_s = sum_ppb/mean_sum_ppb,
         sum_ppb_s_01 = sum_ppb_s + 0.0001)

saveRDS(merged_d2_sum, file = "data/merged_d2_sum.rds")


merged_d2_sum_taxa = read_csv("data/Farmington_PFAS_FWstudy_Datarelease.csv") %>% 
  clean_names() %>% 
  rename(sample_id = field_id,
         pfas_type = compound) %>% 
  left_join(log_kmw %>% select(pfas_type, log_kmw_mean, log_kmw_sd, log_kpw)) %>% 
  mutate(nondetects = case_when(conc == "ND" ~ "ND", # not detected
                                conc == "NM" ~ "NM", # not measured
                                conc == 'na' ~ 'na')) %>%  # not available
  mutate(conc_ppb = parse_number(conc)) %>% 
  mutate(conc_ppb = case_when(nondetects == "ND" ~ 0, TRUE ~ conc_ppb)) %>% 
  # filter(pfas_type == 'X6.2FTS' | pfas_type == 'PFBS' | pfas_type == 'PFDA' | pfas_type == 'PFDoA' | pfas_type == 'PFHpA' | 
  #          pfas_type == 'PFHxA' | pfas_type == 'PFHxS' | pfas_type == 'PFNA' | pfas_type == 'PFOA' | pfas_type == 'PFOS' | 
  #          pfas_type == 'PFUnA' | pfas_type == 'X8.2FTS') %>%
  droplevels() %>% 
  mutate(pfas_type = case_when(pfas_type == "X6.2FTS" ~ "6:2FTS",
                               pfas_type == "X8.2FTS" ~ "8:2FTS",
                               TRUE ~ pfas_type)) %>% 
  filter(!grepl("lank", sample_id)) %>% 
  mutate(log_kmw_mean_s = scale(log_kmw_mean)) %>% 
  mutate(type_taxon = sample_type) %>% 
  separate(sample_type, into = c('type', 'taxon')) %>% 
  mutate(max_conc_ppb = max(conc_ppb, na.rm = T),
         max_conc = max_conc_ppb) %>% 
  mutate(conc_ppb_s = conc_ppb/max(conc_ppb)) %>% 
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8)) %>%
  filter(!is.na(taxon)) %>% 
  group_by(sample_id, site, type, taxon) %>% 
  reframe(sum_ppb = sum(conc_ppb)) %>% 
  mutate(mean_sum_ppb = mean(sum_ppb, na.rm = T),
         sum_ppb_s = sum_ppb/mean_sum_ppb,
         sum_ppb_s_01 = sum_ppb_s + 0.0001)


saveRDS(merged_d2_sum_taxa, file = "data/merged_d2_sum_taxa.rds")
