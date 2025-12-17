library(brms)
library(tidyverse)
library(tidybayes)
library(janitor)

# pfas --------------------------------------------------------------------

merged_d2 = read_csv("data/Farmington_PFAS_FWstudy_Datarelease.csv")  %>% 
  mutate(nondetects = case_when(conc == "ND" ~ "ND", # not detected
                                conc == "NM" ~ "NM", # not measured
                                conc == 'na' ~ 'na')) %>%  # not available
  mutate(conc_ppb = parse_number(conc)) %>% 
  mutate(conc_ppb = case_when(nondetects == "ND" ~ 0, TRUE ~ conc_ppb)) %>%
  droplevels() %>% 
  mutate(pfas_type = case_when(pfas_type == "X6.2FTS" ~ "6:2FTS",
                               pfas_type == "X8.2FTS" ~ "8:2FTS",
                               TRUE ~ pfas_type)) %>% 
  filter(!grepl("lank", sample_id)) %>% 
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
  rename(sample_id = field_id,
         pfas_type = compound) %>% 
  mutate(nondetects = case_when(conc == "ND" ~ "ND", # not detected
                                conc == "NM" ~ "NM", # not measured
                                conc == 'na' ~ 'na')) %>%  # not available
  mutate(conc_ppb = parse_number(conc)) %>% 
  mutate(conc_ppb = case_when(nondetects == "ND" ~ 0, TRUE ~ conc_ppb)) %>% 
  droplevels() %>% 
  mutate(pfas_type = case_when(pfas_type == "X6.2FTS" ~ "6:2FTS",
                               pfas_type == "X8.2FTS" ~ "8:2FTS",
                               TRUE ~ pfas_type)) %>% 
  filter(!grepl("lank", sample_id)) %>% 
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
  mutate(nondetects = case_when(conc == "ND" ~ "ND", # not detected
                                conc == "NM" ~ "NM", # not measured
                                conc == 'na' ~ 'na')) %>%  # not available
  mutate(conc_ppb = parse_number(conc)) %>% 
  mutate(conc_ppb = case_when(nondetects == "ND" ~ 0, TRUE ~ conc_ppb)) %>% 
  droplevels() %>% 
  mutate(pfas_type = case_when(pfas_type == "X6.2FTS" ~ "6:2FTS",
                               pfas_type == "X8.2FTS" ~ "8:2FTS",
                               TRUE ~ pfas_type)) %>% 
  filter(!grepl("lank", sample_id)) %>% 
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

# insect mass -------------------------------------------------------------
insect_mass = read_csv("data/Farmington_InsectMasses.csv") %>% 
  mutate(units = "grams") %>% 
  mutate(gdw = individual_mass_ww) 

saveRDS(insect_mass, file = "data/insect_mass.rds")

# isotopes ----------------------------------------------------------------

isotopes = read_csv("data/Farmington_PFAS_Stable_Isotopes_Datarelease.csv") %>% 
  mutate(original_d15n = d15n) %>% 
  mutate(lifestage = case_when(family != "Spider" ~ lifestage)) %>% 
  # mutate(d15n = case_when(lifestage == "Adult" ~ d15n - 1, TRUE ~ d15n)) %>% # correct adult enrichment due to metamorphosis (Kraus et al. 2014 EST)
  group_by(site) %>% 
  mutate(d15n_s = scale(d15n),
         mean_15n = mean(d15n, na.rm = T),
         mean_13c = mean(d13c, na.rm = T),
         sd_15n = sd(d15n, na.rm = T),
         mean_centered_15n = d15n - mean_15n,
         mean_centered_13c = d13c - mean_13c) %>% 
  ungroup %>% 
  rename(taxon = family,
         type = lifestage) %>% 
  mutate(taxon = case_when(is.na(taxon) ~ media, TRUE ~ taxon),            # make the names match with pfas data for combining later
         type = case_when(type == "Adult" ~ "Emergent", TRUE ~ type)) %>% 
  mutate(type_taxon = paste0(type, "_", taxon),
         type_taxon = str_remove(type_taxon, "NA_"),
         type_taxon = case_when(type_taxon == "Emergent_Spider" ~ "Tetragnathidae", TRUE ~ type_taxon))

saveRDS(isotopes, file = "data/isotopes.rds")
