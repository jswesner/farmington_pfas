library(tidyverse)
library(brms)

hg4_taxon = readRDS(file = "models/hg4_taxon.rds")

pp_check(hg4_taxon) + scale_x_log10()

test = plot(conditional_effects(hg4_taxon), points = T)

test$type_taxon + scale_y_log10()

test$data %>% 
  select(-conc_ppb_s) %>% 
  distinct()
