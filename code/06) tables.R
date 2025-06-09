
pfas_names = read_csv("data/log_kw.csv") %>% 
  mutate(pfas_type = str_remove(pfas_type, " ")) %>% 
  mutate(log_kmw_mean = parse_number(str_sub(log_kmw, 1, 4)),
         log_kmw_sd = parse_number(str_sub(log_kmw, -4, -1)),
         log_kpw_mean = parse_number(log_kpw),
         number_carbons = parse_number(no_carbons)) %>% 
  rename(pfas_category = type) %>% 
  distinct(pfas_category, pfas_type, number_carbons) %>% 
  filter(pfas_type %in% unique(merged_d2$pfas_type)) %>% 
  arrange(pfas_category, number_carbons) %>% 
  mutate(pfas_order = 1:nrow(.))

posts_taxon = readRDS(file = "posteriors/posts_taxon.rds")

raw_baf = as_tibble(mod_dat) %>% 
  select(-conc_ppb_s) %>%
  group_by(type, pfas_type, site) %>% 
  reframe(conc_ppb = mean(conc_ppb, na.rm = T)) %>%
  pivot_wider(names_from = type, values_from = conc_ppb) %>% 
  clean_names()  %>%
  group_by(pfas_type, site) %>% 
  mutate(a_kd_water_to_sediment_ttf = sediment/water,
         b_kd_water_to_biofilm_ttf = biofilm/water,
         a2_kd_water_to_seston_ttf = seston/water,
         b2_kd_water_to_detritus_ttf = detritus/water,
         c_baf_sediment_to_larvae_ttf = larval/sediment,
         d_baf_biofilm_to_larvae_ttf = larval/biofilm,
         e_baf_detritus_to_larvae_ttf = larval/detritus,
         f_baf_seston_to_larvae_ttf = larval/seston,
         g_mtf_larvae_to_emergent_ttf = emergent/larval,
         h_ttf_emergent_to_spider_ttf = tetragnathidae/emergent) %>% 
  pivot_longer(cols = ends_with("ttf"))  %>%
  mutate(across(where(is.numeric), ~ na_if(., 0))) %>%
  mutate(across(where(is.numeric), ~ na_if(., Inf))) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))


raw_data_ids = raw_baf %>% 
  group_by(pfas_type, name, site) %>% 
  mutate(alpha = case_when(is.na(value) ~ "no raw data",
                           TRUE ~ "raw data")) %>% 
  distinct(pfas_type, name, alpha, site)

posts_taxon_type = posts_taxon %>% 
  group_by(type, site, pfas_type, order, .draw) %>% 
  reframe(.epred = median(.epred)) %>% 
  left_join(pfas_names) %>% 
  mutate(pfas_type = reorder(pfas_type, pfas_order)) %>% 
  left_join(pfas_names)

posts_ttf = posts_taxon_type %>% 
  ungroup %>% 
  select(type, site, pfas_type, .draw, .epred, pfas_category) %>% 
  pivot_wider(names_from = type, values_from = .epred) %>% 
  clean_names() %>% 
  mutate(a_kd_water_to_sediment_ttf = sediment/water,
         b_kd_water_to_biofilm_ttf = biofilm/water,
         a2_kd_water_to_seston_ttf = seston/water,
         b2_kd_water_to_detritus_ttf = detritus/water,
         c_baf_sediment_to_larvae_ttf = larval/sediment,
         d_baf_biofilm_to_larvae_ttf = larval/biofilm,
         e_baf_detritus_to_larvae_ttf = larval/detritus,
         f_baf_seston_to_larvae_ttf = larval/seston,
         g_mtf_larvae_to_emergent_ttf = emergent/larval,
         h_ttf_emergent_to_spider_ttf = tetragnathidae/emergent) %>% 
  pivot_longer(cols = ends_with("ttf")) 


# kd tables ---------------------------------------------------------------

kd_water_median_ci = posts_ttf %>% 
  filter(grepl("kd", name)) %>% 
  group_by(pfas_category, pfas_type, name) %>% 
  # reframe(mean = median(value),
          # sd = sd(value)) %>% 
  median_qi(value) %>% 
  filter(pfas_type %in% c("PFHxA", "PFHpA", "PFOA", "PFNA", "PFUnA", "PFOS", "8:2FTS")) %>% 
  mutate(name_short = case_when(grepl("sediment", name) ~ "KdSed",
                                grepl("seston", name) ~ "KdSes",
                                grepl("detritus", name) ~ "KdDet",
                                grepl("biofilm", name) ~ "KdBio"))  %>% 
  mutate(value = round(value, 1),
         .lower = comma(round(.lower, 0)),
         .upper = comma(round(.upper, 0))) %>%
  mutate(median_ci = paste0(value, " (",.lower, " to ", .upper, ")")) %>% 
  select(pfas_category, pfas_type, name_short, median_ci) %>% 
  pivot_wider(names_from = name_short, values_from = median_ci) %>% 
  select(pfas_category, pfas_type, KdSed, everything()) %>% 
  mutate(pfas_type = fct_relevel(pfas_type, c("PFHxA", "PFHpA", "PFOA", "PFNA", "PFUnA", "PFOS", "8:2FTS"))) %>% 
  arrange(pfas_type)

kd_water_mean_sd = posts_ttf %>% 
  filter(grepl("kd", name)) %>% 
  group_by(pfas_category, pfas_type, name) %>% 
  reframe(mean = mean(value),
          sd = sd(value)) %>%
  # median_qi(value) %>% 
  filter(pfas_type %in% c("PFHxA", "PFHpA", "PFOA", "PFNA", "PFUnA", "PFOS", "8:2FTS")) %>% 
  mutate(name_short = case_when(grepl("sediment", name) ~ "KdSed",
                                grepl("seston", name) ~ "KdSes",
                                grepl("detritus", name) ~ "KdDet",
                                grepl("biofilm", name) ~ "KdBio")) %>% 
  mutate(mean_sd = paste0(round(mean,0), " \u00b1 ", round(sd, 0))) %>% 
  select(pfas_category, pfas_type, name_short, mean_sd) %>% 
  pivot_wider(names_from = name_short, values_from = mean_sd) %>% 
  select(pfas_category, pfas_type, KdSed, everything()) %>% 
  mutate(pfas_type = fct_relevel(pfas_type, c("PFHxA", "PFHpA", "PFOA", "PFNA", "PFUnA", "PFOS", "8:2FTS"))) %>% 
  arrange(pfas_type)

write_csv(kd_water_median_ci, file = "tables/kd_water_median_ci.csv")
write_csv(kd_water_mean_sd, file = "tables/kd_water_mean_sd.csv")
# bsaf (i.e., sediment to bug) tables -----------------------------------------------------

posts_taxon_type_taxa = posts_taxon %>% 
  # group_by(type, site, pfas_type, order, .draw) %>% 
  # reframe(.epred = median(.epred)) %>% 
  left_join(pfas_names) %>% 
  mutate(pfas_type = reorder(pfas_type, pfas_order)) %>% 
  left_join(pfas_names)

posts_baf_taxa = posts_taxon_type_taxa %>% 
  ungroup %>% 
  select(type_taxon, site, pfas_type, .draw, .epred, pfas_category) %>% 
  pivot_wider(names_from = type_taxon, values_from = .epred) %>% 
  clean_names() %>% 
  mutate(zbsaf_dip = larval_diptera/sediment,
         zbsaf_eph = larval_ephemeroptera/sediment,
         zbsaf_meg = larval_megaloptera/sediment,
         zbsaf_odo = larval_odonata/sediment,
         zbsaf_ple = larval_plecoptera/sediment,
         zbsaf_tri = larval_trichoptera/sediment,
         zbaf_dip = larval_diptera/water,
         zbaf_eph = larval_ephemeroptera/water,
         zbaf_meg = larval_megaloptera/water,
         zbaf_odo = larval_odonata/water,
         zbaf_ple = larval_plecoptera/water,
         zbaf_tri = larval_trichoptera/water,
         zbafselect_eph = larval_ephemeroptera/biofilm,
         zbafselect_ple = larval_plecoptera/detritus,
         zbafselect_tri = larval_trichoptera/seston) %>% 
  pivot_longer(cols = starts_with("z")) %>% 
  group_by(pfas_type, draw, pfas_category, name) %>% 
  reframe(value = median(value, na.rm = T))

name_pfas_order = tibble(name_pfas = c("bsafPFHxA","bsafPFHpA","bsafPFOA","bsafPFNA", "bsafPFOS", "bsaf6:2FTS",
                                       "bafPFHxA","bafPFHpA","bafPFOA","bafPFNA","bafPFDA","bafPFUnA","bafPFHxS", "bafPFOS", "baf6:2FTS", "baf8:2FTS",
                                       "bafselectPFHxA","bafselectPFHpA","bafselectPFOA","bafselectPFNA", "bafselectPFUnA","bafselectPFHxS", "bafselectPFOS", "bafselect6:2FTS")) %>% 
  mutate(order = 1:nrow(.))

baf_median_ci = posts_baf_taxa %>% 
  select(pfas_type, pfas_category, draw, name, value) %>% 
  group_by(pfas_category, pfas_type, name) %>% 
  median_qi(value, na.rm = T) %>% 
  separate(name, into = c("measure", "taxon"), remove = F) %>% 
  mutate(measure = str_remove(measure, "z")) %>% 
  mutate(name_pfas = paste0(measure, pfas_type)) %>% 
  filter(name_pfas %in% c("bsafPFHxA","bsafPFHpA","bsafPFOA","bsafPFNA", "bsafPFOS", "bsaf6:2FTS",
                          "bafPFHxA","bafPFHpA","bafPFOA","bafPFNA","bafPFDA","bafPFUnA","bafPFHxS", "bafPFOS", "baf6:2FTS", "baf8:2FTS",
                          "bafselectPFHxA","bafselectPFHpA","bafselectPFOA","bafselectPFNA", "bafselectPFUnA","bafselectPFHxS", "bafselectPFOS", "bafselect6:2FTS")) %>% 
  left_join(name_pfas_order) %>% 
  arrange(order) %>% 
  # filter(!is.na(value))  %>% 
  mutate(value = round(value, 1),
         .lower = comma(round(.lower, 0)),
         .upper = comma(round(.upper, 0))) %>%
  mutate(median_ci = paste0(value, " (",.lower, " to ", .upper, ")")) %>% 
  select(pfas_category, measure, pfas_type, taxon, median_ci) %>% 
  pivot_wider(names_from = taxon, values_from = median_ci)

baf_mean_sd = posts_baf_taxa %>% 
  select(pfas_type, pfas_category, draw, name, value) %>% 
  group_by(pfas_category, pfas_type, name) %>% 
  reframe(mean = mean(value, na.rm = T),
          sd = sd(value, na.rm = T)) %>% 
  separate(name, into = c("measure", "taxon"), remove = F) %>% 
  mutate(measure = str_remove(measure, "z")) %>% 
  mutate(name_pfas = paste0(measure, pfas_type)) %>% 
  filter(name_pfas %in% c("bsafPFHxA","bsafPFHpA","bsafPFOA","bsafPFNA", "bsafPFOS", "bsaf6:2FTS",
                          "bafPFHxA","bafPFHpA","bafPFOA","bafPFNA","bafPFDA","bafPFUnA","bafPFHxS", "bafPFOS", "baf6:2FTS", "baf8:2FTS",
                          "bafselectPFHxA","bafselectPFHpA","bafselectPFOA","bafselectPFNA", "bafselectPFUnA","bafselectPFHxS", "bafselectPFOS", "bafselect6:2FTS")) %>% 
  left_join(name_pfas_order) %>% 
  arrange(order) %>% 
  # filter(!is.na(value)) %>% 
  mutate(mean_sd = paste0(round(mean,0), " \u00b1 ", round(sd, 0))) %>% 
  select(pfas_category, measure, pfas_type, taxon, mean_sd) %>% 
  pivot_wider(names_from = taxon, values_from = mean_sd)

write_csv(baf_median_ci, file = "tables/baf_median_ci.csv")
write_csv(baf_mean_sd, file = "tables/baf_mean_sd.csv")

# metamorphic retention -----------------------------------------------------

posts_taxon_type_taxa = posts_taxon %>% 
  # group_by(type, site, pfas_type, order, .draw) %>% 
  # reframe(.epred = median(.epred)) %>% 
  left_join(pfas_names) %>% 
  mutate(pfas_type = reorder(pfas_type, pfas_order)) %>% 
  left_join(pfas_names)

posts_metamorphic_taxa = posts_taxon_type_taxa %>% 
  filter(!is.na(taxon)) %>%
  ungroup %>% 
  select(type, taxon, site, pfas_type, .draw, .epred, pfas_category) %>% 
  pivot_wider(names_from = type, values_from = .epred) %>% 
  clean_names() %>% 
  mutate(mtf = emergent/larval) %>% 
  group_by(pfas_type, draw, pfas_category, taxon) %>% 
  reframe(value = median(mtf, na.rm = T))

mtf_to_keep = tibble(pfas_type = c("PFHxA", "PFHpA", "PFOA", "PFNA", "PFUnA", "PFDA", "PFDoA", "PFHxS",
                        "PFOS", "6:2FTS")) %>% 
  mutate(order = 1:nrow(.))

# overall average (all emerged/larval, averaged over sites)
posts_metamorphic_overall = posts_taxon_type_taxa %>% 
  filter(type == "Emergent" | type == "Larval") %>% 
  ungroup %>% 
  select(type, taxon, site, pfas_type, .draw, .epred, pfas_category) %>% 
  pivot_wider(names_from = type, values_from = .epred) %>% 
  clean_names() %>% 
  mutate(mtf = emergent/larval) %>% 
  group_by(pfas_type, pfas_category) %>%
  mutate(mean = mean(mtf, na.rm = T),
         sd = sd(mtf, na.rm = T)) %>% 
  group_by(pfas_category, pfas_type, mean, sd) %>% 
  median_qi(mtf, na.rm = T) %>% 
  mutate(overall = paste0(round(mtf, 1), " (", round(.lower,0), " to ", round(.upper, 0), ")"))  %>% 
  mutate(mean_sd = paste0(round(mean,0), " \u00b1 ", round(sd, 0)))

mtf_median_ci = posts_metamorphic_taxa %>% 
  group_by(pfas_category, pfas_type, taxon) %>% 
  median_qi(value) %>% 
  mutate(value = round(value, 1),
         .lower = comma(round(.lower, 0)),
         .upper = comma(round(.upper, 0))) %>%
  mutate(median_ci = paste0(value, " (",.lower, " to ", .upper, ")")) %>% 
  select(pfas_category, pfas_type, taxon, median_ci) %>% 
  pivot_wider(names_from = taxon, values_from = median_ci) %>% 
  left_join(mtf_to_keep) %>% 
  arrange(order) %>% 
  filter(pfas_type %in% c(unique(mtf_to_keep$pfas_type))) %>% 
  left_join(posts_metamorphic_overall %>% select(pfas_type, overall))

mtf_mean_sd = posts_metamorphic_taxa %>% 
  group_by(pfas_category, pfas_type, taxon) %>% 
  reframe(mean = mean(value, na.rm = T),
          sd = sd(value, na.rm = T)) %>% 
  mutate(mean_sd = paste0(round(mean,0), " \u00b1 ", round(sd, 0))) %>% 
  select(pfas_category, pfas_type, taxon, mean_sd)  %>% 
  pivot_wider(names_from = taxon, values_from = mean_sd) %>% 
  left_join(mtf_to_keep) %>% 
  arrange(order) %>% 
  filter(pfas_type %in% c(unique(mtf_to_keep$pfas_type)))

write_csv(mtf_median_ci, file = "tables/mtf_median_ci.csv")
write_csv(mtf_mean_sd, file = "tables/mtf_mean_sd.csv")


# trophic transfer --------------------------------------------------------

posts_taxon_type_taxa = posts_taxon %>% 
  # group_by(type, site, pfas_type, order, .draw) %>% 
  # reframe(.epred = median(.epred)) %>% 
  left_join(pfas_names) %>% 
  mutate(pfas_type = reorder(pfas_type, pfas_order)) %>% 
  left_join(pfas_names)

posts_ttf_taxa = posts_taxon_type_taxa %>% 
  mutate(taxon = case_when(type == "Tetragnathidae" ~ "Tetragnathidae",
                           TRUE ~ taxon)) %>%
  filter(!is.na(taxon)) %>%
  filter(type == "Emergent" | type == "Tetragnathidae") %>%
  ungroup %>% 
  select(taxon, site, pfas_type, .draw, .epred, pfas_category) %>% 
  pivot_wider(names_from = taxon, values_from = .epred) %>% 
  clean_names() %>% 
  pivot_longer(cols = c(diptera, ephemeroptera, odonata, trichoptera, plecoptera)) %>% 
  mutate(ttf = tetragnathidae/value) %>% 
  group_by(pfas_type, draw, pfas_category, name) %>% 
  reframe(value = median(ttf, na.rm = T))

ttf_to_keep = tibble(pfas_type = c("PFHxA", "PFHpA", "PFOA", "PFNA", "PFUnA", "PFDA", "PFDoA", "PFHxS",
                                   "PFOS", "6:2FTS")) %>% 
  mutate(order = 1:nrow(.))

ttf_median_ci = posts_ttf_taxa %>% 
  group_by(pfas_category, pfas_type, name) %>% 
  median_qi(value)  %>% 
  mutate(value = round(value, 1),
         .lower = comma(round(.lower, 0)),
         .upper = comma(round(.upper, 0))) %>%
  mutate(median_ci = paste0(value, " (",.lower, " to ", .upper, ")")) %>% 
  select(pfas_category, pfas_type, name, median_ci) %>% 
  pivot_wider(names_from = name, values_from = median_ci) %>% 
  left_join(ttf_to_keep) %>% 
  arrange(order) %>% 
  filter(pfas_type %in% c(unique(ttf_to_keep$pfas_type)))

ttf_mean_sd = posts_ttf_taxa %>% 
  group_by(pfas_category, pfas_type, name) %>% 
  reframe(mean = mean(value, na.rm = T),
          sd = sd(value, na.rm = T)) %>% 
  mutate(mean_sd = paste0(round(mean,0), " \u00b1 ", round(sd, 0))) %>% 
  select(pfas_category, pfas_type, name, mean_sd)  %>% 
  pivot_wider(names_from = name, values_from = mean_sd) %>% 
  left_join(ttf_to_keep) %>% 
  arrange(order) %>% 
  filter(pfas_type %in% c(unique(ttf_to_keep$pfas_type)))

write_csv(ttf_median_ci, file = "tables/ttf_median_ci.csv")
write_csv(ttf_mean_sd, file = "tables/ttf_mean_sd.csv")


# ppbs --------------------------------------------------------------------

posts_taxon_type = posts_taxon %>% 
  group_by(type, site, pfas_type, order, .draw) %>% 
  reframe(.epred = median(.epred)) %>% 
  left_join(pfas_names) %>% 
  mutate(pfas_type = reorder(pfas_type, pfas_order))

posts_taxa_only = posts_taxon %>% 
  filter(!is.na(taxon)) 

posts_taxon_type_summary = posts_taxon_type %>% 
  group_by(type, pfas_type, pfas_category, order, number_carbons) %>% 
  median_qi(.epred) 

ppb_values = posts_taxon_type_summary %>% 
  mutate(.epred = case_when(type == "Water" ~ round(.epred, 4),
                            TRUE ~ round(.epred, 3)),
         .lower = case_when(type == "Water" ~ round(.lower, 4),
                            TRUE ~ round(.lower, 3)),
         .upper = case_when(type == "Water" ~ round(.upper, 4),
                            TRUE ~ round(.upper, 1))) %>%
  mutate(median_ci = paste0(.epred, " (",.lower, "-", .upper, ")")) %>%
  mutate(type = fct_relevel(type, "Water", "Sediment", "Detritus", "Seston", "Biofilm", "Larval", "Emergent",
                            "Tetragnathidae")) %>% 
  select(type, pfas_type, median_ci) %>%
  pivot_wider(names_from = type, values_from = median_ci)

write_csv(ppb_values, file = "tables/ppb_values.csv")  


ppb_average_per_taxon = posts_taxa_only %>% 
  group_by(type, pfas_type, taxon) %>% 
  median_qi(.epred) %>% 
  mutate(.epred = case_when(type == "Water" ~ round(.epred, 4),
                            TRUE ~ round(.epred, 3)),
         .lower = case_when(type == "Water" ~ round(.lower, 4),
                            TRUE ~ round(.lower, 3)),
         .upper = case_when(type == "Water" ~ round(.upper, 4),
                            TRUE ~ round(.upper, 1))) %>%
  mutate(median_ci = paste0(.epred, " (",.lower, "-", .upper, ")")) %>%
  mutate(type = fct_relevel(type, "Larval", "Emergent")) %>% 
  select(taxon, type, pfas_type, median_ci) %>%
  pivot_wider(names_from = type, values_from = median_ci)

write_csv(ppb_average_per_taxon, file = "tables/ppb_values_pertaxon.csv")  

ppb_values_site_taxa = posts_taxon %>% 
  group_by(type, pfas_type,  order, site, taxon) %>% 
  median_qi(.epred) %>% 
  mutate(.epred = case_when(type == "Water" ~ round(.epred, 4),
                            TRUE ~ round(.epred, 3)),
         .lower = case_when(type == "Water" ~ round(.lower, 4),
                            TRUE ~ round(.lower, 3)),
         .upper = case_when(type == "Water" ~ round(.upper, 4),
                            TRUE ~ round(.upper, 1)))  %>% 
  rename(median = .epred,
         low95 = .lower,
         high95 = .upper) %>% 
  select(site, type, pfas_type, taxon, median, low95, high95)

write_csv(ppb_values_site_taxa , file = "tables/ppb_values_site_taxa .csv")  



# sum pfas by type --------------------------------------------------------

mod1 = readRDS(file = "models/mod1.rds")
mod1_taxa = readRDS(file = "models/mod1_taxa.rds")

posts_taxa_sumpfas = mod1_taxa$data %>% 
  distinct(type, taxon) %>%
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8)) %>% 
  add_epred_draws(mod1_taxa, re_formula = ~ (1 + type|taxon)) %>% 
  mutate(.epred = .epred ,
         .epred = .epred*unique(mod1_taxa$data2$raw_data$mean_sum_ppb))

posts_sumpfas = mod1$data %>% 
  distinct(type) %>%
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8)) %>% 
  add_epred_draws(mod1, re_formula = NA) %>% 
  mutate(.epred = .epred ,
         .epred = .epred*unique(mod1$data2$raw_data$mean_sum_ppb)) 

posts_sumpfas_summary = posts_sumpfas %>% 
  group_by(type, order) %>% 
  median_qi(.epred) 

sum_pfas_type = posts_sumpfas_summary %>%
  mutate(.epred = round(.epred, 2),
         .lower = round(.lower, 2),
         .upper = round(.upper, 2)) %>% 
  mutate(median_cri = paste0(.epred, " (",.lower, " to ", .upper, ")")) %>% 
  ungroup %>% 
  select(type, median_cri) %>% 
  mutate(units = "ppb") %>% 
  mutate(description = "Median and CrI of sum pfas ppb.")

write.csv(sum_pfas_type, file = "tables/sum_pfas_type.csv")

# sum pfas by order and lifestage -----------------------------------------

mod1 = readRDS(file = "models/mod1.rds")
mod1_taxa = readRDS(file = "models/mod1_taxa.rds")

posts_taxa_sumpfas = mod1_taxa$data %>% 
  distinct(type, taxon) %>%
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8)) %>% 
  add_epred_draws(mod1_taxa, re_formula = ~ (1 + type|taxon)) %>% 
  mutate(.epred = .epred ,
         .epred = .epred*unique(mod1_taxa$data2$raw_data$mean_sum_ppb))

posts_taxa_sumpfas_summary = posts_taxa_sumpfas %>% 
  group_by(type, order, taxon) %>% 
  median_qi(.epred) 

sum_pfas_order_stage = posts_taxa_sumpfas_summary %>% 
  group_by(type, taxon) %>% 
  median_qi(.epred) %>%
  mutate(.epred = round(.epred, 2),
         .lower = round(.lower, 2),
         .upper = round(.upper, 2)) %>% 
  mutate(median_cri = paste0(.epred, " (",.lower, " to ", .upper, ")")) %>% 
  ungroup %>% 
  select(type, taxon, median_cri) %>% 
  mutate(units = "ppb") %>% 
  mutate(description = "Median and CrI of sum pfas ppb per taxon and life stage. Sums are first calculated per site, taxon, and life stage. Then those are averaged over sites")

write.csv(sum_pfas_order_stage, file = "tables/sum_pfas_order_stage.csv")


# partitioning coefficients sum pfas --------------------------------------

mod1 = readRDS(file = "models/mod1.rds")
mod1_taxa = readRDS(file = "models/mod1_taxa.rds")

posts_sumpfas = mod1$data %>% 
  distinct(type) %>%
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8)) %>% 
  add_epred_draws(mod1, re_formula = NA) %>% 
  mutate(.epred = .epred ,
         .epred = .epred*unique(mod1$data2$raw_data$mean_sum_ppb)) 

sum_partitioning = posts_sumpfas %>% 
  select(-order) %>% 
  filter(.draw <= 500) %>% 
  ungroup %>% 
  select(-.row, -.chain, -.iteration, -order) %>% 
  group_by(type, .draw) %>% 
  pivot_wider(names_from = type, values_from = .epred) %>% 
  clean_names() %>% 
  mutate(a_kd_water_to_sediment_ttf = sediment/water,
         b_kd_water_to_biofilm_ttf = biofilm/water,
         a2_kd_water_to_seston_ttf = seston/water,
         b2_kd_water_to_detritus_ttf = detritus/water,
         c_baf_sediment_to_larvae_ttf = larval/sediment,
         d_baf_biofilm_to_larvae_ttf = larval/biofilm,
         e_baf_detritus_to_larvae_ttf = larval/detritus,
         f_baf_seston_to_larvae_ttf = larval/seston,
         g_mtf_larvae_to_emergent_ttf = emergent/larval,
         h_ttf_emergent_to_spider_ttf = tetragnathidae/emergent) %>% 
  pivot_longer(cols = ends_with("ttf"), values_to = ".epred") %>% 
    select(draw, name, .epred)  %>% 
  group_by(draw, name) %>% 
  reframe(.epred = median(.epred))  %>% 
  group_by(name) %>% 
  median_qi(.epred) %>%
  mutate(.epred = round(.epred, 2),
         .lower = round(.lower, 2),
         .upper = round(.upper, 2)) %>% 
  mutate(median_cri = paste0(.epred, " (",.lower, " to ", .upper, ")")) %>% 
  ungroup %>% 
  select(name, median_cri) %>% 
  mutate(units = "partitioning coefficients") %>% 
  mutate(description = "Median and CrI of partitioning sum pfas")

write.csv(sum_partitioning, file = "tables/sum_partitioning.csv")

