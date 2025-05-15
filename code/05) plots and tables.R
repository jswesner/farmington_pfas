library(tidybayes)
library(tidyverse)
library(brms)
library(janitor)
library(viridis)
library(ggh4x)
library(scales)
theme_set(theme_default())

scale_fill_custom <- function() {
  scale_fill_manual(values = c(
    "PFHxA" = "#C7E9C0", "PFHpA" = "#A1D99B", "PFOA" = "#41AB5D", "PFNA" = "#006600",
    "PFDA" = "darkslategray3", "PFUnA" = "darkslategray4", "PFDoA" = "darkslategray",
    "PFBS" = "#C6DBEF", "PFHxS" = "#9ECAE1", "PFOS" = "#3399ff",
    "6:2FTS" = "grey", "8:2FTS" = "black"
  ))
}

scale_color_custom <- function() {
  scale_color_manual(values = c(
    "PFHxA" = "#C7E9C0", "PFHpA" = "#A1D99B", "PFOA" = "#41AB5D", "PFNA" = "#006600",
    "PFDA" = "darkslategray3", "PFUnA" = "darkslategray4", "PFDoA" = "darkslategray",
    "PFBS" = "#C6DBEF", "PFHxS" = "#9ECAE1", "PFOS" = "#3399ff",
    "6:2FTS" = "grey", "8:2FTS" = "black"
  ))
}

pfas_orders = tibble(pfas_type = c(
  "PFHxA", "PFHpA" , "PFOA" , "PFNA" ,
  "PFDA", "PFUnA", "PFDoA" ,
  "PFBS" , "PFHxS" , "PFOS" ,
  "6:2FTS" , "8:2FTS")) %>% 
  mutate(order = 1:nrow(.))

# load model
hg4_taxon = readRDS(file = "models/hg4_taxon.rds")
merged_d2 = hg4_taxon$data2$merged_d2

max_conc = unique(merged_d2$max_conc)

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

mod_dat = hg4_taxon$data %>% 
  separate(type_taxon, into = c("type", "taxon"), remove = F) %>% 
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8)) %>% 
  mutate(conc_ppb = conc_ppb_s*max_conc) %>%
  left_join(pfas_names) %>% 
  mutate(pfas_type = reorder(pfas_type, pfas_order))

posts_taxon = mod_dat  %>%
  select(-contains("conc_ppb")) %>%
  distinct() %>%
  add_epred_draws(hg4_taxon, re_formula = NULL, dpar = T, ndraws = 500) %>%
  mutate(.epred = .epred*max_conc) 

# saveRDS(posts_taxon, file = "posteriors/posts_taxon.rds")
# posts_taxon = readRDS(file = "posteriors/posts_taxon.rds")

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


pfas_concentration = mod_dat %>% 
  ggplot(aes(x = reorder(type, order), y = conc_ppb + 1)) + 
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, shape = 1,
              size = 0.3) +
  geom_pointrange(data = posts_taxon_type_summary, aes(y = .epred + 1, 
                                                       ymin = .lower + 1,
                                                       ymax = .upper + 1,
                                                       color = pfas_type), 
                  size = 0.3) + 
  scale_y_log10(breaks = c(1, 10, 100), 
                labels = c("0", "10", "100")) +
  # facet_wrap(pfas_category~reorder(pfas_type, number_carbons)) +
  facet_wrap2(~pfas_type) +
  scale_color_custom() +
  labs(y = "PFAS Concentration (ppb)",
       x = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  guides(color = "none") +
  NULL

ggsave(pfas_concentration, file = "plots/pfas_concentration.jpg", width = 6.5, height = 9)

posts_taxon_type_summary_site = posts_taxon_type %>% 
  group_by(type, pfas_type, site, order) %>% 
  median_qi(.epred) %>% 
  left_join(pfas_names) %>% 
  mutate(pfas_type = reorder(pfas_type, pfas_order))


pfas_concentration_bysite = mod_dat %>% 
  ggplot(aes(x = reorder(type, order), y = conc_ppb + 1, fill = site)) + 
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, shape = 20,
              size = 0.3, aes(color = site)) +
  geom_pointrange(data = posts_taxon_type_summary_site, aes(y = .epred + 1, 
                                                            ymin = .lower + 1,
                                                            ymax = .upper + 1), 
                  size = 0.3, shape = 21, color = "black") + 
  scale_y_log10(breaks = c(1, 10, 100), 
                labels = c("0", "10", "100")) +
  facet_wrap(~pfas_type) +
  labs(y = "PFAS Concentration (ppb)",
       x = "",
       color = "") +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  # guides(color = "none") + 
  NULL

ggsave(pfas_concentration_bysite, file = "plots/pfas_concentration_bysite.jpg", width = 8, height = 8)


# get transfer factors ------------------
posts_ttf = posts_taxon_type %>% 
  ungroup %>% 
  select(type, site, pfas_type, .draw, .epred) %>% 
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

saveRDS(posts_ttf, file = "posteriors/posts_ttf.rds")

posts_ttf %>% 
  group_by(pfas_type, name) %>% 
  median_qi(value)

ttf_table_notaxa = posts_ttf %>% 
  group_by(draw, name, pfas_type) %>% 
  reframe(value = median(value)) %>% 
  group_by(name, pfas_type) %>% 
  median_qi(value, .width = 0.75) %>% 
  mutate(value = round(value, 1)) %>% 
  mutate(value_ci = paste0(round(value, 1), " (",
                           round(.lower, 1), " to ",
                           round(.upper, 1), ")")) %>% 
  select(-value, -.lower, -.upper, -.width, -.point, -.interval) %>% 
  pivot_wider(names_from = name, values_from = value_ci) %>% 
  mutate(category = case_when(pfas_type == "6:2FTS" ~ "FTSA",
                              pfas_type %in% c("PFBS", "PFHxS", "PFOS") ~ "PFSAs",
                              T ~ "PFCAs"),
         category = as.factor(category),
         category = fct_relevel(category, "PFCAs", "PFSAs")) %>% 
  select(category, everything()) %>% 
  arrange(category, pfas_type)

write_csv(ttf_table_notaxa, file = "tables/ttf_table_notaxa.csv")

ttf_table_taxa = posts_taxon %>% 
  # filter(site == "Russell Brook") %>%
  # filter(pfas_type == "PFOS") %>%
  filter(.draw <= 500) %>%
  ungroup %>% 
  select(type, site, pfas_type, .draw, .epred, taxon) %>% 
  mutate(type = paste0(type, "_", taxon)) %>% 
  select(-taxon) %>% 
  mutate(type = str_remove(type, "_NA"),
         type = str_to_lower(type)) %>% 
  pivot_wider(names_from = type, values_from = .epred) %>% 
  pivot_longer(cols = c(starts_with("larval"))) %>%
  mutate(c_baf_sediment_ttf = value/sediment,
         d_baf_biofilm_ttf = value/biofilm,
         e_baf_detritus_ttf = value/detritus,
         f_baf_seston_ttf = value/seston,
         g_mtf_ttf = emergent_trichoptera/value,
         h_ttf_ttf = tetragnathidae/value
  ) %>% 
  pivot_longer(cols = ends_with("ttf"), names_to = "transfer_factor",
               values_to = "transfer_value") %>% 
  select(transfer_factor, transfer_value, .draw, site, pfas_type, name) %>% 
  filter(!is.na(transfer_value)) %>% 
  group_by(name, transfer_factor, pfas_type) %>% 
  reframe(ttf = median(transfer_value),
          lower = quantile(transfer_value, probs = 0.025),
          upper = quantile(transfer_value, probs = 0.975)) %>% 
  mutate(ttf = paste0(round(ttf, 1), 
                      " (",
                      round(lower,1),
                      " to ", 
                      round(upper, 1),
                      ")")) %>% 
  select(-lower, -upper) %>% 
  # ungroup %>% distinct(transfer_factor)
  filter(name == "larval_ephemeroptera" & transfer_factor == "d_baf_biofilm_ttf" |
           name == "larval_trichoptera" & transfer_factor == "f_baf_seston_ttf" |
           name == "larval_plecoptera" & transfer_factor == "e_baf_detritus_ttf") %>% 
  mutate(name = str_sub(name, 8, 12),
         transfer_factor = str_sub(transfer_factor, 7, 10)) %>% 
  mutate(baf = paste0("baf_", transfer_factor, "_", name)) %>% 
  select(pfas_type, ttf, baf) %>% 
  pivot_wider(names_from = baf, values_from = ttf) %>% 
  mutate(category = case_when(pfas_type == "6:2FTS" ~ "FTSA",
                              pfas_type %in% c("PFBS", "PFHxS", "PFOS") ~ "PFSAs",
                              T ~ "PFCAs"),
         category = as.factor(category),
         category = fct_relevel(category, "PFCAs", "PFSAs")) %>% 
  select(category, pfas_type, baf_biof_ephem, baf_detr_pleco, baf_sest_trich) %>% 
  arrange(category, pfas_type)

write_csv(ttf_table_taxa, file = "tables/ttf_table_taxa.csv")


# ttf by site

posts_ttf_site = posts_taxon_type %>% 
  ungroup %>% 
  select(type, site, pfas_type, .draw, .epred) %>% 
  pivot_wider(names_from = type, values_from = .epred) %>% 
  clean_names() %>% 
  mutate(bat_to_sed = sediment/water,
         wat_to_biof = biofilm/water,
         wat_to_sest = seston/water,
         wat_to_det = detritus/water,
         sed_to_larv = larval/sediment,
         biof_to_larv = larval/biofilm,
         det_to_larv = larval/detritus,
         sest_to_larv = larval/seston,
         larv_to_emerge = emergent/larval,
         emerg_to_spid = tetragnathidae/emergent) %>% 
  select(site, pfas_type, draw, contains("_to_")) %>% 
  pivot_longer(cols = c(-site, -pfas_type, -draw)) 

ttf_table_notaxa_site = posts_ttf %>% 
  group_by(draw, name, pfas_type, site) %>% 
  reframe(value = median(value)) %>% 
  group_by(name, pfas_type, site) %>% 
  median_qi(value, .width = 0.95) %>% 
  mutate(value = round(value, 1)) %>% 
  mutate(value_ci = paste0(round(value, 1), " (",
                           round(.lower, 1), " to ",
                           round(.upper, 1), ")")) %>% 
  select(-value, -.lower, -.upper, -.width, -.point, -.interval) %>% 
  pivot_wider(names_from = name, values_from = value_ci) %>% 
  mutate(category = case_when(pfas_type == "6:2FTS" ~ "FTSA",
                              pfas_type %in% c("PFBS", "PFHxS", "PFOS") ~ "PFSAs",
                              T ~ "PFCAs"),
         category = as.factor(category),
         category = fct_relevel(category, "PFCAs", "PFSAs")) %>% 
  select(category, everything()) %>% 
  arrange(category, pfas_type)

write_csv(ttf_table_notaxa_site, file = "tables/ttf_table_notaxa_site.csv")


# ttf table taxa by site
ttf_table_taxa_site = posts_taxon %>% 
  # filter(site == "Russell Brook") %>%
  # filter(pfas_type == "PFOS") %>%
  filter(.draw <= 500) %>%
  ungroup %>% 
  select(type, site, pfas_type, .draw, .epred, taxon) %>% 
  mutate(type = paste0(type, "_", taxon)) %>% 
  select(-taxon) %>% 
  mutate(type = str_remove(type, "_NA"),
         type = str_to_lower(type)) %>% 
  pivot_wider(names_from = type, values_from = .epred) %>% 
  pivot_longer(cols = c(starts_with("larval"))) %>%
  mutate(c_baf_sediment_ttf = value/sediment,
         d_baf_biofilm_ttf = value/biofilm,
         e_baf_detritus_ttf = value/detritus,
         f_baf_seston_ttf = value/seston,
         # g_mtf_ttf = emergent_trichoptera/value,
         # h_ttf_ttf = tetragnathidae/value
  ) %>% 
  pivot_longer(cols = ends_with("ttf"), names_to = "transfer_factor",
               values_to = "transfer_value") %>% 
  select(transfer_factor, transfer_value, .draw, site, pfas_type, name) %>% 
  filter(!is.na(transfer_value)) %>% 
  group_by(name, transfer_factor, pfas_type, site) %>% 
  reframe(ttf = median(transfer_value),
          lower = quantile(transfer_value, probs = 0.025),
          upper = quantile(transfer_value, probs = 0.975)) %>% 
  mutate(ttf = paste0(round(ttf, 1), 
                      " (",
                      round(lower,1),
                      " to ", 
                      round(upper, 1),
                      ")")) %>% 
  select(-lower, -upper) %>% 
  # ungroup %>% distinct(transfer_factor)
  filter(name == "larval_ephemeroptera" & transfer_factor == "d_baf_biofilm_ttf" |
           name == "larval_trichoptera" & transfer_factor == "f_baf_seston_ttf" |
           name == "larval_plecoptera" & transfer_factor == "e_baf_detritus_ttf") %>% 
  mutate(name = str_sub(name, 8, 12),
         transfer_factor = str_sub(transfer_factor, 7, 10)) %>% 
  mutate(baf = paste0("baf_", transfer_factor, "_", name)) %>% 
  select(pfas_type, ttf, baf, site) %>% 
  pivot_wider(names_from = baf, values_from = ttf) %>% 
  mutate(category = case_when(pfas_type == "6:2FTS" ~ "FTSA",
                              pfas_type %in% c("PFBS", "PFHxS", "PFOS") ~ "PFSAs",
                              T ~ "PFCAs"),
         category = as.factor(category),
         category = fct_relevel(category, "PFCAs", "PFSAs")) %>% 
  select(category, pfas_type, baf_biof_ephem, baf_detr_pleco, baf_sest_trich, site) %>% 
  arrange(category, pfas_type)

write_csv(ttf_table_taxa_site, file = "tables/ttf_table_taxa_site.csv")

ttf_site_all = left_join(ttf_table_notaxa_site, ttf_table_taxa_site) %>% 
  rename(water_sediment = a_kd_sediment_ttf,
         water_seston = a2_kd_seston_ttf,
         water_biofilm = b_kd_biofilm_ttf,
         water_detritus = b2_kd_detritus_ttf,
         bsaf = c_baf_sediment_ttf,
         biofilm_mayfly = baf_biof_ephem,
         seston_caddis = baf_sest_trich,
         detritus_stone = baf_detr_pleco,
         mtf = g_mtf_ttf,
         ttf = h_ttf_ttf) %>% 
  mutate(measure = "Median(95%CrI)") %>% 
  select(category, pfas_type, site, measure,
         starts_with("water"),
         bsaf, biofilm_mayfly, seston_caddis, detritus_stone,
         mtf, ttf)

write_csv(ttf_site_all, file = "tables/ttf_site_all.csv")


# taxon ttf and mtf -------------------------------------------------------

ttf_full_table = posts_taxon %>% 
  # filter(site == "Russell Brook") %>%
  # filter(pfas_type == "PFOS") %>%
  filter(.draw <= 100) %>%
  ungroup %>% 
  select(type, site, pfas_type, .draw, .epred, taxon) %>% 
  mutate(type = paste0(type, "_", taxon)) %>% 
  select(-taxon) %>% 
  mutate(type = str_remove(type, "_NA"),
         type = str_to_lower(type)) %>% 
  pivot_wider(names_from = type, values_from = .epred) %>% 
  mutate(ttf_diptera = emergent_diptera/tetragnathidae,
         ttf_ephem = emergent_ephemeroptera/tetragnathidae,
         ttf_odonata = emergent_odonata/tetragnathidae,
         ttf_plecoptera = emergent_plecoptera/tetragnathidae,
         ttf_trichoptera = emergent_trichoptera/tetragnathidae) %>% 
  pivot_longer(cols = starts_with("ttf")) %>% 
  group_by(.draw, pfas_type, name) %>% 
  reframe(median = median(value, na.rm = T)) %>% 
  group_by(pfas_type, name) %>% 
  reframe(mean = mean(median),
          sd = sd(median))

ttf_full_table_csv = ttf_full_table %>% 
  mutate(across(where(is.numeric), round, 1)) %>% 
  unite("mean_sd", mean:sd, sep = " +/- ") %>% 
  pivot_wider(names_from = name, values_from = mean_sd) %>% 
  mutate(pfas_type = fct_relevel(pfas_type, "PFHxA", "PFHpA", "PFOA",
                                 "PFNA", "PFUnA", "PFDA", "PFDoA", "PFHxS", "PFOS", 
                                 "6:2 FTS")) %>% 
  arrange(pfas_type)

write_csv(ttf_full_table_csv, file = "tables/ttf_full_table_csv.csv")

posts_taxon %>% 
  # filter(site == "Russell Brook") %>%
  # filter(pfas_type == "PFOS") %>%
  filter(.draw <= 500) %>%
  ungroup %>% 
  select(type, site, pfas_type, .draw, .epred, taxon) %>% 
  mutate(type = paste0(type, "_", taxon)) %>% 
  select(-taxon) %>% 
  mutate(type = str_remove(type, "_NA"),
         type = str_to_lower(type)) %>% 
  pivot_wider(names_from = type, values_from = .epred) %>% 
  mutate(ttf_diptera = emergent_diptera/tetragnathidae,
         ttf_ephem = emergent_ephemeroptera/tetragnathidae,
         ttf_odonata = emergent_odonata/tetragnathidae,
         ttf_plecoptera = emergent_plecoptera/tetragnathidae,
         ttf_trichoptera = emergent_trichoptera/tetragnathidae) %>% 
  pivot_longer(cols = starts_with("ttf")) %>% 
  group_by(.draw, pfas_type, name) %>% 
  reframe(median = median(value, na.rm = T)) %>% 
  group_by(pfas_type, name) %>% 
  reframe(mean = mean(median),
          sd = sd(median))


# plot_ttf ----------------------------------------------------------------

raw_ttf = as_tibble(mod_dat) %>% 
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


posts_ttf_summary = posts_ttf %>% 
  group_by(draw, name, pfas_type) %>% 
  reframe(value = median(value)) %>% 
  group_by(name, pfas_type) %>%
  median_qi(value) %>% 
  left_join(pfas_names) %>% 
  mutate(pfas_type = fct_reorder(pfas_type, pfas_order))

raw_nas_for_color = raw_ttf %>% 
  group_by(pfas_type, name) %>% 
  mutate(alpha = case_when(is.na(value) ~ "no raw data (estimate entirely from the priors/hierarcical effects)",
                           TRUE ~ "raw data (estimate from data and model)")) %>% 
  distinct(pfas_type, name, alpha) %>% 
  left_join(pfas_names) %>% 
  mutate(pfas_type = fct_reorder(pfas_type, pfas_order))
# 
# plot_ttf = posts_ttf_summary %>% 
#   left_join(raw_nas_for_color) %>% 
#   ggplot(aes(x = name, y = value, color = pfas_type)) + 
#   # stat_pointinterval() + 
#   geom_pointrange(aes(ymin = .lower, ymax = .upper, alpha = alpha), size = 0.2) +
#   facet_wrap(~pfas_type) +
#   geom_point(data = raw_ttf, color = "black", size = 0.4) +
#   scale_color_custom() +
#   scale_y_log10(label = c("0.01", "0.1", "1", "10", "100", "1,000", "10,000"), breaks = c(0.01, 0.1, 1, 10, 100, 1000, 10000)) + 
#   guides(color = "none") +
#   labs(y = "Trophic Transer Factor or Partitioning Coefficient",
#        x = "Matrix",
#        caption = "Figure X. Posterior distribution of trophic transfer factors (TTF)") +
#   geom_hline(aes(yintercept = 1), linetype = "dotted") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,
#                                    vjust = 1)) +
#   guides(alpha = "none") +
#   NULL

matrix_new = posts_ttf_summary %>% ungroup %>% distinct(name) %>% 
  separate(name, into = c("sort", "coefficient", "source","delete", "recipient", "ttf"), remove = F) %>% 
  mutate(path = paste0(source," --> ", recipient)) %>% 
  mutate(coefficient_name = case_when(coefficient == "baf" ~ "BAF", 
                                      coefficient == "kd" ~ "kd",
                                      coefficient == "mtf" ~ "MTF",
                                      coefficient == "ttf" ~ "TTF"),
         coefficient_name = as.factor(coefficient_name),
         coefficient_name = fct_relevel(coefficient_name, "kd", "BAF", "MTF", "TTF")) %>% 
  arrange(sort) %>% 
  mutate(path_order = 1:nrow(.)) %>% 
  left_join(pfas_names) 

plot_ttf = posts_ttf_summary %>% 
  left_join(raw_nas_for_color) %>% 
  left_join(matrix_new) %>% 
  left_join(pfas_names) %>% 
  mutate(pfas_type = reorder(pfas_type, pfas_order)) %>% 
  ggplot(aes(x = reorder(path, path_order), y = value, color = pfas_type)) + 
  # stat_pointinterval() + 
  geom_pointrange(aes(ymin = .lower, ymax = .upper, alpha = alpha), size = 0.2) +
  facet_wrap(~pfas_type) +
  geom_point(data = raw_ttf %>% 
               left_join(pfas_names) %>% 
               mutate(pfas_type = reorder(pfas_type, pfas_order)) %>% 
               left_join(matrix_new), color = "black", size = 0.4) +
  scale_color_custom() +
  scale_y_log10(label = c("0.01", "0.1", "1", "10", "100", "1,000", "10,000"), breaks = c(0.01, 0.1, 1, 10, 100, 1000, 10000)) + 
  guides(color = "none") +
  labs(y = "Trophic Transer Factor or Partitioning Coefficient",
       x = "Matrix") +
  geom_hline(aes(yintercept = 1), linetype = "dotted") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   vjust = 0.3)) +
  guides(alpha = "none") +
  NULL

ggsave(plot_ttf, file = "plots/plot_ttf.jpg", width = 7, height = 7)

plot_ttf_remove = plot_ttf + scale_alpha_manual(values = c(0, 1))
ggsave(plot_ttf_remove, file = "plots/plot_ttf_remove.jpg", width = 7, height = 7)


# log_kmw -----------------------------------------------------------------
posts_baf = posts_taxon_type %>% 
  ungroup %>% 
  select(type, site, pfas_type, .draw, .epred) %>% 
  pivot_wider(names_from = type, values_from = .epred) %>% 
  clean_names() %>% 
  mutate(
    baf_sediment_to_larvae = larval/sediment,
    baf_biofilm_to_larvae = larval/biofilm,
    baf_detritus_to_larvae = larval/detritus,
    baf_seston_to_larvae = larval/seston,
    baf_adults_to_spiders = tetragnathidae/emergent,
    baf_sediment_to_adults = emergent/sediment,
    baf_biofilm_to_adults = emergent/biofilm,
    baf_detritus_to_adults = emergent/detritus,
    baf_seston_to_adults = emergent/seston,
    baf_adults_to_spiders = tetragnathidae/emergent) %>% 
  pivot_longer(cols = starts_with("baf")) 

raw_baf = as_tibble(mod_dat) %>% 
  select(-conc_ppb_s) %>%
  group_by(type, pfas_type, site) %>% 
  reframe(conc_ppb = mean(conc_ppb, na.rm = T)) %>%
  pivot_wider(names_from = type, values_from = conc_ppb) %>% 
  clean_names()  %>%
  group_by(pfas_type, site) %>% 
  mutate(baf_sediment_to_larvae = larval/sediment,
         baf_biofilm_to_larvae = larval/biofilm,
         baf_detritus_to_larvae = larval/detritus,
         baf_seston_to_larvae = larval/seston,
         baf_adults_to_spiders = tetragnathidae/emergent,
         baf_sediment_to_adults = emergent/sediment,
         baf_biofilm_to_adults = emergent/biofilm,
         baf_detritus_to_adults = emergent/detritus,
         baf_seston_to_adults = emergent/seston,
         baf_adults_to_spiders = tetragnathidae/emergent) %>% 
  pivot_longer(cols = starts_with("baf"))  %>%
  mutate(across(where(is.numeric), ~ na_if(., 0))) %>%
  mutate(across(where(is.numeric), ~ na_if(., Inf))) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))

log_kmw = read_csv("data/log_kw.csv") %>% 
  mutate(pfas_type = str_remove(pfas_type, " ")) %>% 
  mutate(log_kmw_mean = parse_number(str_sub(log_kmw, 1, 4)),
         log_kmw_sd = parse_number(str_sub(log_kmw, -4, -1)),
         log_kpw_mean = parse_number(log_kpw),
         number_carbons = parse_number(no_carbons)) %>% 
  rename(pfas_category = type)

raw_data_ids = raw_baf %>% 
  group_by(pfas_type, name, site) %>% 
  mutate(alpha = case_when(is.na(value) ~ "no raw data",
                           TRUE ~ "raw data")) %>% 
  distinct(pfas_type, name, alpha, site)

saveRDS(raw_data_ids, file = "data/raw_data_ids.rds")

posts_taxon_kmw = posts_baf %>% 
  left_join(log_kmw %>% select(pfas_type, starts_with("log_kmw_"), log_kpw_mean,
                               number_carbons, pfas_category)) %>% 
  # filter(draw <= 10) %>% 
  left_join(raw_data_ids, by = c("site", "pfas_type", "name")) 

saveRDS(posts_taxon_kmw, file = "posteriors/posts_taxon_kmw.rds")


posts_taxon_kmw_summary = posts_taxon_kmw %>%
  group_by(name, draw, log_kmw_mean,
           log_kpw_mean, number_carbons,
           pfas_category, pfas_type, site) %>% 
  reframe(value = mean(value)) %>%
  group_by(name, log_kmw_mean, log_kpw_mean, site, pfas_type, number_carbons, pfas_category) %>% 
  median_qi(value)

posts_taxon_kmw_adults = posts_baf %>%
  left_join(log_kmw %>% select(pfas_type, starts_with("log_kmw_"), log_kpw_mean,
                               number_carbons, pfas_category, site)) %>% 
  # filter(draw <= 10) %>% 
  left_join(raw_data_ids, by = c("site", "pfas_type", "name")) 

posts_taxon_kmw_summary_adults = posts_taxon_kmw_adults %>%
  group_by(name, draw, log_kmw_mean,
           log_kpw_mean, number_carbons, site) %>% 
  reframe(value = mean(value)) %>%
  group_by(name, log_kmw_mean, log_kpw_mean, number_carbons, site) %>% 
  median_qi(value)

baf_logkmw = posts_taxon_kmw_summary %>% 
  ggplot(aes(x = log_kmw_mean, y = value, color = site)) + 
  # geom_pointrange(aes(ymin = .lower, ymax = .upper)) +
  # geom_line() +
  geom_smooth(se = F) +
  geom_point() +
  # geom_line(data = . %>% filter(draw <= 200),
  #           aes(group = interaction(draw, name)),
  #           alpha = 0.1) +
  # ggthemes::scale_color_colorblind() +
  # scale_color_manual(values = c("FTSA" = "black", "PFSAs" = "#3399ff", "PFCAs" = "green3")) +
  facet_wrap(~ name, scales = "free_y") +
  scale_y_log10() +
  # guides(color= "none") +
  labs(y = "Bioaccumulation Factor", 
       x = expression(log*K[MW])) +
  theme(strip.text.x = element_text(size = 8)) +
  NULL

ggsave(baf_logkmw, file = "plots/baf_logkmw.jpg", width = 8.5, height = 9)

baf_logkpw = posts_taxon_kmw_summary %>% 
  ggplot(aes(x = log_kpw_mean, y = value, color = site)) + 
  # geom_pointrange(aes(ymin = .lower, ymax = .upper)) +
  # geom_line() +
  geom_smooth(se = F) +
  geom_point() +
  # geom_line(data = . %>% filter(draw <= 200),
  #           aes(group = interaction(draw, name)),
  #           alpha = 0.1) +
  # ggthemes::scale_color_colorblind() +
  # scale_color_manual(values = c("FTSA" = "black", "PFSAs" = "#3399ff", "PFCAs" = "green3")) +
  facet_wrap(~ name, scales = "free_y") +
  scale_y_log10() +
  # guides(color= "none") +
  labs(y = "Bioaccumulation Factor", 
       x = expression(log*K[MW])) +
  theme(strip.text.x = element_text(size = 8)) +
  NULL

ggsave(baf_logkpw, file = "plots/baf_logkpw.jpg", width = 8.5, height = 9)


# sum pfas -------------------------------------------------------
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

line_posts = posts_sumpfas %>% 
  group_by(type, order) %>% 
  reframe(.epred = median(.epred)) %>% 
  mutate(group = "lines")

plot_sum_pfas_overall = posts_sumpfas_summary %>% 
  ggplot(aes(x = reorder(type, order),
             y = .epred)) +
  geom_pointrange(aes(ymin = .lower, ymax = .upper)) +
  scale_y_log10(labels = comma) +
  labs(y = "\u2211PFAS (ppb)",
       x = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
        legend.position = "top",
        legend.text = element_text(size = 10),
        legend.title = element_blank()) +
  stat_pointinterval(data = posts_taxa_sumpfas, aes(fill = taxon),
                     shape = 21, .width = 0,
                     position = position_jitter(width = 0.1)) +
  geom_line(data = line_posts %>% filter(type %in% c("Water", "Sediment", "Larval", "Emergent", "Tetragnathidae")),
            aes(group = group),
            linetype = "dashed") +
  geom_line(data = line_posts %>% filter(type %in% c("Water", "Detritus", "Larval", "Emergent", "Tetragnathidae")),
            aes(group = group),
            linetype = "dotted") +
  geom_line(data = line_posts %>% filter(type %in% c("Water", "Seston", "Larval", "Emergent", "Tetragnathidae")),
            aes(group = group),
            linetype = "dotdash") +
  geom_line(data = line_posts %>% filter(type %in% c("Water", "Biofilm", "Larval", "Emergent", "Tetragnathidae")),
            aes(group = group)) + 
  scale_fill_viridis_d(begin = 0.3) +
  NULL

ggsave(plot_sum_pfas_overall, file = "plots/plot_sum_pfas.jpg", width = 8, height = 5)


# sum pfas site -------------------------------------------------------
mod1 = readRDS(file = "models/mod1.rds")
mod1_taxa = readRDS(file = "models/mod1_taxa.rds")

posts_taxa_sumpfas_site = mod1_taxa$data %>% 
  distinct(type, taxon, site) %>%
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8)) %>% 
  add_epred_draws(mod1_taxa, re_formula = NULL) %>% 
  mutate(.epred = .epred ,
         .epred = .epred*unique(mod1_taxa$data2$raw_data$mean_sum_ppb))

posts_sumpfas_site = mod1$data %>% 
  distinct(type, site) %>%
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8)) %>% 
  add_epred_draws(mod1, re_formula = NULL) %>% 
  mutate(.epred = .epred ,
         .epred = .epred*unique(mod1$data2$raw_data$mean_sum_ppb)) 

posts_sumpfas_summary_site = posts_sumpfas_site %>% 
  group_by(type, order, site) %>% 
  median_qi(.epred) 

line_posts_site = posts_sumpfas_site %>% 
  group_by(type, order, site) %>% 
  reframe(.epred = median(.epred)) %>% 
  mutate(group = "lines")

plot_sum_pfas_overall_site = posts_sumpfas_summary_site %>% 
  ggplot(aes(x = reorder(type, order),
             y = .epred)) +
  geom_pointrange(aes(ymin = .lower, ymax = .upper)) +
  scale_y_log10(labels = comma) +
  labs(y = "\u2211PFAS (ppb)",
       x = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
        legend.position = "top",
        legend.text = element_text(size = 10),
        legend.title = element_blank()) +
  stat_pointinterval(data = posts_taxa_sumpfas_site, aes(fill = taxon),
                     shape = 21, .width = 0,
                     size = 0.7,
                     position = position_jitter(width = 0.1)) +
  geom_line(data = line_posts_site %>% filter(type %in% c("Water", "Sediment", "Larval", "Emergent", "Tetragnathidae")),
            aes(group = group),
            linetype = "dashed") +
  geom_line(data = line_posts_site %>% filter(type %in% c("Water", "Detritus", "Larval", "Emergent", "Tetragnathidae")),
            aes(group = group),
            linetype = "dotted") +
  geom_line(data = line_posts_site %>% filter(type %in% c("Water", "Seston", "Larval", "Emergent", "Tetragnathidae")),
            aes(group = group),
            linetype = "dotdash") +
  geom_line(data = line_posts_site %>% filter(type %in% c("Water", "Biofilm", "Larval", "Emergent", "Tetragnathidae")),
            aes(group = group)) + 
  scale_fill_viridis_d(begin = 0.3) +
  facet_wrap(~site, nrow = 1) +
  NULL

ggsave(plot_sum_pfas_overall_site, file = "plots/plot_sum_pfas_site.jpg", width = 8, height = 5)

library(cowplot)
plot_sum = plot_grid(plot_sum_pfas_overall , plot_sum_pfas_overall_site + guides(fill = "none"), 
                     ncol = 1, rel_heights = c(1, 0.8))

ggsave(plot_sum, file = "plots/plot_sum.jpg", width = 9, height = 9)
# sum pfas old ----------------------------------------------------------------

# posts_taxon = readRDS(file = "posteriors/posts_taxon.rds")
# 
# raw_sum_pfas_overall = merged_d2 %>% 
#   as_tibble() %>% 
#   group_by(type, site, order, sample_id) %>% 
#   reframe(conc_ppb = sum(conc_ppb, na.rm = T)) %>% 
#   group_by(type, order, sample_id) %>% 
#   reframe(conc_ppb = median(conc_ppb))
# 
# posts_sumpfas = posts_taxon %>% 
#   # group_by(type, .draw, order, site, taxon) %>% # sum pfas
#   # reframe(.epred = sum(.epred)) %>% 
#   group_by(type, .draw, order, site) %>% # sum pfas
#   reframe(.epred = sum(.epred)) %>% 
#   group_by(type, order) %>% 
#   median_qi(.epred) 
# 
# line_posts = posts_sumpfas %>% 
#   group_by(type, order) %>% 
#   reframe(.epred = median(.epred)) %>% 
#   mutate(group = "lines")
# 
# posts_sumpfas_order = posts_taxon %>% 
#   filter(!is.na(taxon)) %>% 
#   group_by(type, .draw, order, site, taxon) %>% # sum pfas
#   reframe(.epred = sum(.epred)) %>% 
#   group_by(type, order, taxon) %>% 
#   median_qi(.epred) 
# 
# # plot_sum_pfas_overall = posts_sumpfas %>% 
# #   ggplot(aes(x = reorder(type, order),
# #              y = .epred)) + 
# #   geom_pointrange(size = 0.4, aes(ymin = .lower, ymax = .upper)) +
# #   scale_y_log10(labels = comma) +
# #   labs(y = "\u2211PFAS (ppb)",
# #        x = "") +
# #   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
# #         legend.position = "top",
# #         legend.text = element_text(size = 10),
# #         legend.title = element_blank()) +
# #   geom_line(data = line_posts %>% filter(type %in% c("Water", "Sediment", "Larval", "Emergent", "Tetragnathidae")),
# #             aes(group = group),
# #             linetype = "dashed") +
# #   geom_line(data = line_posts %>% filter(type %in% c("Water", "Detritus", "Larval", "Emergent", "Tetragnathidae")),
# #             aes(group = group),
# #             linetype = "dotted") +
# #   geom_line(data = line_posts %>% filter(type %in% c("Water", "Seston", "Larval", "Emergent", "Tetragnathidae")),
# #             aes(group = group),
# #             linetype = "dotdash") +
# #   geom_line(data = line_posts %>% filter(type %in% c("Water", "Biofilm", "Larval", "Emergent", "Tetragnathidae")),
# #             aes(group = group)) +
# #   geom_jitter(data = posts_sumpfas_order, aes(fill = taxon), 
# #               width = 0.08, height = 0, shape = 22) +
# #   scale_fill_viridis_d(begin = 0.3) +
# #   NULL
# # 
# # plot_sum_pfas_overall
# 
# ggsave(plot_sum_pfas_overall, file = "plots/plot_sum_pfas.jpg", width = 8, height = 5)
# 
# sum_pfas_table = posts_sumpfas %>% 
#   mutate(type = fct_reorder(type, order)) %>% 
#   group_by(type) %>% 
#   median_qi(.epred) %>% 
#   select(-.width, -.point, -.interval) %>% 
#   rename(median = .epred,
#          low95 = .lower, 
#          high95 = .upper) %>% 
#   mutate(units = "sum_ppb")
# 
# write_csv(sum_pfas_table, file = "tables/sum_pfas_table.csv")
# 
# posts_sumpfas_site = posts_taxon %>% 
#   # filter(.draw ==222) %>%
#   # group_by(type, site, pfas_type, .draw, order) %>% # average over taxa
#   # reframe(.epred = median(.epred)) %>%
#   # group_by(type, .draw, order) %>% # average over sites
#   # reframe(.epred = median(.epred)) %>%
#   group_by(type, site, .draw, order) %>% # sum pfas
#   reframe(.epred = sum(.epred)) %>% 
#   group_by(type, site, order) %>% 
#   median_qi(.epred) %>%
#   mutate(site = as.factor(site),
#          site = fct_relevel(site, "Hop Brook", "Russell Brook", "Ratlum Brook", "Burr Pond Brook", "Pequabuck River"),
#          type = fct_reorder(type, order))
# 
# 
# posts_sumpfas_order_site = posts_taxon %>% 
#   filter(!is.na(taxon)) %>% 
#   group_by(type, site, .draw, order, taxon) %>% # sum pfas
#   reframe(.epred = sum(.epred)) %>% 
#   group_by(type, site, order, taxon) %>% 
#   median_qi(.epred) %>%
#   mutate(site = as.factor(site),
#          site = fct_relevel(site, "Hop Brook", "Russell Brook", "Ratlum Brook", "Burr Pond Brook", "Pequabuck River"),
#          type = fct_reorder(type, order))
# 
# plot_sum_pfas_bysite = posts_sumpfas_site %>% 
#   ggplot(aes(x = reorder(type, order),
#              y = .epred)) + 
#   geom_pointrange(aes(ymin = .lower,ymax = .upper)) +
#   geom_line(aes(group = site)) +
#   facet_wrap(~site, nrow = 1) +
#   scale_y_log10(labels = comma, breaks = c(0.01, 0.1, 1, 10, 100, 1000)) +
#   labs(y = "\u2211PFAS (ppb)",
#        x = "",
#        color = "") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
#         legend.position = "top",
#         legend.text = element_text(size = 8),
#         legend.title = element_blank()) +
#   geom_jitter(data = posts_sumpfas_order, aes(fill = taxon), 
#               shape = 22, width = 0.08, height = 0) +
#   scale_fill_viridis_d(begin = 0.3) +
#   NULL
# 
# ggsave(plot_sum_pfas_bysite, file = "plots/plot_sum_pfas_bysite.jpg", width = 8, height = 5)
# 
# 
# write_csv(posts_sumpfas_site, file = "tables/sum_pfas_table_bysite.csv")
# 
# library(cowplot)
# 
# 
# plot_sum = plot_grid(plot_sum_pfas_overall , plot_sum_pfas_bysite + guides(fill = "none"), 
#                      ncol = 1, rel_heights = c(1, 0.8))
# 
# ggsave(plot_sum, file = "plots/plot_sum.jpg", width = 9, height = 9)
# 
# 
# posts_taxon %>% 
#   filter(.draw <= 100) %>% 
#   group_by(type, site, taxon, .draw, order) %>% 
#   reframe(.epred = sum(.epred)) %>% 
#   group_by(type, order, .draw) %>% 
#   reframe(.epred = sum(.epred)) %>% 
#   ggplot(aes(x = reorder(type, order),
#              y = .epred)) + 
#   stat_pointinterval() +
#   scale_y_log10() +
#   geom_jitter(data = raw_sum_pfas_overall, aes(y = conc_ppb),
#               size = 0.8, shape = 1, width = 0.1, height = 0) +
#   labs(y = "\u2211PFAS (ppb)",
#        x = "") +
#   NULL


# proportion pfas ---------------------------------------------------------


proportion_by_site = posts_taxon %>% 
  # filter(.draw < 1000) %>% 
  group_by(pfas_type, type, site) %>%
  reframe(.epred = median(.epred)) %>%
  mutate(type = as.factor(type),
         type = fct_relevel(type, "Water", "Sediment", "Biofilm", "Detritus", "Seston", "Larval", "Emergent",
                            "Tetragnathidae"),
         site = as.factor(site),
         site = fct_relevel(site, "Hop Brook", "Russell Brook", "Ratlum Brook", "Burr Pond Brook", "Pequabuck River"),
         pfas_type = fct_reorder(pfas_type, .epred)) %>% 
  left_join(pfas_orders) %>% 
  mutate(pfas_type = fct_reorder(pfas_type, order))

proportion_by_site_plot = proportion_by_site %>% 
  ggplot(aes(x = type, y = .epred)) + 
  geom_bar(position="fill", stat = "identity", aes(fill = pfas_type)) +
  facet_wrap(~site, nrow = 1) +
  # scale_fill_brewer(palette = "BrBG") +
  # scale_fill_viridis_d(direction = -1) +
  scale_fill_custom() +
  labs(x = "Sample Type", 
       fill = "PFAS") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0%", "25%", "50%", "75%", "100%")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
        axis.title.y = element_blank(),
        text = element_text(size = 7),
        legend.text = element_text(size = 7))  +
  theme(legend.key.size = unit(0.2, "cm")) +
  NULL

ggsave(proportion_by_site_plot, file = "plots/proportion_by_site.jpg", width = 6.5, height = 2)


proportion_overall = posts_taxon %>% 
  # filter(.draw < 1000) %>% 
  group_by(pfas_type, type) %>%
  reframe(.epred = median(.epred)) %>%
  mutate(type = as.factor(type),
         type = fct_relevel(type, "Water", "Sediment", "Biofilm", "Detritus", "Seston", "Larval", "Emergent",
                            "Tetragnathidae"),
         # pfas_type = fct_reorder(pfas_type, .epred)
         ) %>% 
  left_join(pfas_orders) %>% 
  mutate(pfas_type = fct_reorder(pfas_type, order))

proportion_overall_plot = proportion_overall %>% 
  ggplot(aes(x = type, y = .epred)) + 
  geom_bar(position="fill", stat = "identity", aes(fill = pfas_type)) +
  # scale_fill_brewer(palette = "BrBG") +
  # scale_fill_viridis_d(direction = -1) +
  scale_fill_custom() +
  labs(x = "Sample Type", 
       fill = "PFAS") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0%", "25%", "50%", "75%", "100%")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        text = element_text(size = 10)) +
  theme(legend.key.size = unit(0.4, "cm")) +
  NULL

ggsave(proportion_overall_plot, file = "plots/proportion_overall.jpg", width = 4, height = 3.5)


# detection ---------------------------------------------------------------

post_hus = posts_taxon %>% 
  ungroup %>% 
  filter(.draw <= 1000) %>% 
  group_by(type, pfas_type, .draw) %>% 
  reframe(hu = mean(hu)) %>% 
  pivot_wider(names_from = type, values_from = hu)

prob_labels = post_hus %>% 
  group_by(pfas_type) %>% 
  reframe(Water = mean(Water),
          Larval = mean(Larval))

prob_detect_average = post_hus %>% 
  ggplot(aes(x = 1 - Water, y = 1 - Larval, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  # facet_wrap(~pfas_type) +
  scale_color_custom() +
  # scale_y_log10() +
  # scale_x_log10() +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = prob_labels, aes(label = pfas_type, fill = pfas_type), 
             color = "white", size = 3) +
  scale_fill_custom() +
  labs(x = "Probability of Detection in Water",
       y = "Probability of Detection in Larval Insects")




# prob_detect_by site -----------------------------------------------------

post_hus_site = posts_taxon %>% 
  ungroup %>% 
  filter(.draw <= 1000) %>% 
  group_by(type, pfas_type, .draw, site) %>% 
  reframe(hu = mean(hu)) %>% 
  pivot_wider(names_from = type, values_from = hu)

prob_labels_site = post_hus_site %>% 
  group_by(pfas_type, site) %>% 
  reframe(Water = mean(Water),
          Larval = mean(Larval))

prob_detect_site = post_hus_site %>% 
  ggplot(aes(x = 1 - Water, y = 1 - Larval, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  facet_wrap(~site) +
  scale_color_custom() +
  # scale_y_log10() +
  # scale_x_log10() +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = prob_labels_site, aes(label = pfas_type, fill = pfas_type), 
             color = "white", 
             size = 2) +
  scale_fill_custom() +
  labs(x = "Probability of Detection in Water",
       y = "Probability of Detection in Larval Insects")


library(cowplot)
prob_detect_water_larvae = plot_grid(prob_detect_average, prob_detect_site, ncol = 1)

ggsave(prob_detect_water_larvae, file = "plots/prob_detect_water_larvae.jpg", width = 6.5, height = 8)


# prob_detect larvae to emerger -------------------------------------------

post_hus = posts_taxon %>% 
  ungroup %>% 
  filter(.draw <= 1000) %>% 
  group_by(type, pfas_type, .draw) %>% 
  reframe(hu = mean(hu)) %>% 
  pivot_wider(names_from = type, values_from = hu)

prob_labels_larv_emerge = post_hus %>% 
  group_by(pfas_type) %>% 
  reframe(Emergent = mean(Emergent),
          Larval = mean(Larval))

prob_detect_average_larv_emerge = post_hus %>% 
  ggplot(aes(x = 1 - Larval, y = 1 - Emergent, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  # facet_wrap(~pfas_type) +
  scale_color_custom() +
  # scale_y_log10() +
  # scale_x_log10() +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = prob_labels_larv_emerge, aes(label = pfas_type, fill = pfas_type), 
             color = "white", size = 3) +
  scale_fill_custom() +
  labs(x = "Probability of Detection in Larval Insects",
       y = "Probability of Detection in Emergent Insects")


# prob_detect_by site 

post_hus_site = posts_taxon %>% 
  ungroup %>% 
  filter(.draw <= 1000) %>% 
  group_by(type, pfas_type, .draw, site) %>% 
  reframe(hu = mean(hu)) %>% 
  pivot_wider(names_from = type, values_from = hu)

prob_labels_site_larv_emerge = post_hus_site %>% 
  group_by(pfas_type, site) %>% 
  reframe(Emergent = mean(Emergent),
          Larval = mean(Larval))

prob_detect_site_larv_emerge = post_hus_site %>% 
  ggplot(aes(x = 1 - Larval, y = 1 - Emergent, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  facet_wrap(~site) +
  scale_color_custom() +
  # scale_y_log10() +
  # scale_x_log10() +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = prob_labels_site_larv_emerge, aes(label = pfas_type, fill = pfas_type), 
             color = "white", 
             size = 2) +
  scale_fill_custom() +
  labs(x = "Probability of Detection in Larval Insects",
       y = "Probability of Detection in Emergent Insects")


library(cowplot)
prob_detect_larv_emerge = plot_grid(prob_detect_average_larv_emerge, prob_detect_site_larv_emerge, ncol = 1)

ggsave(prob_detect_larv_emerge, file = "plots/prob_detect_larv_emerge.jpg", width = 6.5, height = 8)



# prob_detect emerger to spider -------------------------------------------

post_hus = posts_taxon %>% 
  ungroup %>% 
  filter(.draw <= 1000) %>% 
  group_by(type, pfas_type, .draw) %>% 
  reframe(hu = mean(hu)) %>% 
  pivot_wider(names_from = type, values_from = hu)

prob_labels_spider_emerge = post_hus %>% 
  group_by(pfas_type) %>% 
  reframe(Emergent = mean(Emergent),
          Tetragnathidae = mean(Tetragnathidae))

prob_detect_average_spider_emerge = post_hus %>% 
  ggplot(aes(y = 1 - Tetragnathidae, x = 1 - Emergent, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  # facet_wrap(~pfas_type) +
  scale_color_custom() +
  # scale_y_log10() +
  # scale_x_log10() +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = prob_labels_spider_emerge, aes(label = pfas_type, fill = pfas_type), 
             color = "white", size = 3) +
  scale_fill_custom() +
  labs(y = "Probability of Detection in Tetragnathid Spiders",
       x = "Probability of Detection in Emergent Insects")


# prob_detect_by site 

post_hus_site = posts_taxon %>% 
  ungroup %>% 
  filter(.draw <= 1000) %>% 
  group_by(type, pfas_type, .draw, site) %>% 
  reframe(hu = mean(hu)) %>% 
  pivot_wider(names_from = type, values_from = hu)

prob_labels_site_spider_emerge = post_hus_site %>% 
  group_by(pfas_type, site) %>% 
  reframe(Emergent = mean(Emergent),
          Tetragnathidae = mean(Tetragnathidae))

prob_detect_site_spider_emerge = post_hus_site %>% 
  ggplot(aes(y = 1 - Tetragnathidae, x = 1 - Emergent, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  facet_wrap(~site) +
  scale_color_custom() +
  # scale_y_log10() +
  # scale_x_log10() +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = prob_labels_site_spider_emerge, aes(label = pfas_type, fill = pfas_type), 
             color = "white", 
             size = 2) +
  scale_fill_custom() +
  labs(y = "Probability of Detection in Tetragnathid Spiders",
       x = "Probability of Detection in Emergent Insects")


library(cowplot)
prob_detect_spider_emerge = plot_grid(prob_detect_average_spider_emerge, prob_detect_site_spider_emerge, ncol = 1)

ggsave(prob_detect_spider_emerge, file = "plots/prob_detect_spider_emerge.jpg", width = 6.5, height = 8)


