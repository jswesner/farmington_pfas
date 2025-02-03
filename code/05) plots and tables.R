library(tidybayes)
library(tidyverse)
library(brms)
library(janitor)
library(viridis)
library(scales)
theme_set(theme_default())

# load data
merged_d2 = readRDS("data/merged_d2.rds") %>% 
  mutate(type_taxon = sample_type) %>% 
  separate(sample_type, into = c('type', 'taxon')) %>% 
  mutate(max_conc = max(conc_ppb, na.rm = T)) %>% 
  mutate(conc_ppb_s = conc_ppb/max_conc) %>% 
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8))

max_conc = unique(merged_d2$max_conc)

# load model
hg4_taxon = readRDS(file = "models/hg4_taxon.rds")

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
  mutate(conc_ppb = conc_ppb_s*max_conc)

posts_taxon = mod_dat  %>% 
  select(-contains("conc_ppb")) %>% 
  distinct() %>% 
  add_epred_draws(hg4_taxon, re_formula = NULL, dpar = T) %>%
  mutate(.epred = .epred*max_conc)

saveRDS(posts_taxon, file = "posteriors/posts_taxon.rds")

posts_taxon_type = posts_taxon %>% 
  group_by(type, site, pfas_type, order, .draw) %>% 
  reframe(.epred = mean(.epred))

posts_taxa_only = posts_taxon %>% 
  filter(!is.na(taxon)) 

pfas_concentration = mod_dat %>% 
  ggplot(aes(x = reorder(type, order), y = conc_ppb + 1)) + 
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, shape = 1,
              size = 0.3) +
  stat_pointinterval(data = posts_taxon_type, aes(y = .epred + 1), 
                     size = 0.3) + 
  scale_y_log10(breaks = c(1, 10, 100), 
                labels = c("0", "10", "100")) +
  facet_wrap(~pfas_type) +
  labs(y = "PFAS Concentration (ppb)",
       x = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  NULL

ggsave(pfas_concentration, file = "plots/pfas_concentration.jpg", width = 7, height = 5)


pfas_concentration_bysite = mod_dat %>% 
  ggplot(aes(x = reorder(type, order), y = conc_ppb + 1, color = site)) + 
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, shape = 1,
              size = 0.3) +
  stat_pointinterval(data = posts_taxon_type, aes(y = .epred + 1, color = site), 
                     size = 0.3) + 
  scale_y_log10(breaks = c(1, 10, 100), 
                labels = c("0", "10", "100")) +
  facet_wrap(~pfas_type) +
  labs(y = "PFAS Concentration (ppb)",
       x = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  NULL

ggsave(pfas_concentration_bysite, file = "plots/pfas_concentration_bysite.jpg", width = 8, height = 8)



# partitioning coefficients -----------------------------------------------



# get transfer factors ------------------
posts_ttf = posts_taxon_type %>% 
  ungroup %>% 
  select(type, site, pfas_type, .draw, .epred) %>% 
  pivot_wider(names_from = type, values_from = .epred) %>% 
  clean_names() %>% 
  mutate(a_kd_sediment_ttf = sediment/water,
         b_kd_biofilm_ttf = biofilm/water,
         a2_kd_seston_ttf = seston/water,
         b2_kd_detritus_ttf = detritus/water,
         c_baf_sediment_ttf = larval/sediment,
         d_baf_biofilm_ttf = larval/biofilm,
         e_baf_detritus_ttf = larval/detritus,
         f_baf_seston_ttf = larval/seston,
         g_mtf_ttf = emergent/larval,
         h_ttf_ttf = tetragnathidae/emergent) %>% 
  pivot_longer(cols = ends_with("ttf")) 

posts_ttf %>% 
  group_by(pfas_type, name) %>% 
  median_qi(value)

ttf_table_notaxa = posts_ttf %>% 
  group_by(draw, name, pfas_type) %>% 
  reframe(value = mean(value)) %>% 
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
         # g_mtf_ttf = emergent_trichoptera/value,
         # h_ttf_ttf = tetragnathidae/value
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
  reframe(value = mean(value)) %>% 
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



# plot_ttf ----------------------------------------------------------------


plot_ttf = posts_ttf %>% 
  group_by(draw, name, pfas_type) %>% 
  reframe(value = mean(value)) %>% 
  ggplot(aes(y = name, x = value, color = pfas_type)) + 
  stat_pointinterval() + 
  facet_wrap(~pfas_type) +
  # scale_color_brewer() +
  scale_x_log10(label = comma,
                limits = c(0.1, 1e5)) + 
  guides(color = "none") +
  labs(x = "Trophic Transer Factor",
       y = "Matrix",
       caption = "Figure X. Posterior distribution of trophic transfer factors (TTF)") +
  geom_vline(aes(xintercept = 1), linetype = "dotted") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   vjust = 1)) +
  coord_flip() +
  NULL

ggsave(plot_ttf, file = "plots/plot_ttf.jpg", width = 7, height = 8)


# sum pfas ----------------------------------------------------------------

posts_taxon = readRDS(file = "posteriors/posts_taxon.rds")

raw_sum_pfas_overall = merged_d2 %>% 
  as_tibble() %>% 
  group_by(type, site, order, sample_id) %>% 
  reframe(conc_ppb = sum(conc_ppb, na.rm = T)) %>% 
  group_by(type, order, sample_id) %>% 
  reframe(conc_ppb = median(conc_ppb))

posts_sumpfas = posts_taxon %>% 
  # filter(.draw ==222) %>%
  # group_by(type, site, pfas_type, .draw, order) %>% # average over taxa
  # reframe(.epred = mean(.epred)) %>%
  # group_by(type, .draw, order) %>% # average over sites
  # reframe(.epred = mean(.epred)) %>%
  group_by(type, site, .draw, order) %>% # sum pfas
  reframe(.epred = sum(.epred)) %>%
  group_by(type, .draw, order) %>% 
  reframe(.epred = median(.epred))

plot_sum_pfas_overall = posts_sumpfas %>% 
  ggplot(aes(x = reorder(type, order),
             y = .epred)) + 
  stat_pointinterval(size = 0.4) +
  scale_y_log10(labels = comma) +
  labs(y = "\u2211PFAS (ppb)",
       x = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
  NULL

plot_sum_pfas_overall

ggsave(plot_sum_pfas_overall, file = "plots/plot_sum_pfas.jpg", width = 8, height = 5)

sum_pfas_table = posts_sumpfas %>% 
  mutate(type = fct_reorder(type, order)) %>% 
  group_by(type) %>% 
  median_qi(.epred) %>% 
  select(-.width, -.point, -.interval) %>% 
  rename(median = .epred,
         low95 = .lower, 
         high95 = .upper) %>% 
  mutate(units = "sum_ppb")

write_csv(sum_pfas_table, file = "tables/sum_pfas_table.csv")

posts_sumpfas_site = posts_taxon %>% 
  # filter(.draw ==222) %>%
  # group_by(type, site, pfas_type, .draw, order) %>% # average over taxa
  # reframe(.epred = mean(.epred)) %>%
  # group_by(type, .draw, order) %>% # average over sites
  # reframe(.epred = mean(.epred)) %>%
  group_by(type, site, .draw, order) %>% # sum pfas
  reframe(.epred = sum(.epred)) %>% 
  group_by(type, site, order) %>% 
  median_qi(.epred) %>%
  mutate(site = as.factor(site),
         site = fct_relevel(site, "Hop Brook", "Russell Brook", "Ratlum Brook", "Burr Pond Brook", "Pequabuck River"),
         type = fct_reorder(type, order))

plot_sum_pfas_bysite = posts_sumpfas_site %>% 
  ggplot(aes(x = reorder(type, order),
             y = .epred)) + 
  geom_pointrange(aes(ymin = .lower,ymax = .upper)) +
  geom_line(aes(group = site)) +
  facet_wrap(~site, nrow = 1) +
  scale_y_log10(labels = comma, breaks = c(0.01, 0.1, 1, 10, 100, 1000)) +
  labs(y = "\u2211PFAS (ppb)",
       x = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
  NULL

ggsave(plot_sum_pfas_bysite, file = "plots/plot_sum_pfas_bysite.jpg", width = 8, height = 5)


write_csv(posts_sumpfas_site, file = "tables/sum_pfas_table_bysite.csv")

library(patchwork)

plot_sum = plot_sum_pfas_overall/plot_sum_pfas_bysite


ggsave(plot_sum, file = "plots/plot_sum.jpg", width = 8, height = 8)


posts_taxon %>% 
  filter(.draw <= 100) %>% 
  group_by(type, site, taxon, .draw, order) %>% 
  reframe(.epred = sum(.epred)) %>% 
  group_by(type, order, .draw) %>% 
  reframe(.epred = sum(.epred)) %>% 
  ggplot(aes(x = reorder(type, order),
             y = .epred)) + 
  stat_pointinterval() +
  scale_y_log10() +
  geom_jitter(data = raw_sum_pfas, aes(y = conc_ppb),
              size = 0.8, shape = 1, width = 0.1, height = 0) +
  labs(y = "\u2211PFAS (ppb)",
       x = "") +
  NULL

proportion_by_site = posts_taxon %>% 
  # filter(.draw < 1000) %>% 
  group_by(pfas_type, type, site) %>%
  reframe(.epred = median(.epred)) %>%
  mutate(type = as.factor(type),
         type = fct_relevel(type, "Water", "Sediment", "Biofilm", "Detritus", "Seston", "Larval", "Emergent",
                            "Tetragnathidae"),
         site = as.factor(site),
         site = fct_relevel(site, "Hop Brook", "Russell Brook", "Ratlum Brook", "Burr Pond Brook", "Pequabuck River"))

proportion_by_site_plot = proportion_by_site %>% 
  ggplot(aes(x = type, y = .epred)) + 
  geom_bar(position="fill", stat = "identity", aes(fill = pfas_type)) +
  facet_wrap(~site, nrow = 1) +
  scale_fill_brewer(palette = "BrBG") +
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
                            "Tetragnathidae"))

proportion_overall_plot = proportion_overall %>% 
  ggplot(aes(x = type, y = .epred)) + 
  geom_bar(position="fill", stat = "identity", aes(fill = pfas_type)) +
  scale_fill_brewer(palette = "BrBG") +
  labs(x = "Sample Type", 
       fill = "PFAS") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0%", "25%", "50%", "75%", "100%")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8)) +
  theme(legend.key.size = unit(0.4, "cm")) +
  NULL

ggsave(proportion_overall_plot, file = "plots/proportion_overall.jpg", width = 4, height = 3.5)
