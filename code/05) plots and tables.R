library(tidybayes)
library(tidyverse)
library(brms)
library(janitor)
library(viridis)
library(ggh4x)
library(scales)
library(patchwork)
library(cowplot)
library(ggridges)

# wrangle models, data, functions -----------------------------------------

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


make_summary_table <- function(df, center = ".epred", lower = ".lower", upper = ".upper",
                               center_interval = "median_cri", digits = 1) {
  df %>%
    mutate(
      .center_val = round(.data[[center]], digits),
      .lower_val = round(.data[[lower]], digits),
      .upper_val = round(.data[[upper]], digits),
      center_interval = paste0(.center_val, " (", .lower_val, " to ", .upper_val, ")")) %>% 
    select(-.center_val, -.lower_val, -.upper_val)
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
  reframe(.epred = mean(.epred)) %>% 
  left_join(pfas_names) %>% 
  mutate(pfas_type = reorder(pfas_type, pfas_order))

posts_taxa_only = posts_taxon %>% 
  filter(!is.na(taxon)) 

# for lines
posts_taxon_type_summary = posts_taxon_type %>% 
  group_by(type, pfas_type, pfas_category, order, number_carbons) %>% 
  median_qi(.epred) 


# ppb (fig 1) -----------------------------------
posts_concentrations = mod_dat %>%
  select(-contains("conc_ppb")) %>%
  distinct() %>%
  # mutate(site = "new") %>% 
  add_epred_draws(seed = 20202, hg4_taxon, re_formula = NULL, dpar = T, ndraws = 500) %>%
  mutate(.epred = .epred*unique(hg4_taxon$data2$max_conc_ppb)) 

posts_concentrations_summary = posts_concentrations %>% 
  group_by(type, site, pfas_type, order, .draw) %>% # average over taxa
  reframe(.epred = mean(.epred)) %>%
  left_join(pfas_names) %>%
  mutate(pfas_type = reorder(pfas_type, pfas_order)) %>%
  group_by(pfas_type, order, type, .draw) %>%  # average over sites
  reframe(.epred = mean(.epred)) %>%
  group_by(pfas_type, order, type) %>%
  median_qi(.epred) %>% 
  left_join(pfas_names) %>% 
  mutate(pfas_type = reorder(pfas_type, pfas_order))
  
pfas_concentration = mod_dat %>% 
  ggplot(aes(x = reorder(type, order), y = conc_ppb + 1)) + 
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, shape = 1,
              size = 0.3) +
  geom_pointrange(data = posts_concentrations_summary, aes(y = .epred + 1, 
                                                       ymin = .lower + 1,
                                                       ymax = .upper + 1,
                                                       color = pfas_type), 
                  size = 0.3) + 
  scale_y_log10(breaks = c(1, 10, 100), 
                labels = c("0", "10", "100")) +
  # facet_wrap(pfas_category~reorder(pfas_type, number_carbons)) +
  facet_wrap2(~reorder(pfas_type, pfas_order)) +
  scale_color_custom() +
  labs(y = "PFAS Concentration (ppb)",
       x = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  guides(color = "none") +
  geom_line(data = posts_concentrations_summary %>% 
              mutate(group = "lines") %>% 
              filter(type %in% c("Water", "Sediment", "Larval", "Emergent", "Tetragnathidae")),
            aes(group = group, y = .epred + 1),
            linetype = "dashed") +
  geom_line(data = posts_concentrations_summary %>% 
              mutate(group = "lines") %>% filter(type %in% c("Water", "Detritus", "Larval", "Emergent", "Tetragnathidae")),
            aes(group = group, y = .epred + 1),
            linetype = "dotted") +
  geom_line(data = posts_concentrations_summary %>% 
              mutate(group = "lines") %>% filter(type %in% c("Water", "Seston", "Larval", "Emergent", "Tetragnathidae")),
            aes(group = group, y = .epred + 1),
            linetype = "dotdash") +
  geom_line(data = posts_concentrations_summary %>% 
              mutate(group = "lines") %>% filter(type %in% c("Water", "Biofilm", "Larval", "Emergent", "Tetragnathidae")),
            aes(group = group, y = .epred + 1)) +
  NULL

ggsave(pfas_concentration, file = "plots/ms_plots_tables/pfas_concentration_fig1.jpg", width = 6.5, height = 9)

# make tables
pfas_concentrations_fig1_summary = posts_concentrations_summary %>% 
  make_summary_table(digits = 3) %>% 
  select(pfas_type, type, center_interval) %>% 
  pivot_wider(names_from = type, values_from = center_interval) %>% 
  mutate(units = "ppb (posterior median and 95% CrI)")

write_csv(pfas_concentrations_fig1_summary, file = "plots/ms_plots_tables/pfas_concentration_fig1_summary.csv")

pfas_concentrations_bysite_summary = posts_concentrations %>% 
  group_by(type, site, pfas_type, order, .draw) %>% # average over taxa
  reframe(.epred = mean(.epred)) %>%
  left_join(pfas_names) %>%
  mutate(pfas_type = reorder(pfas_type, pfas_order)) %>%
  # group_by(pfas_type, order, type, .draw) %>%  # average over sites 
  # reframe(.epred = mean(.epred)) %>%
  group_by(pfas_type, order, type, site) %>%
  median_qi(.epred) %>% 
  left_join(pfas_names) %>% 
  mutate(pfas_type = reorder(pfas_type, pfas_order)) %>% 
  make_summary_table(digits = 3) %>% 
  select(site, pfas_type, type, center_interval) %>% 
  pivot_wider(names_from = type, values_from = center_interval) %>% 
  mutate(units = "ppb (posterior median and 95% CrI)")

write_csv(pfas_concentrations_bysite_summary, file = "plots/ms_plots_tables/pfas_concentration_bysite_summary.csv")

# sum ppb (fig 1) -------------------------------------------------------
mod1 = readRDS(file = "models/mod1.rds")
mod1_taxa = readRDS(file = "models/mod1_taxa.rds")

posts_sumpfas = mod1$data %>% 
  distinct(type, site) %>%
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8)) %>% 
  add_epred_draws(seed = 20202, mod1, re_formula = NULL) %>% 
  mutate(sum_ppb = (.epred - 0.0001)*unique(mod1_taxa$data2$mean_sum_ppb)) %>% 
  mutate(sum_ppb = case_when(sum_ppb < 0.000101 ~ 0, TRUE ~ sum_ppb - 0.0001)) %>% 
  group_by(type, .draw, order) %>% 
  reframe(sum_ppb = mean(sum_ppb))


posts_taxa_sumpfas = mod1_taxa$data2 %>% 
  distinct(site, type, taxon) %>%
  filter(taxon != "Megaloptera") %>% 
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8)) %>% 
  add_epred_draws(mod1_taxa) %>% 
  mutate(sum_ppb = (.epred - 0.0001)*unique(mod1_taxa$data2$mean_sum_ppb)) %>% 
  mutate(.epred = case_when(sum_ppb < 0.000101 ~ 0, TRUE ~ sum_ppb - 0.0001)) %>% 
  group_by(order, type, taxon, .draw) %>% 
  reframe(sum_ppb = mean(.epred)) 

posts_sumpfas_summary = posts_sumpfas %>% 
  group_by(type, order) %>% 
  median_qi(sum_ppb) 

line_posts = posts_sumpfas %>% 
  group_by(type, order) %>% 
  reframe(sum_ppb = median(sum_ppb)) %>% 
  mutate(group = "lines")

sum_pfas_overall = posts_sumpfas_summary %>% 
  ggplot(aes(x = reorder(type, order),
             y = sum_ppb)) +
  geom_pointrange(aes(ymin = .lower, ymax = .upper)) +
  scale_y_log10(labels = comma) +
  labs(y = expression("\u2211"[12]~PFAS~ (ppb)),
       x = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        legend.position = "top",
        legend.text = element_text(size = 10),
        legend.title = element_blank()) +
  stat_pointinterval(data = posts_taxa_sumpfas, aes(fill = taxon),
                     shape = 21, .width = 0,
                     position = position_jitter(width = 0.1, height = 0)) +
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

ggsave(sum_pfas_overall, file = "plots/ms_plots_tables/sum_pfas_fig1.jpg", width = 8, height = 5)

# make tables
sumpfas_type = posts_sumpfas_summary %>% 
  make_summary_table(center = "sum_ppb", digits = 1) %>% 
  select(type, center_interval) %>% 
  mutate(taxon = case_when(type == "Emergent" ~ "Overall",
                           type == "Larval" ~ "Overall"),
         group_order = 1) %>% 
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8))


sum_pfas_fig1_summary = sumpfas_type %>% 
  arrange(order, type, taxon) %>% 
  select(-order) %>% 
  mutate(units = "\u221112 ppb (posterior median and 95% CrI)")

write_csv(sum_pfas_fig1_summary, file = "plots/ms_plots_tables/sum_pfas_fig1_summary.csv")

sum_pfas_fig1_summary_taxa = posts_taxa_sumpfas %>% 
  group_by(type, taxon) %>% 
  median_qi(sum_ppb) %>% 
  make_summary_table(center = "sum_ppb") %>% 
  select(type, taxon, center_interval) %>% 
  pivot_wider(names_from = type, values_from = center_interval) %>% 
  select(taxon, Larval, Emergent) %>% 
  mutate(units = "\u221112 ppb (posterior median and 95% CrI)")

write_csv(sum_pfas_fig1_summary_taxa, file = "plots/ms_plots_tables/sum_pfas_fig1_summary_taxa.csv")


# proportion ppb (fig 1) -------------------------------------------------------------
posts_concentrations = mod_dat %>%
  select(-contains("conc_ppb")) %>%
  distinct() %>%
  # mutate(site = "new") %>% 
  add_epred_draws(seed = 20202, hg4_taxon, re_formula = NULL, dpar = T, ndraws = 500) %>%
  mutate(.epred = .epred*unique(hg4_taxon$data2$max_conc_ppb)) 

proportion_posts = posts_concentrations %>%
  filter(.draw <= 100) %>% 
  group_by(type, site, pfas_type, order, .draw) %>% # average over taxa
  reframe(.epred = mean(.epred)) %>%
  left_join(pfas_names) %>%
  mutate(pfas_type = reorder(pfas_type, pfas_order))  %>%
  group_by(pfas_type, order, type, .draw) %>%  # average over sites
  reframe(.epred = mean(.epred)) %>%
  group_by(pfas_type, type, .draw) %>%
  reframe(.epred = mean(.epred)) %>% 
  group_by(.draw, type) %>%
  mutate(total = sum(.epred)) %>% 
  mutate(.epred = .epred/total) %>%
  group_by(pfas_type, type) %>% 
  mean_qi(.epred) %>% 
  mutate(type = as.factor(type),
         type = fct_relevel(type, "Water", "Sediment", "Detritus", "Seston", "Biofilm", "Larval", "Emergent",
                            "Tetragnathidae"))

proportion_plot_fig1 = proportion_posts %>% 
  ggplot(aes(x = type, y = .epred*100, fill = pfas_type)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_custom() +
  labs(x = "Sample Type", 
       fill = "PFAS") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0%", "25%", "50%", "75%", "100%")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        axis.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        text = element_text(size = 10)) +
  theme(legend.key.size = unit(0.4, "cm")) +
  NULL

ggsave(proportion_plot_fig1, file = "plots/ms_plots_tables/proportion_plot_fig1.jpg", width = 4, height = 3.5)

# tables
proportion_plot_fig1_summary = proportion_posts %>% 
  make_summary_table() %>%
  rename(mean_cri = center_interval) %>% 
  select(pfas_type, type, mean_cri) %>% 
  pivot_wider(names_from = type, values_from = mean_cri) %>% 
  select(pfas_type, Water, Sediment, Detritus, Seston, Biofilm, Larval, Emergent, Tetragnathidae) %>% 
  mutate(units = "proportion of \u221112 pfas (posterior mean and 95% CrI)",
         notes = "Summarized as mean cri to ensure that columns sum to 1. Using medians, they don't sum to 1") 

write_csv(proportion_plot_fig1_summary , file = "plots/ms_plots_tables/proportion_plot_fig1_summary.csv")


# fig1_combined -----------------------------------------------------------
fig1a = sum_pfas_overall + labs(subtitle = "a)")
fig1b = pfas_concentration  + labs(subtitle = "b)")
fig1c = proportion_plot_fig1 + labs(subtitle = "c)")

fig1_left = fig1a/fig1c

fig1_combined = plot_grid(fig1_left, fig1b, ncol =2 )

ggsave(fig1_combined, file = "plots/ms_plots_tables/fig1_combined.jpg", width = 11, height = 9)

# sum ppb by site (fig s1a) ------------------------------------------------
# sumpfas by type and site
mod1 = readRDS(file = "models/mod1.rds")
mod1_taxa = readRDS(file = "models/mod1_taxa.rds")

posts_sum_pfas_site = mod1$data %>% 
  distinct(type, site) %>%
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8)) %>% 
  add_epred_draws(seed = 20202, mod1, re_formula = NULL) %>% 
  mutate(sum_ppb = (.epred - 0.0001)*unique(mod1_taxa$data2$mean_sum_ppb)) %>% 
  mutate(sum_ppb = case_when(sum_ppb < 0.000101 ~ 0, TRUE ~ sum_ppb - 0.0001)) %>% 
  group_by(type, site, order) %>% 
  median_qi(sum_ppb) %>% 
  mutate(site = as.factor(site),
         site = fct_relevel(site, "Hop Brook", "Russell Brook", "Ratlum Brook", "Burr Pond Brook", "Pequabuck River"))

line_posts_site = posts_sum_pfas_site  %>% 
  mutate(group = "lines") %>% 
  mutate(site = as.factor(site),
         site = fct_relevel(site, "Hop Brook", "Russell Brook", "Ratlum Brook", "Burr Pond Brook", "Pequabuck River"))

# sumpfas by order and life stage
posts_sum_pfas_taxa_site = mod1_taxa$data %>% 
  distinct(type, site, taxon) %>%
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 3,
                           type == "Detritus" ~ 4,
                           type == "Seston" ~ 5,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8)) %>% 
  add_epred_draws(seed = 20202, mod1_taxa, re_formula = NULL) %>% 
  mutate(sum_ppb = (.epred - 0.0001)*unique(mod1_taxa$data2$mean_sum_ppb)) %>% 
  mutate(sum_ppb = case_when(sum_ppb < 0.000101 ~ 0, TRUE ~ sum_ppb - 0.0001)) %>% 
  group_by(type, taxon, site, order) %>% 
  median_qi(sum_ppb) %>% 
  mutate(site = as.factor(site),
         site = fct_relevel(site, "Hop Brook", "Russell Brook", "Ratlum Brook", "Burr Pond Brook", "Pequabuck River"))

sum_ppb_by_site_fig_s1a = posts_sum_pfas_site %>% 
  ggplot(aes(x = reorder(type, order),
             y = sum_ppb)) +
  geom_pointrange(aes(ymin = .lower, ymax = .upper),
                  size = 0.2) +
  scale_y_log10(labels = comma) +
  labs(y = expression("\u2211"[12]~PFAS~ (ppb)),
       x = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        legend.text = element_text(size = 10),
        legend.title = element_blank()) +
  geom_point(data = posts_sum_pfas_taxa_site, aes(color = taxon, shape = taxon),
                     size = 1) +
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
  scale_color_viridis_d(begin = 0.2, end = 0.8) +
  facet_wrap(~site, nrow = 1) +
  NULL

ggsave(sum_ppb_by_site_fig_s1a, file = "plots/ms_plots_tables/sum_ppb_by_site_fig_s1a.jpg", width = 10, height = 3)


sum_ppb_by_site_summary = posts_sum_pfas_site %>% 
  make_summary_table(center = "sum_ppb", digits = 1) %>% 
  select(type, site, center_interval) %>% 
  pivot_wider(names_from = site, values_from = center_interval) %>% 
  mutate(units = "\u221112 ppb (posterior median and 95% CrI)") %>% 
  mutate(taxon = "Overall") %>% 
  select(type, taxon, everything())

sum_ppb_by_site_taxa_summary = posts_sum_pfas_taxa_site %>% 
  make_summary_table(center = "sum_ppb", digits = 1) %>% 
  select(site, type, taxon, center_interval) %>% 
  pivot_wider(names_from = site, values_from = center_interval) %>% 
  mutate(units = "\u221112 ppb (posterior median and 95% CrI)")

sum_ppb_by_site_figs1a_summary = bind_rows(sum_ppb_by_site_summary, sum_ppb_by_site_taxa_summary)

write_csv(sum_ppb_by_site_figs1a_summary, file = "plots/ms_plots_tables/sum_ppb_by_site_figs1a_summary.csv")


# sum partitioning (fig 2) ------------------------------------------------
mod1 = readRDS(file = "models/mod1.rds")

# raw ttf sum averaged over sites
raw_ttf_sum = as_tibble(mod1$data) %>% 
  mutate(sum_ppb = (sum_ppb_s_01 - 0.0001)*unique(mod1$data2$mean_sum_ppb)) %>% 
  # mutate(sum_ppb = case_when(sum_ppb < 0 ~ 0, TRUE ~ sum_ppb)) %>% # fixes the issue with back-transforming using -0.01 b/c some posterior estimates are less than that
  # select(-sum_ppb_s) %>%
  group_by(type, site) %>% 
  reframe(sum_ppb = mean(sum_ppb, na.rm = T)) %>%
  pivot_wider(names_from = type, values_from = sum_ppb) %>% 
  clean_names()  %>%
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

posts_ttf_sum = mod1$data %>% 
  distinct(type, site) %>%
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8)) %>% 
  add_epred_draws(seed = 20202, mod1, re_formula = NULL) %>% 
  mutate(sum_ppb = (.epred - 0.0001)*unique(mod1_taxa$data2$mean_sum_ppb)) %>% 
  mutate(sum_ppb = case_when(sum_ppb < 0.000101 ~ 0, TRUE ~ sum_ppb - 0.0001)) %>%
  group_by(type, .draw) %>% 
  reframe(sum_ppb = mean(sum_ppb)) %>% 
  pivot_wider(names_from = type, values_from = sum_ppb) %>%
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
  pivot_longer(cols = ends_with("ttf")) %>% 
  select(draw, name, value) %>% 
  group_by(draw, name) %>% 
  reframe(value = mean(value))

posts_ttf_summary_sum = posts_ttf_sum %>% 
  group_by(name) %>%
  median_qi(value) 

raw_nas_for_color_sum = raw_ttf_sum %>% 
  group_by(name) %>% 
  mutate(alpha = case_when(is.na(value) ~ "no raw data (estimate entirely from the priors/hierarcical effects)",
                           TRUE ~ "raw data (estimate from data and model)")) %>% 
  distinct(name, alpha) 

matrix_new_sum = posts_ttf_summary_sum %>% ungroup %>% distinct(name) %>% 
  separate(name, into = c("sort", "coefficient", "source","delete", "recipient", "ttf"), remove = F) %>% 
  mutate(path = paste0(source," --> ", recipient)) %>% 
  mutate(coefficient_name = case_when(coefficient == "baf" ~ "BAF", 
                                      coefficient == "kd" ~ "kd",
                                      coefficient == "mtf" ~ "MTF",
                                      coefficient == "ttf" ~ "TTF"),
         coefficient_name = as.factor(coefficient_name),
         coefficient_name = fct_relevel(coefficient_name, "kd", "BAF", "MTF", "TTF")) %>% 
  arrange(sort) %>% 
  mutate(path_order = 1:nrow(.)) 

sum_partitioning_fig2 = posts_ttf_summary_sum %>% 
  left_join(raw_nas_for_color_sum) %>% 
  left_join(matrix_new_sum) %>% 
  ggplot(aes(x = reorder(path, path_order), y = value)) + 
  geom_pointrange(aes(ymin = .lower, ymax = .upper), size = 0.2) +
  geom_point(data = raw_ttf_sum %>% 
               left_join(matrix_new_sum), color = "black", size = 0.4,
             shape = 21) +
  scale_y_log10(label = c("0.01", "0.1", "1", "10", "100", "1,000", "10,000"),
                breaks = c(0.01, 0.1, 1, 10, 100, 1000, 10000)) +
  guides(color = "none") +
  labs(y = "Trophic Transer Factor or Partitioning Coefficient",
       x = "Matrix") +
  geom_hline(aes(yintercept = 1), linetype = "dotted") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   vjust = 0.3)) +
  guides(alpha = "none") +
  NULL

ggsave(sum_partitioning_fig2, file = "plots/ms_plots_tables/sum_partitioning_fig2.jpg", width = 5, height = 5)

sum_partitioning_fig2_summary = posts_ttf_summary_sum %>% 
  make_summary_table(center = "value", digits = 1) %>% 
  separate(name, into = c("sort", "coefficient", "source","delete", "recipient", "ttf"), remove = F) %>% 
  mutate(path = paste0(source," --> ", recipient)) %>% 
  select(path, center_interval) %>% 
  mutate(units = "partitioning coefficients (median and 95% CrI)") %>% 
  rename(median_cri = center_interval)

write_csv(sum_partitioning_fig2_summary, file = "plots/ms_plots_tables/sum_partitioning_fig2_summary.csv")


# mef by taxon -----------------------------------------------

mod1_taxa = readRDS(file = "models/mod1_taxa.rds")

posts_ttf_sum_taxon = mod1_taxa$data %>% 
  distinct(type, site, taxon) %>%
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8)) %>% 
  add_epred_draws(seed = 20202, mod1_taxa, re_formula = NULL) %>% 
  mutate(sum_ppb = (.epred - 0.0001)*unique(mod1_taxa$data2$mean_sum_ppb)) %>% 
  mutate(sum_ppb = case_when(sum_ppb < 0.000101 ~ 0, TRUE ~ sum_ppb - 0.0001)) %>%
  group_by(type, taxon, .draw) %>% 
  reframe(sum_ppb = mean(sum_ppb)) %>% 
  pivot_wider(names_from = type, values_from = sum_ppb) %>%
  clean_names() %>% 
  mutate(g_mtf_larvae_to_emergent_ttf = emergent/larval) %>% 
  pivot_longer(cols = ends_with("ttf")) %>% 
  select(taxon, draw, larval, emergent, value) 

sum_partitioning_taxon = posts_ttf_sum_taxon %>%
  filter(taxon != "Megaloptera") %>% 
  rename(mef = value) %>% 
  pivot_longer(cols = c(larval, emergent, mef)) %>% 
  group_by(taxon, name) %>% 
  median_qi(value) %>% 
  make_summary_table(center = "value", digits = 1) %>% 
  select(taxon, name, center_interval) %>% 
  pivot_wider(names_from = name, values_from = center_interval) %>% 
  select(taxon, larval, emergent, mef)

write_csv(sum_partitioning_taxon, file = "plots/ms_plots_tables/sum_partitioning_taxon.csv")


# ppb perpfas by taxon ------------------------------------------------------------

perpfas_taxon_summary = posts_taxon %>% 
  filter(type %in% c("Emergent", "Larval")) %>% 
  group_by(type, taxon, pfas_type, pfas_order, .draw) %>% 
  reframe(.epred = mean(.epred)) %>%  # average over sites
  group_by(type, taxon, pfas_type, pfas_order) %>% 
  median_qi(.epred) %>% 
  make_summary_table(center = ".epred", digits = 2) %>% 
  select(pfas_type, type, taxon, center_interval) %>% 
  pivot_wider(names_from = type, values_from = center_interval) %>% 
  mutate(units = "ppb (median and 95% CrI)")

write_csv(perpfas_taxon_summary, file = "plots/ms_plots_tables/perpfas_taxon_summary.csv")

# compound partitioning (fig 2) -----------------------------------

posts_concentrations = mod_dat %>%
  select(-contains("conc_ppb")) %>%
  distinct() %>%
  # mutate(site = "new") %>% 
  add_epred_draws(seed = 20202, hg4_taxon, re_formula = NULL, dpar = T, ndraws = 500) %>%
  mutate(.epred = .epred*unique(hg4_taxon$data2$max_conc_ppb)) 

posts_ttf = posts_concentrations %>% 
  group_by(type, site, pfas_type, .draw) %>% 
  reframe(.epred = mean(.epred)) %>%  
  group_by(type, pfas_type, .draw) %>% 
  reframe(.epred = mean(.epred)) %>% 
  ungroup %>% 
  select(type, pfas_type, .draw, .epred) %>% 
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
  reframe(value = mean(value)) %>% 
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

matrix_new = posts_ttf_summary %>% 
  ungroup %>% 
  distinct(name, pfas_type) %>% 
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

compound_partitioning_fig2 = posts_ttf_summary %>% 
  left_join(raw_nas_for_color) %>% 
  left_join(matrix_new) %>% 
  left_join(pfas_names) %>% 
  mutate(pfas_type = reorder(pfas_type, pfas_order)) %>% 
  ggplot(aes(x = reorder(path, path_order), y = value, color = pfas_type)) + 
  # stat_pointinterval() + 
  geom_pointrange(aes(ymin = .lower, ymax = .upper, alpha = alpha), size = 0.2) +
  facet_wrap(~reorder(pfas_type, pfas_order)) +
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

ggsave(compound_partitioning_fig2, file = "plots/ms_plots_tables/compound_partitioning_fig2.jpg", width = 7, height = 7)

compound_partitioning_fig2_summary = posts_ttf_summary %>% 
  left_join(matrix_new) %>% 
  make_summary_table(center = "value", digits = 1) %>% 
  select(pfas_type, path, center_interval) %>% 
  pivot_wider(names_from = pfas_type, values_from = center_interval) %>% 
  mutate(units = "Partitioning (posterior median and 95% CrI)")

write_csv(compound_partitioning_fig2_summary, file = "plots/ms_plots_tables/compound_partitioning_fig2_summary.csv")

# fig2_combined -----------------------------------------------------------
library(cowplot)
partitioning_fig2_combined = plot_grid(sum_partitioning_fig2 + theme(axis.title.x = element_blank()),
          compound_partitioning_fig2 + theme(axis.title.x = element_blank()),
          nrow = 1)

ggsave(partitioning_fig2_combined, file = "plots/ms_plots_tables/partitioning_fig2_combined.jpg",
       dpi = 400, width = 12, height = 7)


# mef partitioning (fig 3) --------------------------------------------------

raw_mef_pertaxon = as_tibble(mod_dat) %>% 
  select(-conc_ppb_s) %>%
  group_by(type, pfas_type, site, taxon) %>% 
  reframe(conc_ppb = mean(conc_ppb, na.rm = T)) %>%
  pivot_wider(names_from = type, values_from = conc_ppb) %>% 
  clean_names()  %>%
  group_by(pfas_type, site, taxon) %>% 
  mutate(g_mtf_larvae_to_emergent_ttf = emergent/larval) %>% 
  pivot_longer(cols = ends_with("ttf"))  %>%
  mutate(across(where(is.numeric), ~ na_if(., 0))) %>%
  mutate(across(where(is.numeric), ~ na_if(., Inf))) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .))) %>% 
  filter(!is.na(taxon)) %>% 
  filter(taxon != 'Megaloptera') %>%
  mutate(alpha = case_when(is.na(value) ~ "no raw data (estimate entirely from the priors/hierarcical effects)",
                           TRUE ~ "raw data (estimate from data and model)")) %>% 
  mutate(taxon = as.factor(taxon),
         taxon = fct_relevel(taxon, "Diptera", "Ephemeroptera", "Odonata", "Plecoptera", "Trichoptera"))

posts_concentrations = mod_dat %>%
  select(-contains("conc_ppb")) %>%
  distinct() %>%
  # mutate(site = "new") %>% 
  add_epred_draws(seed = 20202, hg4_taxon, re_formula = NULL, dpar = T, ndraws = 500) %>%
  mutate(.epred = .epred*unique(hg4_taxon$data2$max_conc_ppb)) 

posts_mef = posts_concentrations %>% 
  filter(!is.na(taxon)) %>% 
  group_by(taxon, type, site, pfas_type, .draw) %>% 
  reframe(.epred = mean(.epred)) %>%  
  group_by(taxon, type, pfas_type, .draw) %>% 
  reframe(.epred = mean(.epred)) %>% 
  ungroup %>% 
  select(type, taxon, pfas_type, .draw, .epred) %>% 
  pivot_wider(names_from = type, values_from = .epred) %>% 
  clean_names() %>% 
  mutate(g_mtf_larvae_to_emergent_ttf = emergent/larval) %>% 
  pivot_longer(cols = ends_with("ttf")) %>% 
  left_join(pfas_names) %>% 
  mutate(pfas_type = reorder(pfas_type, pfas_order)) %>% 
  left_join(pfas_names) %>% 
  mutate(taxon = as.factor(taxon),
         taxon = fct_relevel(taxon, "Diptera", "Ephemeroptera", "Odonata", "Plecoptera", "Trichoptera"))

raw_nas_for_mef = raw_mef_pertaxon %>%
  mutate(alpha = case_when(is.na(value) ~ "no raw data (estimate entirely from the priors/hierarcical effects)",
                           TRUE ~ "raw data (estimate from data and model)")) %>% 
  ungroup %>% 
  distinct(pfas_type, taxon, alpha) %>% 
  mutate(taxon = as.factor(taxon),
         taxon = fct_relevel(taxon, "Diptera", "Ephemeroptera", "Odonata", "Plecoptera", "Trichoptera")) 

mef_by_taxon_and_pfas = posts_mef %>% 
  group_by(taxon, draw, pfas_category, pfas_type) %>% 
  reframe(mef = median(value, na.rm = T)) %>% 
  group_by(taxon, pfas_category, pfas_type) %>% 
  median_qi(mef) %>% 
  left_join(raw_nas_for_mef)

mef_by_taxon_fig3 = mef_by_taxon_and_pfas %>% 
  filter(taxon != "Megaloptera") %>% 
  ggplot(aes(x = taxon, y = mef, color = pfas_type)) + 
  geom_pointrange(aes(x = taxon, ymin = .lower, ymax = .upper, alpha = alpha),
                  size = 0.2) +
  geom_point(data = raw_mef_pertaxon, aes(x = taxon, y = value),
             size = 0.2,
             color = "black") +
  facet_wrap(~pfas_type) +
  scale_color_custom() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2,
                                   hjust = 1),
        axis.title.x = element_blank()) +
  scale_y_log10(label = c( "0.01", "0.1", "1", "10", "100"), 
                breaks = c(0.01,0.1, 1, 10, 100)) + 
  guides(color = "none",
         alpha = "none") +
  labs(y = "Metamorphic Enrichment Factor") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  coord_cartesian(ylim = c(0.01, 100))

ggsave(mef_by_taxon_fig3, file = "plots/ms_plots_tables/mef_by_taxon_fig3.jpg", width = 7, height = 7)


# body burden per pfas by taxon (fig 3) --------------------------------------------

mean_insect_mass = readRDS("posteriors/mean_insect_mass.rds") %>% clean_names() %>% 
  mutate(units = "grams")

insect_posts = merged_d2 %>% 
  distinct(site, type_taxon, taxon, type, max_conc, pfas_type) %>% 
  left_join(mean_insect_mass) %>% 
  mutate(type = as.factor(type),
         type = fct_relevel(type, "Larval", "Emergent")) %>% 
  add_epred_draws(hg4_taxon, ndraws = 500) %>% 
  mutate(.epred = .epred*max_conc) 

diffs = insect_posts %>% 
  filter(type %in% c("Emergent", "Larval")) %>% 
  mutate(ng_per_bug = (.epred/mean_gdw)) %>% 
  ungroup %>%
  select(site, taxon, type, ng_per_bug, .draw, pfas_type) %>%
  pivot_wider(names_from = type, values_from = ng_per_bug) %>% 
  group_by(pfas_type) %>% 
  mutate(diff = Emergent - Larval,
         ratio = Emergent / Larval) %>% 
  group_by(pfas_type) %>% 
  filter(!is.na(diff)) %>% 
  # filter(pfas_type == "PFOS") %>%
  group_by(taxon, .draw, pfas_type) %>% 
  reframe(diff = mean(diff),
          ratio = mean(ratio))  

diff_summary = diffs %>% group_by(taxon, pfas_type) %>% 
  median_qi(diff, .width = 0.5) %>% 
  mutate(taxon = as.factor(taxon),
         taxon = fct_relevel(taxon, "Plecoptera", "Odonata", "Ephemeroptera", "Trichoptera", "Diptera"))

change_in_body_burden_fig3 = diff_summary %>% 
  ggplot(aes(x = diff, y = pfas_type)) + 
  geom_point(aes(color = taxon), size = 0.8) + 
  geom_errorbarh(aes(xmin = .lower, xmax = .upper,
                     color = taxon), height = 0,
                 linewidth = 0.2) +
  # facet_grid2( pfas_type ~ .) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  xlim(-2000, 2000) +
  # guides(color = "none") +
  scale_color_brewer(type = "qual", palette = 7) +
  labs(x = "Change in body burden during metamorphosis (ng/individual)",
       y = "PFAS Compound",
       color = "") +
  geom_text(data = tibble(diff = 0, pfas_type = "8:2FTS"),
            label = "Larvae higher           Adults higher",
            aes(y = 12.5),
            size = 2) +
  theme(legend.text = element_text(size = 8),
        text = element_text(size = 8)) +
  NULL

ggsave(change_in_body_burden_fig3, file = "plots/ms_plots_tables/change_in_body_burden_fig3.jpg", 
       width = 6.5, height = 4)

change_in_body_burden_fig3_summary = diff_summary %>% 
  make_summary_table(center = "diff", digits = 1) %>% 
  select(taxon, pfas_type, center_interval) %>% 
  pivot_wider(names_from = taxon, values_from = center_interval) %>% 
  left_join(pfas_orders) %>% 
  arrange(order) %>% 
  select(-order)

write_csv(change_in_body_burden_fig3_summary, file = "plots/ms_plots_tables/change_in_body_burden_fig3_summary.csv")


ratio_burden_plot = diffs %>% 
  filter(!is.infinite(ratio)) %>% 
  mutate(taxon = as.factor(taxon),
         taxon = fct_relevel(taxon, "Plecoptera", "Odonata", "Ephemeroptera", "Trichoptera", "Diptera")) %>%
  group_by(taxon, pfas_type) %>% 
  mutate(taxon_median = median(diff)) %>% 
  mutate(ratio = ratio + 0.0001) %>% 
  ggplot(aes(x = ratio, y = pfas_type)) + 
  # ggridges::geom_density_ridges_gradient(aes(fill = after_stat(x)), scale = 0.7,
  #                               color = NA) +
  ggridges::geom_density_ridges(aes(fill = taxon)) +
  facet_wrap2(~taxon) +
  geom_vline(xintercept = 1) +
  scale_x_log10(limits = c(1e-04, 1e04),
                breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000")) +
  guides(fill = "none") +
  scale_fill_brewer(type = "qual", palette = 7) +
  labs(x = "Ratio of adult:larval PFAS body burden",
       y = "PFAS Compound",
       color = "") +
  geom_text(data = tibble(ratio = 1, pfas_type = "8:2FTS"),
            label = "Larvae higher                        Adults higher",
            aes(y = 13),
            size = 2) +
  theme(legend.text = element_text(size = 8),
        text = element_text(size = 8)) +
  NULL

ggsave(ratio_burden_plot, file = "plots/ms_plots_tables/ratio_burden_plot_fig3.jpg", 
       width = 6.5, height = 8)

ratio_burden_plot_fig3_summary = ratio_summary %>% 
  make_summary_table(center = "ratio", digits = 1) %>% 
  select(taxon, pfas_type, center_interval) %>% 
  pivot_wider(names_from = taxon, values_from = center_interval) %>% 
  left_join(pfas_orders) %>% 
  arrange(order) %>% 
  select(-order)

write_csv(ratio_burden_plot_fig3_summary, file = "plots/ms_plots_tables/ratio_burden_plot_fig3_summary.csv")


# body burden sum by taxon ---------------------------------------------------------
mean_insect_mass = readRDS("posteriors/mean_insect_mass.rds") %>% clean_names() %>% 
  mutate(units = "grams")

mod1_taxa = readRDS(file = "models/mod1_taxa.rds")

posts_ttf_sum_taxon = mod1_taxa$data %>% 
  distinct(type, site, taxon) %>%
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8)) %>% 
  left_join(mean_insect_mass %>% select(-type) %>% separate(type_taxon, into = c("type", "taxon"))) %>% 
  add_epred_draws(seed = 20202, mod1_taxa, re_formula = NULL) %>% 
  mutate(sum_ppb = (.epred - 0.0001)*unique(mod1_taxa$data2$mean_sum_ppb)) %>% 
  mutate(sum_ppb = case_when(sum_ppb < 0.000101 ~ 0, TRUE ~ sum_ppb - 0.0001)) %>% 
  mutate(sum_burden = sum_ppb/mean_gdw) %>% 
  group_by(type, taxon, .draw) %>% 
  reframe(sum_burden = mean(sum_burden, na.rm = T)) %>% 
  ungroup %>% 
  select(type, taxon, .draw, sum_burden) %>% 
  pivot_wider(names_from = type, values_from = sum_burden) %>%
  clean_names() 
  

diff_sum_burden = posts_ttf_sum_taxon %>% 
  mutate(g_mtf_larvae_to_emergent_ttf = emergent - larval) %>% 
  pivot_longer(cols = ends_with("ttf")) %>% 
  select(taxon, draw, larval, emergent, value) 

change_in_sum_body_burden = diff_sum_burden %>% 
  group_by(taxon) %>% 
  median_qi(value) %>% 
  make_summary_table(center = "value", digits = 0) %>% 
  filter(taxon != "Megaloptera") %>% 
  select(taxon, center_interval) %>% 
  mutate(metric = "Change in Body Burden adult:larval ng/individual")

write_csv(change_in_sum_body_burden, file = "plots/ms_plots_tables/change_in_sum_body_burden.csv")

ratio_sum_burden = posts_ttf_sum_taxon %>% 
  mutate(g_mtf_larvae_to_emergent_ttf = emergent / larval) %>% 
  pivot_longer(cols = ends_with("ttf")) %>% 
  select(taxon, draw, larval, emergent, value) 

ratio_in_sum_body_burden = ratio_sum_burden %>% 
  group_by(taxon) %>% 
  median_qi(value) %>% 
  make_summary_table(center = "value", digits = 1) %>% 
  filter(taxon != "Megaloptera") %>% 
  select(taxon, center_interval) %>% 
  mutate(metric = "Change in Body Burden adult:larval ng/individual")

write_csv(ratio_in_sum_body_burden, file = "plots/ms_plots_tables/ratio_in_sum_body_burden.csv")


# fig3_combined -----------------------------------------------------------
fig3a = mef_by_taxon_fig3

fig3b = change_in_body_burden_fig3 + theme(legend.position = c(0.8, 0.2))

fig3_combined = plot_grid(fig3a, fig3b, ncol = 2)

ggsave(fig3_combined, file = "plots/ms_plots_tables/fig3_combined.jpg", 
       width = 12, height = 6, dpi = 400)


# proportion ppb by site (fig s1b) -------------------------------------------------
proportion_posts_by_site = posts_concentrations %>%
  filter(.draw <= 100) %>% 
  group_by(type, site, pfas_type, order, .draw) %>% # average over taxa
  reframe(.epred = mean(.epred)) %>%
  left_join(pfas_names) %>%
  mutate(pfas_type = reorder(pfas_type, pfas_order))  %>%
  # group_by(pfas_type, order, type, .draw) %>%  # average over sites 
  # reframe(.epred = mean(.epred)) %>%
  group_by(pfas_type, type, site, .draw) %>%
  reframe(.epred = mean(.epred)) %>% 
  group_by(.draw, type, site) %>%
  mutate(total = sum(.epred)) %>% 
  mutate(.epred = .epred/total) %>%
  group_by(pfas_type, type, site) %>% 
  mean_qi(.epred) %>% 
  mutate(type = as.factor(type),
         type = fct_relevel(type, "Water", "Sediment", "Detritus", "Seston", "Biofilm", "Larval", "Emergent",
                            "Tetragnathidae"),
         site = as.factor(site),
         site = fct_relevel(site, "Hop Brook", "Russell Brook", "Ratlum Brook", "Burr Pond Brook", "Pequabuck River"),
         pfas_type = fct_reorder(pfas_type, .epred)) %>% 
  left_join(pfas_orders) %>% 
  mutate(pfas_type = fct_reorder(pfas_type, order))

proportion_plot_figs1b = proportion_posts_by_site %>% 
  ggplot(aes(x = type, y = .epred)) + 
  geom_bar(position="fill", stat = "identity", aes(fill = pfas_type)) +
  facet_wrap(~site, nrow = 1) +
  # scale_fill_brewer(palette = "BrBG") +
  # scale_fill_viridis_d(direction = -1) +
  scale_fill_custom() +
  labs(x = "Sample Type",
       y = "Percent Contribution",
       fill = "PFAS") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0%", "25%", "50%", "75%", "100%")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 10))  +
  theme(legend.key.size = unit(0.2, "cm")) +
  NULL

ggsave(proportion_plot_figs1b, file = "plots/ms_plots_tables/proportion_plot_figs1b.jpg", width = 6.5, height = 2)

proportion_plot_figs1b_summary = proportion_posts_by_site %>% 
  make_summary_table() %>%
  rename(mean_cri = center_interval) %>% 
  select(site, pfas_type, type, mean_cri) %>% 
  pivot_wider(names_from = type, values_from = mean_cri) %>% 
  mutate(units = "proportion of \u221112 pfas (posterior mean and 95% CrI)",
         notes = "Summarized as mean cri to ensure that columns sum to 1. Using medians, they don't sum to 1") %>% 
  arrange(site, pfas_type)

write_csv(proportion_plot_figs1b_summary, file = "plots/ms_plots_tables/proportion_plot_figs1b_summary.csv")

# combine figs1a and figs1b
library(patchwork)
figs1_combined = sum_ppb_by_site_fig_s1a/proportion_plot_figs1b
ggsave(figs1_combined, file = "plots/ms_plots_tables/figs1.jpg", width = 10, height = 5)

figs1_combined_alternate = (sum_ppb_by_site_fig_s1a + 
    theme(axis.text.x = element_blank()) + labs(subtitle = "a)"))/(proportion_plot_figs1b + 
                                                                     theme(strip.text = element_blank()) + 
                                                                     labs(subtitle = "b)"))
ggsave(figs1_combined_alternate, file = "plots/ms_plots_tables/figs1_alternate.jpg", width = 10, height = 6)


# probability of detection ------------------------------------------------
post_hus = posts_concentrations %>% 
  ungroup %>% 
  filter(.draw <= 1000) %>% 
  group_by(type, pfas_type, site, .draw) %>% 
  reframe(hu = mean(hu)) %>% 
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



# prob detect by site
post_hus_site = posts_concentrations %>% 
  ungroup %>% 
  filter(.draw <= 1000) %>% 
  group_by(type, pfas_type, site, .draw) %>% 
  # reframe(hu = mean(hu)) %>% 
  # group_by(type, pfas_type, .draw) %>% 
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


# prob_detect larvae to emerger
post_hus = posts_concentrations %>% 
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

post_hus_site = posts_concentrations %>% 
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



# prob_detect emerger to spider

post_hus = posts_concentrations %>% 
  ungroup %>% 
  filter(.draw <= 1000) %>% 
  group_by(type, pfas_type, site, .draw) %>% 
  reframe(hu = mean(hu)) %>% 
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

post_hus_site = posts_concentrations %>% 
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


# fig4_combined -----------------------------------------------------------

fig4a = prob_detect_average + labs(subtitle = "a)")
fig4b = prob_detect_average_larv_emerge + labs(subtitle = "b)")
fig4c = prob_detect_average_spider_emerge + labs(subtitle = "c)")

fig4_combined = plot_grid(fig4a, fig4b, fig4c, ncol = 1)

ggsave(fig4_combined, file = "plots/ms_plots_tables/fig4_combined.jpg", width = 6, height = 11,
       dpi = 400)

