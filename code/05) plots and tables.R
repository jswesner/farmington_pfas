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
    "PFHxA"="#c6dbee", "PFHpA" = "#a0cae0", "PFOA" = "#6aaed5", 
    "PFNA" = "#3e91c5","PFDA"="#1e71b5", "PFUnA" = "#166cb7", 
    "PFDoA"= "#0e2c6a", "PFBS" ="#fff7bc", "PFHxS" = "#fec34d",
    "PFOS" = "#FF9A28", "6:2FTS" = "#A1D99B", "8:2FTS" = "#41AB5D"
  ))
}

scale_color_custom <- function() {
  scale_color_manual(values = c(
    "PFHxA"="#c6dbee", "PFHpA" = "#a0cae0", "PFOA" = "#6aaed5", 
    "PFNA" = "#3e91c5","PFDA"="#1e71b5", "PFUnA" = "#166cb7", 
    "PFDoA"= "#0e2c6a", "PFBS" ="#fff7bc", "PFHxS" = "#fec34d",
    "PFOS" = "#FF9A28", "6:2FTS" = "#A1D99B", "8:2FTS" = "#41AB5D"
  ))
}

make_summary_table <- function(df, center = ".epred", lower = ".lower", upper = ".upper",
                               center_interval = "median_cri", digits = 2) {
  df %>%
    mutate(
      .center_val = signif(.data[[center]], digits),
      .lower_val = signif(.data[[lower]], digits),
      .upper_val = signif(.data[[upper]], digits),
      center_interval = paste0(.center_val, " (", .lower_val, " to ", .upper_val, ")")) %>% 
    select(-.center_val, -.lower_val, -.upper_val)
}

# load model
hg4_taxon = readRDS(file = "models/hg4_taxon.rds")
merged_d2 = hg4_taxon$data2

max_conc = unique(merged_d2$max_conc)

pfas_names = read_csv("data/pfas_names.csv")

mod_dat = hg4_taxon$data2 %>% 
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
  left_join(pfas_names %>% select(pfas_type, pfas_order)) %>% 
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
  group_by(type, pfas_type, pfas_category, order) %>% 
  median_qi(.epred) 


# ppb (fig 1) -----------------------------------
posts_concentrations = mod_dat %>%
  distinct(type_taxon, site, pfas_type, order, type) %>% 
  add_epred_draws(seed = 20202, hg4_taxon, re_formula = NULL, dpar = T, ndraws = 500) %>%
  mutate(.epred = .epred*unique(hg4_taxon$data2$max_conc_ppb)) 

posts_concentrations_summary = posts_concentrations %>% 
  group_by(type, site, pfas_type, order, .draw) %>% # average over taxa
  reframe(.epred = mean(.epred)) %>%
  left_join(pfas_names %>% select(pfas_type, pfas_order)) %>%
  mutate(pfas_type = reorder(pfas_type, pfas_order)) %>%
  group_by(pfas_type, order, type, .draw) %>%  # average over sites
  reframe(.epred = mean(.epred)) %>%
  group_by(pfas_type, order, type) %>%
  median_qi(.epred) %>% 
  left_join(pfas_names) %>% 
  mutate(pfas_type = reorder(pfas_type, pfas_order)) %>% 
  mutate(type = case_when(type == "Emergent" ~ "Adult", T ~ type)) %>%
  filter(!type %in% c("Detritus", "Seston", "Sediment"))
  
pfas_concentration = mod_dat %>%
  filter(!type %in% c("Detritus", "Seston", "Sediment")) %>% 
  mutate(type = case_when(type == "Emergent" ~ "Adult", T ~ type)) %>% 
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
  facet_wrap2(~reorder(pfas_type, pfas_order)) +
  scale_color_custom() +
  labs(y = "PFAS Concentration (ppb)",
       x = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  guides(color = "none") +
  geom_line(data = posts_concentrations_summary %>% 
              mutate(group = "lines") %>% filter(type %in% c("Water", "Biofilm", "Larval", "Adult", "Tetragnathidae")),
            aes(group = group, y = .epred + 1)) +
  NULL

ggsave(pfas_concentration, file = "plots/ms_plots_tables/pfas_concentration_fig1.jpg", width = 6.5, height = 9)

# make tables
pfas_concentrations_fig1_summary = posts_concentrations_summary %>% 
  make_summary_table(digits = 2) %>% 
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
  make_summary_table(digits = 2) %>% 
  select(site, pfas_type, type, center_interval) %>% 
  pivot_wider(names_from = type, values_from = center_interval) %>% 
  mutate(units = "ppb (posterior median and 95% CrI)")

write_csv(pfas_concentrations_bysite_summary, file = "plots/ms_plots_tables/pfas_concentration_bysite_summary.csv")

# sum ppb (fig 1) -------------------------------------------------------
mod1 = readRDS(file = "models/mod1.rds")
mod1_taxa = readRDS(file = "models/mod1_taxa.rds")

posts_sumpfas = mod1$data2 %>% 
  distinct(type, site, mean_sum_ppb) %>%
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8)) %>% 
  add_epred_draws(seed = 20202, mod1, re_formula = NULL) %>% 
  mutate(sum_ppb = (.epred - 0.0001)*mean_sum_ppb) %>% 
  mutate(sum_ppb = case_when(sum_ppb < 0.000101 ~ 0, TRUE ~ sum_ppb - 0.0001)) %>% 
  group_by(type, .draw, order) %>% 
  reframe(sum_ppb = mean(sum_ppb)) %>% 
  mutate(type = case_when(type == "Emergent" ~ "Adult", T ~ type))


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
  mutate(.epred = (.epred - 0.0001)*unique(mod1_taxa$data2$mean_sum_ppb)) %>% 
  mutate(.epred = case_when(.epred < 0.000101 ~ 0, TRUE ~ .epred - 0.0001)) %>% 
  group_by(order, type, taxon, .draw) %>% 
  reframe(sum_ppb = mean(.epred)) %>% 
  mutate(type = case_when(type == "Emergent" ~ "Adult", T ~ type))

posts_sumpfas_summary = posts_sumpfas %>% 
  group_by(type, order) %>% 
  median_qi(sum_ppb) 

posts_sumpfas_summary_taxa = posts_taxa_sumpfas %>% 
  group_by(type, order, taxon) %>% 
  median_qi(sum_ppb) 

# prob differences

posts_sumpfas_diff_taxa <- posts_taxa_sumpfas %>%
  group_by(taxon, type) %>% 
  mutate(taxon_id = cur_group_id()) %>% 
  inner_join(posts_taxa_sumpfas, 
             by = c("type", ".draw"),
             suffix = c("_a", "_b"),
             relationship = "many-to-many") %>%
  filter(taxon_a != taxon_b) %>%
  group_by(type, taxon_a, taxon_b) %>%
  reframe(
    prob_a_gt_b = mean(sum_ppb_a > sum_ppb_b)
  ) %>% 
  arrange(prob_a_gt_b) %>% 
  filter(prob_a_gt_b > 0.4999) # ensures a unidirectional comparison

line_posts = posts_sumpfas %>% 
  group_by(type, order) %>% 
  reframe(sum_ppb = median(sum_ppb)) %>% 
  mutate(group = "lines")

sum_pfas_overall = posts_sumpfas_summary %>%
  filter(!type %in% c("Detritus", "Seston", "Sediment")) %>% 
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
  geom_point(data = posts_sumpfas_summary_taxa %>%
               filter(!type %in% c("Detritus", "Seston", "Sediment")), aes(fill = taxon, shape = taxon),
             position = position_jitter(width = 0.1, height = 0),
             color = "black") +
  # stat_pointinterval(data = posts_taxa_sumpfas, aes(color = taxon,
  #                                                   shape = taxon), .width = 0.0,
  #                    position = position_jitter(width = 0.1, height = 0)) +
  # geom_line(data = line_posts %>% filter(type %in% c("Water", "Sediment", "Larval", "Emergent", "Tetragnathidae")),
  #           aes(group = group),
  #           linetype = "dashed") +
  # geom_line(data = line_posts %>% filter(type %in% c("Water", "Detritus", "Larval", "Emergent", "Tetragnathidae")),
  #           aes(group = group),
  #           linetype = "dotted") +
  # geom_line(data = line_posts %>% filter(type %in% c("Water", "Seston", "Larval", "Emergent", "Tetragnathidae")),
  #           aes(group = group),
  #           linetype = "dotdash") +
  geom_line(data = line_posts %>% filter(type %in% c("Water", "Biofilm", "Larval", "Adult", "Tetragnathidae")),
            aes(group = group)) + 
  scale_fill_brewer(type = "qual", palette = 7) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  NULL

ggsave(sum_pfas_overall, file = "plots/ms_plots_tables/sum_pfas_fig1.jpg", width = 8, height = 5)

# make tables
sumpfas_type = posts_sumpfas_summary %>% 
  make_summary_table(center = "sum_ppb", digits = 2) %>% 
  select(type, center_interval) %>% 
  mutate(taxon = case_when(type == "Adult" ~ "Overall",
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
  select(taxon, Larval, Adult) %>% 
  mutate(units = "\u221112 ppb (posterior median and 95% CrI)")

write_csv(sum_pfas_fig1_summary_taxa, file = "plots/ms_plots_tables/sum_pfas_fig1_summary_taxa.csv")


# magnitudes and probabilities
post_sumpfas_ratios = posts_sumpfas %>% 
  select(-order) %>% 
  pivot_wider(names_from = type, values_from = sum_ppb) %>% 
  mutate(a_bw = Biofilm/Water,
         b_lb = Larval/Biofilm,
         c_al = Adult/Larval,
         d_ta = Tetragnathidae/Adult) %>% 
  select(.draw, a_bw, b_lb, c_al, d_ta) %>% 
  pivot_longer(cols = -.draw)


post_sumpfas_ratios %>% 
  group_by(name) %>% 
  median_qi(value)

post_sumpfas_ratios %>% 
  group_by(name) %>% 
  reframe(prob = sum(value>1)/max(.draw))



# proportion ppb (fig 1) -------------------------------------------------------------
posts_concentrations = mod_dat %>%
  select(-contains("conc_ppb")) %>%
  distinct(type, taxon, type_taxon, site, pfas_type, order, pfas_order) %>%
  # mutate(site = "new")  %>% 
  add_epred_draws(seed = 20202, hg4_taxon, re_formula = NULL, dpar = T, ndraws = 500) %>%
  mutate(.epred = .epred*unique(hg4_taxon$data2$max_conc_ppb)) %>% 
  mutate(type = case_when(type == "Emergent" ~ "Adult", T ~ type))

proportion_posts = posts_concentrations %>%
  filter(.draw <= 1000) %>% 
  group_by(type, site, pfas_type, order, .draw) %>% # average over taxa
  reframe(.epred = mean(.epred)) %>%
  left_join(pfas_names %>% select(pfas_type, pfas_order))  %>%
  group_by(pfas_type, pfas_order, type, .draw) %>%  # average over sites
  reframe(.epred = mean(.epred)) %>%
  group_by(pfas_type, type, pfas_order, .draw) %>%
  reframe(.epred = mean(.epred)) %>% 
  group_by(.draw, type) %>%
  mutate(total = sum(.epred)) %>% 
  mutate(.epred = .epred/total) %>%
  group_by(pfas_type, pfas_order, type) %>% 
  mean_qi(.epred) %>% 
  mutate(type = as.factor(type),
         type = fct_relevel(type, "Water", "Biofilm", "Larval", "Adult",
                            "Tetragnathidae")) %>% 
  mutate(pfas_type = reorder(pfas_type, pfas_order)) 

proportion_plot_fig1 = proportion_posts %>%
  filter(!type %in% c("Detritus", "Seston", "Sediment")) %>% 
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
        legend.title = element_text(size = 8)) +
  theme(legend.key.size = unit(0.4, "cm")) +
  NULL

ggsave(proportion_plot_fig1, file = "plots/ms_plots_tables/proportion_plot_fig1.jpg", width = 4, height = 3.5)

# tables
proportion_plot_fig1_summary = proportion_posts %>% 
  mutate(.epred = .epred*100,
         .lower = .lower*100,
         .upper = .upper*100) %>% 
  make_summary_table(digits = 2) %>%
  rename(mean_cri = center_interval) %>% 
  select(pfas_type, type, mean_cri) %>% 
  pivot_wider(names_from = type, values_from = mean_cri) %>% 
  select(pfas_type, Water, Sediment, Detritus, Seston, Biofilm, Larval, Adult, Tetragnathidae) %>% 
  mutate(units = "proportion of \u221112 pfas (posterior mean and 95% CrI)",
         notes = "Summarized as mean cri to ensure that columns sum to 1. Using medians, they don't sum to 1") %>% 
  left_join(pfas_names) %>% 
  arrange(pfas_order) %>% 
  select(-pfas_order)

write_csv(proportion_plot_fig1_summary , file = "plots/ms_plots_tables/proportion_plot_fig1_summary.csv")


# fig1_combined -----------------------------------------------------------
fig1a = sum_pfas_overall + labs(subtitle = "a)")
fig1b = pfas_concentration  + labs(subtitle = "b)")
fig1c = proportion_plot_fig1 + labs(subtitle = "c)")

fig1_left = fig1a/fig1c

fig1_combined = plot_grid(fig1_left, fig1b, ncol =2 )

ggsave(fig1_combined, file = "plots/ms_plots_tables/fig1_combined.jpg", width = 11, height = 9)
ggsave(fig1b, file = "plots/ms_plots_tables/fig1b.jpg", width = 6.5, height = 9)


# proportion ppb per taxon ------------------------------------------------
posts_concentrations_taxon = mod_dat %>%
  select(-contains("conc_ppb")) %>%
  distinct(type, taxon, type_taxon, site, pfas_type, order, pfas_order) %>%
  filter(!is.na(taxon)) %>% 
  filter(taxon != "Megaloptera") %>% 
  mutate(taxon = case_when(type == "Tetragnathidae" ~ "Spider", TRUE ~ taxon))  %>%
  add_epred_draws(seed = 20202, hg4_taxon, re_formula = NULL, dpar = T, ndraws = 500) %>%
  mutate(.epred = .epred*unique(hg4_taxon$data2$max_conc_ppb)) %>% 
  mutate(type = case_when(type == "Emergent" ~ "Adult", T ~ type))

proportion_posts_taxon = posts_concentrations_taxon %>%
  filter(.draw <= 1000) %>% 
  # group_by(type, site, pfas_type, order, .draw) %>% # average over taxa
  # reframe(.epred = mean(.epred)) %>%
  left_join(pfas_names) %>%
  mutate(pfas_type = reorder(pfas_type, pfas_order))  %>%
  group_by(pfas_type,order, type, .draw, taxon) %>%  # average over sites
  reframe(.epred = mean(.epred)) %>%
  group_by(pfas_type, type, .draw, taxon) %>%
  reframe(.epred = mean(.epred)) %>% 
  group_by(.draw, type, taxon) %>%
  mutate(total = sum(.epred)) %>% 
  mutate(.epred = .epred/total) %>%
  group_by(pfas_type, type, taxon) %>% 
  mean_qi(.epred) %>% 
  mutate(type = as.factor(type))

proportion_plot_fig1_taxon = proportion_posts_taxon %>%
  mutate(type_label = case_when(type == "Larval" ~ "a) Larval",
                                TRUE ~ "b) Adult")) %>%
  left_join(pfas_names) %>%
  mutate(pfas_type = reorder(pfas_type, pfas_order)) %>% 
  ggplot(aes(x = taxon, y = .epred*100, fill = pfas_type)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_custom() +
  facet_wrap(~type_label, ncol = 2) +
  labs(x = "Sample Type", 
       fill = "PFAS") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0%", "25%", "50%", "75%", "100%")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        axis.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        text = element_text(size = 10),
        strip.text = element_text(hjust = 0)) +
  theme(legend.key.size = unit(0.4, "cm")) +
  NULL

ggsave(proportion_plot_fig1_taxon, file = "plots/ms_plots_tables/proportion_plot_fig1_taxon.jpg", width = 6.5, height = 4)

# tables
proportion_plot_fig1_summary_taxon = proportion_posts_taxon %>% 
  mutate(.epred = .epred*100,
         .lower = .lower*100,
         .upper = .upper*100) %>% 
  make_summary_table(digits = 2) %>%
  rename(mean_cri = center_interval) %>% 
  select(pfas_type, type, taxon, mean_cri) %>% 
  pivot_wider(names_from = taxon, values_from = mean_cri) %>% 
  mutate(units = "proportion of \u221112 pfas (posterior mean and 95% CrI)") %>% 
  left_join(pfas_names) %>% 
  arrange(pfas_order) %>% 
  select(-pfas_order)

write_csv(proportion_plot_fig1_summary_taxon, file = "plots/ms_plots_tables/proportion_plot_fig1_summary_taxon.csv")


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
  mutate(sum_ppb = (.epred - 0.0001)*unique(mod1$data2$mean_sum_ppb)) %>% 
  mutate(sum_ppb = case_when(sum_ppb < 0.000101 ~ 0, TRUE ~ sum_ppb - 0.0001)) %>% 
  group_by(type, site, order) %>% 
  median_qi(sum_ppb) %>% 
  mutate(site = as.factor(site),
         site = fct_relevel(site, "Hop Brook", "Russell Brook", "Ratlum Brook", "Burr Pond Brook", "Pequabuck River")) %>% 
  mutate(type = case_when(type == "Emergent" ~ "Adult", T ~ type))

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
         site = fct_relevel(site, "Hop Brook", "Russell Brook", "Ratlum Brook", "Burr Pond Brook", "Pequabuck River")) %>% 
  mutate(type = case_when(type == "Emergent" ~ "Adult", T ~ type))

sum_ppb_by_site_fig_s1a = posts_sum_pfas_site %>% 
  filter(type %in% c("Water", "Biofilm", "Larval", "Adult", "Tetragnathidae")) %>% 
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
  # geom_line(data = line_posts_site %>% filter(type %in% c("Water", "Sediment", "Larval", "Adult", "Tetragnathidae")),
  #           aes(group = group),
  #           linetype = "dashed") +
  # geom_line(data = line_posts_site %>% filter(type %in% c("Water", "Detritus", "Larval", "Adult", "Tetragnathidae")),
  #           aes(group = group),
  #           linetype = "dotted") +
  # geom_line(data = line_posts_site %>% filter(type %in% c("Water", "Seston", "Larval", "Adult", "Tetragnathidae")),
  #           aes(group = group),
            # linetype = "dotdash") +
  geom_line(data = line_posts_site %>% filter(type %in% c("Water", "Biofilm", "Larval", "Adult", "Tetragnathidae")),
            aes(group = group)) + 
  scale_color_viridis_d(begin = 0.2, end = 0.8) +
  facet_wrap(~site, nrow = 1) +
  NULL

ggsave(sum_ppb_by_site_fig_s1a, file = "plots/ms_plots_tables/sum_ppb_by_site_fig_s1a.jpg", width = 10, height = 3)


sum_ppb_by_site_summary = posts_sum_pfas_site %>% 
  make_summary_table(center = "sum_ppb", digits = 2) %>% 
  select(type, site, center_interval) %>% 
  pivot_wider(names_from = site, values_from = center_interval) %>% 
  mutate(units = "\u221112 ppb (posterior median and 95% CrI)") %>% 
  mutate(taxon = "Overall") %>% 
  select(type, taxon, everything())

sum_ppb_by_site_taxa_summary = posts_sum_pfas_taxa_site %>% 
  make_summary_table(center = "sum_ppb", digits = 2) %>% 
  select(site, type, taxon, center_interval) %>% 
  pivot_wider(names_from = site, values_from = center_interval) %>% 
  mutate(units = "\u221112 ppb (posterior median and 95% CrI)")

sum_ppb_by_site_figs1a_summary = bind_rows(sum_ppb_by_site_summary, sum_ppb_by_site_taxa_summary)

write_csv(sum_ppb_by_site_figs1a_summary, file = "plots/ms_plots_tables/sum_ppb_by_site_figs1a_summary.csv")


# sum partitioning (fig 2) ------------------------------------------------
mod1 = readRDS(file = "models/mod1.rds")
mod1_taxa = readRDS(file = "models/mod1_taxa.rds")

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
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .))) %>% 
  mutate(name = str_replace(name, "emergent", "adults"))

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
  mutate(sum_ppb = (.epred - 0.0001)*unique(mod1$data2$mean_sum_ppb)) %>% 
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
  mutate(name = str_replace(name, "emergent", "adults")) %>% 
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
  mutate(name = str_replace(name, "emergent", "adults")) %>% 
  left_join(raw_nas_for_color_sum %>% 
              mutate(name = str_replace(name, "emergent", "adults"))) %>% 
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
  make_summary_table(center = "value", digits = 2) %>% 
  separate(name, into = c("sort", "coefficient", "source","delete", "recipient", "ttf"), remove = F) %>% 
  mutate(path = paste0(source," --> ", recipient)) %>% 
  select(path, center_interval) %>% 
  mutate(units = "partitioning coefficients (median and 95% CrI)") %>% 
  rename(median_cri = center_interval)

write_csv(sum_partitioning_fig2_summary, file = "plots/ms_plots_tables/sum_partitioning_fig2_summary.csv")

sum_partitioning_fig2_summary_log10 = posts_ttf_summary_sum %>% 
  mutate(value = log10(value),
         .lower = log10(.lower),
         .upper = log10(.upper)) %>% 
  make_summary_table(center = "value", digits = 2) %>% 
  separate(name, into = c("sort", "coefficient", "source","delete", "recipient", "ttf"), remove = F) %>% 
  mutate(path = paste0(source," --> ", recipient)) %>% 
  select(path, center_interval) %>% 
  mutate(units = "partitioning coefficients (median and 95% CrI)") %>% 
  rename(median_cri = center_interval)

write_csv(sum_partitioning_fig2_summary_log10, file = "plots/ms_plots_tables/sum_partitioning_fig2_summary_log10.csv")


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
  make_summary_table(center = "value", digits = 2) %>% 
  select(taxon, name, center_interval) %>% 
  pivot_wider(names_from = name, values_from = center_interval) %>% 
  select(taxon, larval, emergent, mef)

write_csv(sum_partitioning_taxon, file = "plots/ms_plots_tables/sum_partitioning_taxon.csv")


# ppb perpfas by taxon ------------------------------------------------------------

perpfas_taxon_summary = posts_taxon %>% 
  filter(type %in% c("Adult", "Larval")) %>% 
  group_by(type, taxon, pfas_type, pfas_order, .draw) %>% 
  reframe(.epred = mean(.epred)) %>%  # average over sites
  group_by(type, taxon, pfas_type, pfas_order) %>% 
  median_qi(.epred) %>% 
  make_summary_table(center = ".epred", digits = 2) %>% 
  select(pfas_type, type, taxon, center_interval) %>% 
  pivot_wider(names_from = type, values_from = center_interval) %>% 
  mutate(units = "ppb (median and 95% CrI)")

write_csv(perpfas_taxon_summary, file = "plots/ms_plots_tables/perpfas_taxon_summary.csv")

perpfas_taxon_summary_lifestage = posts_taxon %>% 
  filter(!is.na(taxon)) %>% 
  ungroup %>% 
  group_by(taxon, type, pfas_type, .draw) %>%
  summarize(.epred = mean(.epred)) %>% 
  median_qi(.epred)

write_csv(perpfas_taxon_summary_lifestage, file = "plots/ms_plots_tables/perpfas_taxon_summary_lifestage.csv")


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
  mutate(path = str_replace(path, "emergent", "adults")) %>% 
  mutate(pfas_type = reorder(pfas_type, pfas_order)) %>% 
  ggplot(aes(x = reorder(path, path_order), y = value, color = pfas_type)) + 
  # stat_pointinterval() + 
  geom_pointrange(aes(ymin = .lower, ymax = .upper, alpha = alpha), size = 0.2) +
  facet_wrap(~reorder(pfas_type, pfas_order)) +
  geom_point(data = raw_ttf %>% 
               left_join(pfas_names) %>% 
               mutate(pfas_type = reorder(pfas_type, pfas_order)) %>% 
               left_join(matrix_new) %>% 
               mutate(path = str_replace(path, "emergent", "adults")), color = "black", size = 0.4) +
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
  make_summary_table(center = "value", digits = 2) %>% 
  select(pfas_type, path, center_interval) %>% 
  pivot_wider(names_from = pfas_type, values_from = center_interval) %>% 
  mutate(units = "Partitioning (posterior median and 95% CrI)")

write_csv(compound_partitioning_fig2_summary, file = "plots/ms_plots_tables/compound_partitioning_fig2_summary.csv")

compound_partitioning_fig2_summary_log10 = posts_ttf_summary %>% 
  mutate(value = log10(value),
         .lower = log10(.lower),
         .upper = log10(.upper)) %>% 
  left_join(matrix_new) %>% 
  make_summary_table(center = "value", digits = 2) %>% 
  select(pfas_type, path, center_interval) %>% 
  pivot_wider(names_from = pfas_type, values_from = center_interval) %>% 
  mutate(units = "Partitioning (posterior median and 95% CrI)")

write_csv(compound_partitioning_fig2_summary_log10, file = "plots/ms_plots_tables/compound_partitioning_fig2_summary_log10.csv")

# combine sum and per pfas partitioning
sum_partitioning_fig2_summary_log10 = read_csv(file = "plots/ms_plots_tables/sum_partitioning_fig2_summary_log10.csv")
compound_partitioning_fig2_summary_log10 = read_csv(file = "plots/ms_plots_tables/compound_partitioning_fig2_summary_log10.csv")

partitioning_table = compound_partitioning_fig2_summary_log10 %>% left_join(sum_partitioning_fig2_summary_log10 %>% select(path, median_cri) %>% 
                                                         rename(sum_pfas = median_cri)) %>% 
  select(-units) %>% 
  pivot_longer(cols = -path) %>% 
  filter(grepl("water", path)) %>% 
  pivot_wider(names_from = path, values_from = value) %>%
  rename(pfas_type = name) %>% 
  left_join(pfas_names %>% add_row(pfas_type = "sum_pfas", pfas_order = 0)) %>% 
  arrange(pfas_order) %>% 
  select(pfas_type, `water --> sediment`, `water --> detritus`, `water --> seston`, everything(), -pfas_order)

write_csv(partitioning_table, file = "plots/ms_plots_tables/partitioning_table.csv")

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

mef_by_taxon_fig3 = mef_by_taxon_and_pfas  %>% 
  filter(pfas_type %in% c("PFOA", "PFNA","PFUnA", "PFOS")) %>% 
  mutate(pfas_type = fct_relevel(pfas_type, "PFOA", "PFNA", "PFUnA")) %>% 
  filter(taxon != "Megaloptera") %>% 
  ggplot(aes(x = taxon, y = mef, color = pfas_type)) + 
  geom_pointrange(aes(x = taxon, ymin = .lower, ymax = .upper, alpha = alpha),
                  size = 0.2) +
  geom_point(data = raw_mef_pertaxon%>% 
               filter(pfas_type %in% c("PFOA", "PFNA","PFUnA", "PFOS")) %>% 
               mutate(pfas_type = fct_relevel(pfas_type, "PFOA", "PFNA", "PFUnA")), 
             aes(x = taxon, y = value),
             size = 0.2,
             color = "black") +
  facet_wrap(~pfas_type, nrow = 1) +
  scale_color_custom() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2,
                                   hjust = 1),
        axis.title.x = element_blank(),
        text = element_text(size = 8)) +
  scale_y_log10(label = c( "0.01", "0.1", "1", "10", "100"), 
                breaks = c(0.01,0.1, 1, 10, 100)) + 
  guides(color = "none",
         alpha = "none") +
  labs(y = "Bioamplification Factor") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  coord_cartesian(ylim = c(0.01, 100))

ggsave(mef_by_taxon_fig3, 
       file = "plots/ms_plots_tables/mef_by_taxon_fig3.jpg", 
       width = 7, height = 4)


# body burden per pfas by taxon (fig 3) --------------------------------------------

insect_posts = merged_d2 %>% 
  distinct(site, type_taxon, taxon, type, max_conc, pfas_type) %>% 
  mutate(type = as.factor(type),
         type = fct_relevel(type, "Larval", "Emergent")) %>% 
  add_epred_draws(hg4_taxon, ndraws = 500) %>% 
  mutate(.epred = .epred*max_conc) 

insect_mass_posts = readRDS(file = "posteriors/insect_mass_posts.rds")

diffs = insect_posts %>% 
  filter(type %in% c("Emergent", "Larval")) %>% 
  left_join(insect_mass_posts %>% ungroup %>% select(type, taxon, site, mean_gww, .draw)) %>% 
  mutate(ng_per_bug = (.epred*mean_gww)) %>% 
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
  filter(pfas_type %in% c("PFOA", "PFNA","PFUnA", "PFOS")) %>% 
  ggplot(aes(x = diff, y = pfas_type)) + 
  geom_point(aes(color = taxon), size = 0.8) + 
  geom_errorbarh(aes(xmin = .lower, xmax = .upper,
                     color = taxon), height = 0,
                 linewidth = 0.2) +
  # facet_grid2( pfas_type ~ .) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  xlim(-1, 1) +
  # guides(color = "none") +
  scale_color_brewer(type = "qual", palette = 7) +
  labs(x = "Change in body burden during metamorphosis (ng/individual)",
       y = "PFAS Compound",
       color = "") +
  geom_text(data = tibble(diff = 0, pfas_type = "8:2FTS"),
            label = "Larvae higher           Adults higher",
            aes(y = 4.5),
            size = 2) +
  theme(legend.text = element_text(size = 8),
        text = element_text(size = 8)) +
  NULL

ggsave(change_in_body_burden_fig3, file = "plots/ms_plots_tables/change_in_body_burden_fig3.jpg", 
       width = 5, height = 3)

change_in_body_burden_fig3_summary = diff_summary %>% 
  make_summary_table(center = "diff", digits = 2) %>% 
  select(taxon, pfas_type, center_interval) %>% 
  pivot_wider(names_from = taxon, values_from = center_interval) %>% 
  left_join(pfas_names) %>% 
  arrange(pfas_order) %>% 
  select(-pfas_order)

write_csv(change_in_body_burden_fig3_summary, file = "plots/ms_plots_tables/change_in_body_burden_fig3_summary.csv")


ratio_burden_plot = diffs  %>% 
  filter(pfas_type %in% c("PFOA", "PFNA","PFUnA", "PFOS"))  %>% 
  mutate(pfas_type = fct_relevel(pfas_type, "PFOA", "PFNA", "PFUnA", "PFOS")) %>%  
  filter(!is.infinite(ratio)) %>% 
  mutate(taxon = as.factor(taxon),
         # taxon = fct_relevel(taxon, "Plecoptera", "Odonata", "Ephemeroptera", "Trichoptera", "Diptera")
         ) %>%
  group_by(taxon, pfas_type) %>% 
  mutate(taxon_median = median(diff)) %>% 
  mutate(ratio = ratio + 0.0001) %>% 
  ggplot(aes(x = ratio, y = taxon)) + 
  # ggridges::geom_density_ridges_gradient(aes(fill = after_stat(x)), scale = 0.7,
  #                               color = NA) +
  ggridges::geom_density_ridges(fill = "grey", scale = 1) +
  facet_wrap2(~fct_relevel(pfas_type, "PFOA", "PFNA", "PFUnA", "PFOS"), nrow = 1) +
  geom_vline(xintercept = 1) +
  scale_x_log10(limits = c(1e-04, 1e04),
                breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000")) +
  guides(fill = "none") +
  # scale_fill_grey(start = 0.6, end = 0.4) +
  labs(x = "Ratio of adult:larval PFAS body burden",
       y = "PFAS Compound",
       color = "") +
  geom_text(data = tibble(ratio = 1) %>% 
              expand_grid(pfas_type =  c("PFOA", "PFNA", "PFUnA", "PFOS")),
            label = "Larvae higher                Adults higher",
            aes(y =6),
            size = 2) +
  theme(legend.text = element_text(size = 8),
        text = element_text(size = 8)) +
  NULL

ggsave(ratio_burden_plot, file = "plots/ms_plots_tables/ratio_burden_plot_fig3.jpg", 
       width = 6.5, height = 4)

ratio_burden_plot_fig3_summary = diffs %>% 
  group_by(taxon, pfas_type) %>% 
  median_qi(ratio) %>% 
  make_summary_table(center = "ratio", digits = 2) %>% 
  select(taxon, pfas_type, center_interval) %>% 
  pivot_wider(names_from = taxon, values_from = center_interval) %>% 
  left_join(pfas_names) %>% 
  arrange(pfas_order) %>% 
  select(-pfas_order)

write_csv(ratio_burden_plot_fig3_summary, file = "plots/ms_plots_tables/ratio_burden_plot_fig3_summary.csv")


# body burden sum by taxon ---------------------------------------------------------

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
  left_join(readRDS(file = "posteriors/insect_mass_posts.rds") %>% ungroup %>% select(type, taxon, site, mean_gww, .draw)) %>% 
  mutate(sum_burden = sum_ppb*mean_gww) %>% 
  group_by(type, taxon, .draw) %>% 
  reframe(sum_burden = mean(sum_burden, na.rm = T)) %>% 
  ungroup %>% 
  select(type, taxon, .draw, sum_burden) %>% 
  pivot_wider(names_from = type, values_from = sum_burden) %>%
  clean_names() 
  
sum_burden_taxa_lifestage = posts_ttf_sum_taxon %>% 
  pivot_longer(cols = c(emergent, larval)) %>% 
  group_by(taxon, name) %>% 
  median_qi(value) %>% 
  make_summary_table(center = "value", digits = 2) %>% 
  select(taxon, name, center_interval) %>% 
  pivot_wider(names_from = name, values_from = center_interval) %>% 
  filter(taxon != "Megaloptera")

write_csv(sum_burden_taxa_lifestage, file = "plots/ms_plots_tables/sum_burden_taxa_lifestage.csv")

diff_sum_burden = posts_ttf_sum_taxon %>% 
  mutate(g_mtf_larvae_to_emergent_ttf = emergent - larval) %>% 
  pivot_longer(cols = ends_with("ttf")) %>% 
  select(taxon, draw, larval, emergent, value) 

change_in_sum_body_burden = diff_sum_burden %>% 
  group_by(taxon) %>% 
  median_qi(value) %>% 
  make_summary_table(center = "value", digits = 2) %>% 
  filter(taxon != "Megaloptera") %>% 
  select(taxon, center_interval) %>% 
  mutate(metric = "Change in Body Burden adult:larval ng/individual")

write_csv(change_in_sum_body_burden, file = "plots/ms_plots_tables/change_in_sum_body_burden.csv")

ratio_sum_burden = posts_ttf_sum_taxon %>% 
  mutate(g_mtf_larvae_to_emergent_ttf = emergent / larval) %>% 
  pivot_longer(cols = ends_with("ttf")) %>% 
  select(taxon, draw, larval, emergent, value) 

ratio_sum_burden %>% 
  group_by(taxon) %>% 
  reframe(prob_negative = sum(value < 1)/max(draw))

ratio_in_sum_body_burden = ratio_sum_burden %>% 
  group_by(taxon) %>% 
  median_qi(value) %>% 
  make_summary_table(center = "value", digits = 2) %>% 
  filter(taxon != "Megaloptera") %>% 
  select(taxon, center_interval) %>% 
  mutate(metric = "Change in Body Burden adult:larval ng/individual")

write_csv(ratio_in_sum_body_burden, file = "plots/ms_plots_tables/ratio_in_sum_body_burden.csv")


# fig3_combined -----------------------------------------------------------
fig3a = mef_by_taxon_fig3

fig3b = ratio_burden_plot 

fig3_combined = plot_grid(fig3a, fig3b, nrow = 2, align = "v")

ggsave(fig3_combined, file = "plots/ms_plots_tables/fig3_combined.jpg", 
       width = 8.5, height = 6, dpi = 400)


# proportion ppb by site (fig s1b) -------------------------------------------------
proportion_posts_by_site = posts_concentrations %>%
  filter(!type %in% c("Detritus", "Seston", "Sediment")) %>% 
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
  mutate(type = case_when(type == "Emergent" ~ "Adult", TRUE ~ type)) %>% 
  mutate(type = as.factor(type),
         type = fct_relevel(type, "Water", "Biofilm", "Larval", "Adult",
                            "Tetragnathidae"),
         site = as.factor(site),
         site = fct_relevel(site, "Hop Brook", "Russell Brook", "Ratlum Brook", "Burr Pond Brook", "Pequabuck River")) %>% 
  left_join(pfas_names) %>% 
  mutate(pfas_type = fct_reorder(pfas_type, pfas_order)) 

proportion_plot_figs1b = proportion_posts_by_site %>% 
  ggplot(aes(x = type, y = .epred)) + 
  geom_bar(position="fill", stat = "identity", aes(fill = pfas_type)) +
  facet_wrap(~site, nrow = 1) +
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
  left_join(pfas_names) %>%
  group_by(type, pfas_type, .draw) %>% 
  reframe(hu = mean(hu)) %>%
  pivot_wider(names_from = type, values_from = hu) 

prob_labels = post_hus %>% 
  group_by(pfas_type) %>% 
  reframe(Water = median(Water),
          Larval = median(Larval),
          Biofilm = median(Biofilm),
          Seston = median(Seston),
          Detritus = median(Detritus))

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
             color = "white", size = 2.3) +
  scale_fill_custom() +
  labs(x = "P(Detect) in Water",
       y = "P(Detect) in Larval Insects")

prob_detect_water_biofilm = post_hus %>% 
  ggplot(aes(x = 1 - Water, y = 1 - Biofilm, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  # facet_wrap(~pfas_type) +
  scale_color_custom() +
  # scale_y_log10() +
  # scale_x_log10() +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = prob_labels, aes(label = pfas_type, fill = pfas_type), 
             color = "white", size = 2.3) +
  scale_fill_custom() +
  labs(x = "P(Detect) in Water",
       y = "P(Detect) in Biofilm")

prob_detect_water_spider = post_hus %>% 
  ggplot(aes(x = 1 - Water, y = 1 - Tetragnathidae, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  # facet_wrap(~pfas_type) +
  scale_color_custom() +
  # scale_y_log10() +
  # scale_x_log10() +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = prob_labels, aes(label = pfas_type, fill = pfas_type), 
             color = "white", size = 2.3) +
  scale_fill_custom() +
  labs(x = "P(Detect) in Water",
       y = "P(Detect) in Tetragnathidae")

ggsave(prob_detect_water_spider, file = "plots/ms_plots_tables/prob_detect_water_spider.jpg", width = 5, height = 5)


prob_detect_biofilm_larvae = post_hus %>% 
  ggplot(aes(x = 1 - Biofilm, y = 1 - Larval, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  # facet_wrap(~pfas_type) +
  scale_color_custom() +
  # scale_y_log10() +
  # scale_x_log10() +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = prob_labels, aes(label = pfas_type, fill = pfas_type), 
             color = "white", size = 2.3) +
  scale_fill_custom() +
  labs(x = "P(Detect) in Biofilm",
       y = "P(Detect) in Larvae")

prob_detect_seston_larvae = post_hus %>% 
  ggplot(aes(x = 1 - Seston, y = 1 - Larval, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  # facet_wrap(~pfas_type) +
  scale_color_custom() +
  # scale_y_log10() +
  # scale_x_log10() +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = prob_labels, aes(label = pfas_type, fill = pfas_type), 
             color = "white", size = 2.3) +
  scale_fill_custom() +
  labs(x = "P(Detect) in Seston",
       y = "P(Detect) in Larvae")

prob_detect_detritus_larvae = post_hus %>% 
  ggplot(aes(x = 1 - Detritus, y = 1 - Larval, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  # facet_wrap(~pfas_type) +
  scale_color_custom() +
  # scale_y_log10() +
  # scale_x_log10() +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = prob_labels, aes(label = pfas_type, fill = pfas_type), 
             color = "white", size = 2.3) +
  scale_fill_custom() +
  labs(x = "P(Detect) in Detritus",
       y = "P(Detect) in Larvae")

# prob detect by site
post_hus_site = posts_concentrations %>% 
  ungroup %>% 
  filter(.draw <= 1000) %>% 
  group_by(type, pfas_type, site, .draw) %>% 
  # reframe(hu = mean(hu)) %>% 
  # group_by(type, pfas_type, .draw) %>% 
  reframe(hu = mean(hu)) %>% 
  pivot_wider(names_from = type, values_from = hu) %>% 
  mutate(Adult = Emergent)

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
             size = 2.3) +
  scale_fill_custom() +
  labs(x = "P(Detect) in Water",
       y = "P(Detect) in Larval Insects")


library(cowplot)
prob_detect_water_larvae = plot_grid(prob_detect_average, prob_detect_site, ncol = 1)

ggsave(prob_detect_water_larvae, file = "plots/prob_detect_water_larvae.jpg", width = 6.5, height = 8)


prob_detect_rev = plot_grid(prob_detect_biofilm_larvae, prob_detect_seston_larvae, prob_detect_detritus_larvae)

ggsave(prob_detect_rev, file = "plots/temporary/prob_detect_rev.jpg", width = 6.5,
       height = 7)

# prob_detect larvae to emerger

prob_labels_larv_emerge = post_hus %>% 
  group_by(pfas_type) %>% 
  reframe(Adult = mean(Emergent),
          Larval = mean(Larval)) %>% 
  mutate(Emergent = Adult)

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
             color = "white", size = 2.3) +
  scale_fill_custom() +
  labs(x = "P(Detect) in Larval Insects",
       y = "P(Detect) in Adult Insects")


# prob_detect_by site 

post_hus_site = posts_concentrations %>% 
  ungroup %>% 
  filter(.draw <= 1000) %>% 
  group_by(type, pfas_type, .draw, site) %>% 
  reframe(hu = mean(hu)) %>% 
  pivot_wider(names_from = type, values_from = hu) %>% 
  mutate(Adult = Emergent)

prob_labels_site_larv_emerge = post_hus_site %>% 
  group_by(pfas_type, site) %>% 
  reframe(Adult = mean(Emergent),
          Larval = mean(Larval)) %>% 
  mutate(Emergent = Adult)

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
             size = 2.3) +
  scale_fill_custom() +
  labs(x = "P(Detect) in Larval Insects",
       y = "P(Detect) in Adult Insects")


library(cowplot)
prob_detect_larv_emerge = plot_grid(prob_detect_average_larv_emerge, prob_detect_site_larv_emerge, ncol = 1)

ggsave(prob_detect_larv_emerge, file = "plots/prob_detect_larv_emerge.jpg", width = 6.5, height = 8)



# prob_detect emerger to spider

prob_labels_spider_emerge = post_hus %>% 
  group_by(pfas_type) %>% 
  reframe(Adult = mean(Emergent),
          Tetragnathidae = mean(Tetragnathidae))%>% 
  mutate(Emergent = Adult)

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
             color = "white", size = 2.3) +
  scale_fill_custom() +
  labs(y = "P(Detect) in Tetragnathid Spiders",
       x = "P(Detect) in Adult Insects")


# prob_detect_by site 

post_hus_site = posts_concentrations %>% 
  ungroup %>% 
  filter(.draw <= 1000) %>% 
  group_by(type, pfas_type, .draw, site) %>% 
  reframe(hu = mean(hu)) %>% 
  pivot_wider(names_from = type, values_from = hu) %>% 
  mutate(Adult = Emergent)

prob_labels_site_spider_emerge = post_hus_site %>% 
  group_by(pfas_type, site) %>% 
  reframe(Adult = mean(Emergent),
          Tetragnathidae = mean(Tetragnathidae)) %>% 
  mutate(Emergent = Adult)

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
  geom_label(data = prob_labels_site_spider_emerge, 
             aes(label = pfas_type, fill = pfas_type), 
             color = "white", 
             size = 2.3) +
  scale_fill_custom() +
  labs(y = "P(Detect) in Tetragnathid Spiders",
       x = "P(Detect) in Adult Insects")


library(cowplot)
prob_detect_spider_emerge = plot_grid(prob_detect_average_spider_emerge, prob_detect_site_spider_emerge, ncol = 1)

ggsave(prob_detect_spider_emerge, file = "plots/prob_detect_spider_emerge.jpg", width = 6.5, height = 8)

# prob detect summary 

prob_posts_summary = post_hus %>% 
  pivot_longer(cols = c(-pfas_type, -.draw)) %>% 
  group_by(name, pfas_type) %>%
  mutate(value = 1 - value) %>% 
  median_qi(value) %>% 
  make_summary_table(center = "value") %>% 
  select(name, pfas_type, center_interval) %>% 
  pivot_wider(names_from = name, values_from = center_interval)

write_csv(prob_posts_summary, file = "plots/ms_plots_tables/prob_posts_summary.csv")

# fig4_combined -----------------------------------------------------------

fig4a = prob_detect_water_biofilm + labs(subtitle = "a)")
fig4b = prob_detect_biofilm_larvae + labs(subtitle = "b)")
fig4c = prob_detect_average_larv_emerge + labs(subtitle = "c)")
fig4d = prob_detect_average_spider_emerge + labs(subtitle = "d)")

fig4_combined = plot_grid(fig4a, fig4b, fig4c, fig4d, ncol = 2)

ggsave(fig4_combined, file = "plots/ms_plots_tables/fig4_combined.jpg", width = 8, height = 8,
       dpi = 400)


# fig4_concentrations -----------------------------------------------------

# concentration biplots ------------------------------------------------
posts_concentrations = mod_dat %>%
  select(-contains("conc_ppb")) %>%
  distinct(type, taxon, type_taxon, site, pfas_type, order, pfas_order) %>%
  # mutate(site = "new")  %>% 
  add_epred_draws(seed = 20202, hg4_taxon, re_formula = NULL, dpar = T, ndraws = 500) %>%
  mutate(.epred = .epred*unique(hg4_taxon$data2$max_conc_ppb)) %>% 
  mutate(type = case_when(type == "Emergent" ~ "Adult", T ~ type))

post_mus = posts_concentrations %>% 
  group_by(type, site, pfas_type, .draw) %>% # average over taxa
  reframe(.epred = mean(.epred)) %>%
  left_join(pfas_names) %>%
  mutate(pfas_type = reorder(pfas_type, pfas_order)) %>%
  group_by(pfas_type, type, .draw) %>%  # average over sites
  reframe(.epred = mean(.epred)) %>%
  pivot_wider(names_from = type, values_from = .epred)

cons_labels = post_mus %>% 
  group_by(pfas_type) %>% 
  reframe(Water = median(Water),
          Larval = median(Larval),
          Biofilm = median(Biofilm),
          Seston = median(Seston),
          Detritus = median(Detritus))

cons_average = post_mus %>% 
  ggplot(aes(x = Water, y = Larval, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  # facet_wrap(~pfas_type) +
  scale_color_custom() +
  scale_y_log10(breaks = c(0.00001, 0.001, 0.1, 10),
                labels = c("0.00001", "0.001", "0.1", "10"),limits = c(0.00001, 20)) +
  scale_x_log10(breaks = c(0.00001, 0.001, 0.1, 10),
                labels = c("0.00001", "0.001", "0.1", "10"), limits = c(0.00001, 20)) +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = cons_labels, aes(label = pfas_type, fill = pfas_type), 
             color = "white", size = 2.3) +
  scale_fill_custom() +
  labs(x = "Concentration in Water (ppb)",
       y = "Concentration in Larval Insects (ppb)")

# concentration_water_spiders = post_mus %>% 
#   ggplot(aes(x = Water, y = Tetragnathidae, color = pfas_type)) + 
#   geom_point(shape = 20, alpha = 0.5, size = 0.2) +
#   # facet_wrap(~pfas_type) +
#   scale_color_custom() +
#   scale_y_log10(breaks = c(0.00001, 0.001, 0.1, 10),
#                 labels = c("0.00001", "0.001", "0.1", "10"),limits = c(0.00001, 20)) +
#   scale_x_log10(breaks = c(0.00001, 0.001, 0.1, 10),
#                 labels = c("0.00001", "0.001", "0.1", "10"), limits = c(0.00001, 20)) +
#   guides(color = "none",
#          fill = "none") +
#   geom_abline(linetype = "dashed") +
#   geom_label(data = cons_labels, aes(label = pfas_type, 
#                                      y = Spiders,
#                                      fill = pfas_type), 
#              color = "white", size = 2.3) +
#   scale_fill_custom() +
#   labs(x = "Concentration in Water (ppb)",
#        y = "Concentration in Spiders (ppb)")
# 
# 
# ggsave(concentration_water_spiders, file = "plots/ms_plots_tables/concentration_water_spiders.jpg", width = 5, height = 5)


cons_water_biofilm = post_mus %>% 
  ggplot(aes(x = Water, y = Biofilm, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  # facet_wrap(~pfas_type) +
  scale_color_custom() +
  scale_y_log10(breaks = c(0.00001, 0.001, 0.1, 10),
                labels = c("0.00001", "0.001", "0.1", "10"),  limits = c(0.00001, 20)) +
  scale_x_log10(breaks = c(0.00001, 0.001, 0.1, 10),
                labels = c("0.00001", "0.001", "0.1", "10"),  limits = c(0.00001, 20)) +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = cons_labels, aes(label = pfas_type, fill = pfas_type), 
             color = "white", size = 2.3) +
  scale_fill_custom() +
  labs(x = "Concentration in Water (ppb)",
       y = "Concentration in Biofilm (ppb)")

cons_biofilm_larvae = post_mus %>% 
  ggplot(aes(x = Biofilm, y = Larval, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  # facet_wrap(~pfas_type) +
  scale_color_custom() +
  scale_y_log10(breaks = c(0.0001, 0.001,0.01, 0.1, 1, 10, 100),
                labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),  limits = c(0.0001, 100)) +
  scale_x_log10(breaks = c(0.0001, 0.001,0.01, 0.1, 1, 10, 100),
                labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),  limits = c(0.0001, 100)) +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = cons_labels, aes(label = pfas_type, fill = pfas_type), 
             color = "white", size = 2.3) +
  scale_fill_custom() +
  labs(x = "Concentration in Biofilm (ppb)",
       y = "Concentration in Larvae (ppb)")

cons_seston_larvae = post_mus %>% 
  ggplot(aes(x = Seston, y = Larval, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  # facet_wrap(~pfas_type) +
  scale_color_custom() +
  scale_y_log10(breaks = c(0.0001, 0.001,0.01, 0.1, 1, 10, 100),
                labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),  limits = c(0.0001, 100)) +
  scale_x_log10(breaks = c(0.0001, 0.001,0.01, 0.1, 1, 10, 100),
                labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),  limits = c(0.0001, 100)) +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = cons_labels, aes(label = pfas_type, fill = pfas_type), 
             color = "white", size = 2.3) +
  scale_fill_custom() +
  labs(x = "Concentration in Seston (ppb)",
       y = "Concentration in Larvae (ppb)")

cons_detritus_larvae = post_mus %>% 
  ggplot(aes(x = Detritus, y = Larval, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  # facet_wrap(~pfas_type) +
  scale_color_custom() +
  scale_y_log10(breaks = c(0.0001, 0.001,0.01, 0.1, 1, 10, 100),
                labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),  limits = c(0.0001, 100)) +
  scale_x_log10(breaks = c(0.0001, 0.001,0.01, 0.1, 1, 10, 100),
                labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),  limits = c(0.0001, 100)) +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = cons_labels, aes(label = pfas_type, fill = pfas_type), 
             color = "white", size = 2.3) +
  scale_fill_custom() +
  labs(x = "Concentration in Detritus (ppb)",
       y = "Concentration in Larvae (ppb)")

# cons detect by site
post_mus_site = posts_concentrations %>% 
  ungroup %>% 
  filter(.draw <= 1000) %>% 
  group_by(type, pfas_type, site, .draw) %>% 
  # reframe(.epred = median(.epred)) %>% 
  # group_by(type, pfas_type, .draw) %>% 
  reframe(.epred = median(.epred)) %>% 
  pivot_wider(names_from = type, values_from = .epred)

cons_labels_site = post_mus_site %>% 
  group_by(pfas_type, site) %>% 
  reframe(Water = median(Water),
          Larval = median(Larval))

cons_site = post_mus_site %>% 
  ggplot(aes(x = Water, y = Larval, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  facet_wrap(~site) +
  scale_color_custom() +
  scale_y_log10(breaks = c(0.00001, 0.001, 0.1, 10),
                labels = c("0.00001", "0.001", "0.1", "10"),  limits = c(0.00001, 20)) +
  scale_x_log10(breaks = c(0.00001, 0.001, 0.1, 10),
                labels = c("0.00001", "0.001", "0.1", "10"),  limits = c(0.00001, 20)) +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = cons_labels_site, aes(label = pfas_type, fill = pfas_type), 
             color = "white", 
             size = 2.3) +
  scale_fill_custom() +
  labs(x = "Concentration in Water (ppb)",
       y = "Concentration in Larval Insects (ppb)")


concentration_seston_larvae = post_mus %>% 
  ggplot(aes(x = Seston, y = Larval, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  # facet_wrap(~pfas_type) +
  scale_color_custom() +
  scale_y_log10(breaks = c(0.00001, 0.001, 0.1, 10),
                labels = c("0.00001", "0.001", "0.1", "10"),limits = c(0.00001, 20)) +
  scale_x_log10(breaks = c(0.00001, 0.001, 0.1, 10),
                labels = c("0.00001", "0.001", "0.1", "10"), limits = c(0.00001, 20)) +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = cons_labels, aes(label = pfas_type, fill = pfas_type), 
             color = "white", size = 2.3) +
  scale_fill_custom() +
  labs(x = "Concentration in Seston (ppb)",
       y = "Concentration in Larvae (ppb)")


ggsave(concentration_seston_larvae, file = "plots/ms_plots_tables/concentration_seston_larvae.jpg", width = 5, height = 5)


concentration_detritus_larvae = post_mus %>% 
  ggplot(aes(x = Detritus, y = Larval, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  # facet_wrap(~pfas_type) +
  scale_color_custom() +
  scale_y_log10(breaks = c(0.00001, 0.001, 0.1, 10),
                labels = c("0.00001", "0.001", "0.1", "10"),limits = c(0.00001, 20)) +
  scale_x_log10(breaks = c(0.00001, 0.001, 0.1, 10),
                labels = c("0.00001", "0.001", "0.1", "10"), limits = c(0.00001, 20)) +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = cons_labels, aes(label = pfas_type, fill = pfas_type), 
             color = "white", size = 2.3) +
  scale_fill_custom() +
  labs(x = "Concentration in Detritus (ppb)",
       y = "Concentration in Larvae (ppb)")


ggsave(concentration_detritus_larvae, 
       file = "plots/ms_plots_tables/concentration_detritus_larvae.jpg", width = 5, height = 5)


library(cowplot)
cons_water_larvae = plot_grid(cons_average, cons_site, ncol = 1)

ggsave(cons_water_larvae, file = "plots/cons_water_larvae.jpg", width = 6.5, height = 8)



# cons_detect larvae to emerger

cons_labels_larv_emerge = post_mus %>% 
  group_by(pfas_type) %>% 
  reframe(Adult = median(Adult),
          Larval = median(Larval)) %>% 
  mutate(Emergent = Adult)

cons_average_larv_emerge = post_mus %>% 
  ggplot(aes(x = Larval, y = Adult, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  # facet_wrap(~pfas_type) +
  scale_color_custom() +
  scale_y_log10(breaks = c(0.0001, 0.001,0.01, 0.1, 1, 10, 100),
                labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),  limits = c(0.0001, 100)) +
  scale_x_log10(breaks = c(0.0001, 0.001,0.01, 0.1, 1, 10, 100),
                labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),  limits = c(0.0001, 100))+
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = cons_labels_larv_emerge, aes(label = pfas_type, fill = pfas_type), 
             color = "white", size = 2.3) +
  scale_fill_custom() +
  labs(x = "Concentration in Larval Insects (ppb)",
       y = "Concentration in Adult Insects (ppb)")


# cons_by site 

post_mus_site = posts_concentrations %>% 
  ungroup %>% 
  filter(.draw <= 1000) %>% 
  group_by(type, pfas_type, .draw, site) %>% 
  reframe(.epred = median(.epred)) %>% 
  pivot_wider(names_from = type, values_from = .epred)

cons_labels_site_larv_emerge = post_mus_site %>% 
  group_by(pfas_type, site) %>% 
  reframe(Adult = median(Adult),
          Larval = median(Larval)) %>% 
  mutate(Emergent = Adult)

cons_site_larv_emerge = post_mus_site %>% 
  ggplot(aes(x = Larval, y = Adult, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  facet_wrap(~site) +
  scale_color_custom() +
  scale_y_log10(breaks = c(0.0001, 0.001,0.01, 0.1, 1, 10, 100),
                labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),  limits = c(0.0001, 100)) +
  scale_x_log10(breaks = c(0.0001, 0.001,0.01, 0.1, 1, 10, 100),
                labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),  limits = c(0.0001, 100)) +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = cons_labels_site_larv_emerge, aes(label = pfas_type, fill = pfas_type), 
             color = "white", 
             size = 2.3) +
  scale_fill_custom() +
  labs(x = "Concentration in Larval Insects (ppb)",
       y = "Concentration in Adult Insects (ppb)")


library(cowplot)
cons_larv_emerge = plot_grid(cons_average_larv_emerge, cons_site_larv_emerge, ncol = 1)

ggsave(cons_larv_emerge, file = "plots/cons_larv_emerge.jpg", width = 6.5, height = 8)



# cons_detect emerger to spider

cons_labels_spider_emerge = post_mus %>% 
  group_by(pfas_type) %>% 
  reframe(Adult = median(Adult),
          Tetragnathidae = median(Tetragnathidae))%>% 
  mutate(Emergent = Adult)

cons_average_spider_emerge = post_mus %>% 
  ggplot(aes(y = Tetragnathidae, x = Adult, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  # facet_wrap(~pfas_type) +
  scale_color_custom() +
  scale_y_log10(breaks = c(0.0001, 0.001,0.01, 0.1, 1, 10, 100),
                labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),  limits = c(0.0001, 100)) +
  scale_x_log10(breaks = c(0.0001, 0.001,0.01, 0.1, 1, 10, 100),
                labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),  limits = c(0.0001, 100)) +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = cons_labels_spider_emerge, aes(label = pfas_type, fill = pfas_type), 
             color = "white", size = 2.3) +
  scale_fill_custom() +
  labs(y = "Concentration in Tetragnathid Spiders (ppb)",
       x = "Concentration in Adult Insects (ppb)") 


# cons_by site 

post_mus_site = posts_concentrations %>% 
  ungroup %>% 
  filter(.draw <= 1000) %>% 
  group_by(type, pfas_type, .draw, site) %>% 
  reframe(.epred = median(.epred)) %>% 
  pivot_wider(names_from = type, values_from = .epred)

cons_labels_site_spider_emerge = post_mus_site %>% 
  group_by(pfas_type, site) %>% 
  reframe(Adult = median(Adult),
          Tetragnathidae = median(Tetragnathidae)) %>% 
  mutate(Emergent = Adult)

cons_site_spider_emerge = post_mus_site %>% 
  ggplot(aes(y = Tetragnathidae, x = Adult, color = pfas_type)) + 
  geom_point(shape = 20, alpha = 0.5, size = 0.2) +
  facet_wrap(~site) +
  scale_color_custom() +
  scale_y_log10(breaks = c(0.0001, 0.001,0.01, 0.1, 1, 10, 100),
                labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),  limits = c(0.0001, 100)) +
  scale_x_log10(breaks = c(0.0001, 0.001,0.01, 0.1, 1, 10, 100),
                labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),  limits = c(0.0001, 100)) +
  guides(color = "none",
         fill = "none") +
  geom_abline(linetype = "dashed") +
  geom_label(data = cons_labels_site_spider_emerge, 
             aes(label = pfas_type, fill = pfas_type), 
             color = "white", 
             size = 2.3) +
  scale_fill_custom() +
  labs(y = "Concentration in Tetragnathid Spiders (ppb)",
       x = "Concentration in Adult Insects (ppb)")


library(cowplot)
cons_spider_emerge = plot_grid(cons_average_spider_emerge, cons_site_spider_emerge, ncol = 1)

ggsave(cons_spider_emerge, file = "plots/cons_spider_emerge.jpg", width = 6.5, height = 8)

# fig4_combined -----------------------------------------------------------

fig4a_cons = cons_water_biofilm + labs(subtitle = "a)")
fig4b_cons = cons_biofilm_larvae + labs(subtitle = "b)")
fig4c_cons = cons_average_larv_emerge + labs(subtitle = "c)")
fig4d_cons = cons_average_spider_emerge + labs(subtitle = "d)")

fig4_combined_cons = plot_grid(fig4a_cons, fig4b_cons, fig4c_cons, fig4d_cons, ncol = 2)

ggsave(fig4_combined_cons, file = "plots/ms_plots_tables/fig4_combined_cons.jpg", width = 8, height = 8,
       dpi = 400)

# for reviewer response - compare biofilm/larvae to seston/larvae and detritus/larvae
fig4a_cons_rev = cons_biofilm_larvae + labs(subtitle = "a)")
fig4b_cons_rev = cons_seston_larvae + labs(subtitle = "b)")
fig4c_cons_rev = cons_detritus_larvae + labs(subtitle = "c)")

fig4_rev = plot_grid(fig4a_cons_rev, fig4b_cons_rev, fig4c_cons_rev)

ggsave(fig4_rev, file = "plots/temporary/fig4_rev.jpg", width = 6.5,
       height = 7)

# figx_TMF_sum ------------------------------------------------------------

brm_tmf_iso_sum = readRDS("models/brm_tmf_iso_sum.rds")

posts_tmf_iso_sum = brm_tmf_iso_sum$data2 %>% 
  filter(site != "Russell Brook") %>% 
  distinct(site) %>% 
  mutate(center_sum = attributes(brm_tmf_iso_sum$data2$log10_median_sum_s)$`scaled:center`,
         scale_sum = attributes(brm_tmf_iso_sum$data2$log10_median_sum_s)$`scaled:scale`) %>% 
  expand_grid(mean_n15 = seq(min(brm_tmf_iso_sum$data$mean_n15),
                             max(brm_tmf_iso_sum$data$mean_n15), 
                             length.out = 30)) %>% 
  add_epred_draws(brm_tmf_iso_sum, re_formula = NULL) %>% 
  group_by(mean_n15, .draw) %>% 
  reframe(.epred = mean(.epred),
          .epred_log10 = (.epred*scale_sum) + center_sum)


posts_tmf_iso_sum_summary = posts_tmf_iso_sum %>% 
  group_by(mean_n15) %>% 
  median_qi(.epred_log10)

sum_insect_data = brm_tmf_iso_sum$data2 %>% 
  mutate(type = case_when(taxon == "Biofilm" ~ "Biofilm",
                                 taxon == "Spider" ~ "Spiders",
                                 TRUE ~ type)) %>%
  mutate(type = str_replace(type, "Emergent", "Adult")) %>% 
  mutate(type = case_when(type == "Adult" ~ paste(type, "Insects"),
                          type == "Larval" ~ paste(type, "Insects"),
                          TRUE ~ type)) %>% 
  mutate(type = case_when(type == "Tetragnathidae" ~ "Spiders",
                          TRUE ~ type)) %>% 
  mutate(type = fct_relevel(type,"Biofilm", "Larval Insects", "Adult Insects"))

plot_tmf_sum = posts_tmf_iso_sum_summary %>%
  ggplot(aes(x = mean_n15, y = .epred_log10)) +
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper), alpha = 0.3) +
  geom_line() +
  guides(fill = "none") +
  geom_point(data = sum_insect_data,
             aes(y = log10_median_sum,
                 shape = type),
             size = 1) +
  theme(strip.text = element_text(size = 6),
        axis.text.x = element_text(size = 7),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
scale_shape_manual(values = c(16, 15, 0, 17)) +
  labs(x = expression(paste(delta^{15}, "N (centered)")),
       y = "\u221112 PFAS Concentration (log10 tissue ppb)",
       subtitle = "Slope (95%CrI): 0.62 (0.49 to 0.72)",
       shape = "") +
  NULL

# plot_tmf_sum = posts_tmf_iso_sum_summary %>%
#   ggplot(aes(x = mean_n15, y = .epred_log10)) +
#   geom_ribbon(aes(ymin = .lower,
#                   ymax = .upper), alpha = 0.3) +
#   geom_line() +
#   guides(fill = "none") +
#   geom_point(data = brm_tmf_iso_sum$data2,
#              aes(y = log10_median_sum,
#                  shape = taxon),
#              size = 1,
#              position = position_jitter(width = 0.02)) +
#   theme(strip.text = element_text(size = 6),
#         axis.text.x = element_text(size = 7),
#         legend.text = element_text(size = 8),
#         legend.title = element_text(size = 8)) +
#   scale_shape_manual(values = c(15, 16, 17, 18, 19, 7, 10, 11)) +
#   labs(x = expression(paste(delta^{15}, "N (centered)")),
#        y = "\u221112 PFAS Concentration (log10 tissue ppb)",
#        subtitle = "Slope (95%CrI): 0.62 (0.49 to 0.72)",
#        shape = "") +
#   NULL

ggsave(plot_tmf_sum, file = "plots/ms_plots_tables/plot_tmf_sum.jpg", width = 6.5, height = 5)

# tmf sum slopes
mean_sum = attributes(brm_tmf_iso_sum$data2$log10_median_sum_s)$`scaled:center`
sd_sum = attributes(brm_tmf_iso_sum$data2$log10_median_sum_s)$`scaled:scale`

tmf_sum_slopes = brm_tmf_iso_sum$data2 %>%
  distinct(site, center_15n) %>%
  expand_grid(mean_n15 = c(-2,0.4, 1.4)) %>% # 3.4 per mill difference = 1 trophic level. Use -2 and 1.4 to keep within the values of centered n15
  # expand_grid(mean_n15 = seq(0, 1, length.out = 2)) %>%
  add_epred_draws(brm_tmf_iso_sum, re_formula = NULL) %>%
  ungroup %>%
  select(-.row, -.chain, -.iteration) %>% 
  group_by(mean_n15, .draw) %>% 
  reframe(.epred = mean(.epred)) %>% 
  mutate(.epred = (.epred*sd_sum) + mean_sum) %>% 
  pivot_wider(names_from = mean_n15, values_from = .epred) %>%
  # mutate(slope = `1.4`/`-2`) # multiplicative change per trophic level (same as this: 10^(`1.4` - `-2`))
  mutate(slope = `1.4` - `0.4`)

tmf_sum_slope_summary = tmf_sum_slopes %>% 
  median_qi(slope) %>% 
  mutate(pfas_type = "sum_pfas")

sum(tmf_sum_slopes$slope > 0)/nrow(tmf_sum_slopes)

as_tibble(bayes_R2(brm_tmf_iso_sum))

# figx_TMF per pfas----------------------------------------------------------------
pfas_conc_isotopes = readRDS("data/pfas_conc_isotopes_notadjustedformetamorphosis.rds")
brm_tmf_iso_notadjustedformetamorphosis = readRDS("models/brm_tmf_iso_notadjustedformetamorphosis.rds")

mean_conc = attributes(brm_tmf_iso_notadjustedformetamorphosis$data2$log10_median_conc_s)$`scaled:center`
sd_conc = attributes(brm_tmf_iso_notadjustedformetamorphosis$data2$log10_median_conc_s)$`scaled:scale`

brm_tmf_iso_raw_notadjustedformetamorphosis = brm_tmf_iso_notadjustedformetamorphosis$data2

# plot same but culled to PFOA, PFNA, PFUnA, PFOS
posts_tmf_iso_overall_filtered = brm_tmf_iso_notadjustedformetamorphosis$data %>% 
  filter(pfas_type %in% c("PFOA", "PFNA", "PFUnA", "PFOS")) %>% 
  distinct(pfas_type, site) %>% 
  expand_grid(mean_n15 = seq(-1.5, 1, length.out = 30)) %>% 
  add_epred_draws(brm_tmf_iso_notadjustedformetamorphosis, re_formula = NULL) %>% 
  mutate(.epred_log10 = (.epred*sd_conc) + mean_conc) %>% 
  group_by(pfas_type, mean_n15, .draw) %>%
  reframe(.epred_log10 = mean(.epred_log10)) %>% # average over sites
  group_by(pfas_type, mean_n15) %>% 
  median_qi(.epred_log10) %>% 
  mutate(pfas_type = fct_relevel(pfas_type, "PFOA", "PFNA", "PFUnA", "PFOS"))

insect_data = brm_tmf_iso_raw_notadjustedformetamorphosis %>% 
  filter(pfas_type %in% c("PFOA", "PFNA", "PFUnA", "PFOS")) %>% 
  mutate(pfas_type = fct_relevel(pfas_type, "PFOA", "PFNA", "PFUnA", "PFOS")) %>% 
  mutate(taxon_group = case_when(taxon == "Biofilm" ~ "Biofilm", 
                                 taxon == "Spider" ~ "Spider",
                                 TRUE ~ "Insects")) %>% 
  # filter(taxon_group == "Insects") %>% 
  mutate(type = str_replace(type, "Emergent", "Adult")) %>% 
  mutate(type = case_when(type == "Adult" ~ paste(type, "Insects"),
                          type == "Larval" ~ paste(type, "Insects"),
                          TRUE ~ type))


plot_tmf_overall_filtered = posts_tmf_iso_overall_filtered  %>%
  ggplot(aes(x = mean_n15, y = .epred_log10, fill = pfas_type)) +
  geom_point(data = insect_data  %>% 
               mutate(type = case_when(type == "Tetragnathidae" ~ "Spiders",
                                       TRUE ~ type)) %>% 
               mutate(type = fct_relevel(type,"Biofilm", "Larval Insects", "Adult Insects")),
             aes(y = log10_median_conc,
                 shape = type,
                 color = pfas_type
             ),
             size = 1,
             # position = position_jitter(width = 0.02)
  ) +
  geom_ribbon(aes(ymin = .lower, 
                  ymax = .upper), alpha = 0.8) +
  geom_line() +
  scale_fill_custom() + 
  scale_color_custom() +
  facet_wrap2(~ pfas_type, nrow = 1) +
  guides(fill = "none",
         color = "none") +
  # geom_point(data = noninsect_data,
  #            aes(color = pfas_type, 
  #                y = log10_median_conc,
  #                shape = type)) +
  scale_x_continuous(breaks = c(-1.5, 0, 1.5),
                     limits = c(-1.5, 1.5)) +
  scale_shape_manual(values = c(16, 15, 0, 17)) +
  theme(strip.text = element_text(size = 6),
        axis.text.x = element_text(size = 7),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  labs(x = expression(paste(delta^{15}, "N (centered)")),
       y = "PFAS Concentration\n(log10 tissue ppb)",
       shape = "") +
  NULL

ggsave(plot_tmf_overall_filtered, file = "plots/ms_plots_tables/plot_tmf_overall_filtered.jpg",
       width = 8.5, height = 3)

tmf_iso_slopes_bycompound = brm_tmf_iso_notadjustedformetamorphosis$data2 %>%
  distinct(site, pfas_type, center_15n) %>%
  expand_grid(mean_n15 = c(-2, 0.4, 1.4)) %>% # 3.4 per mill difference = 1 trophic level. Use -2 and 1.4 to keep within the values of centered n15
  # expand_grid(mean_n15 = seq(0, 1, length.out = 2)) %>%
  add_epred_draws(brm_tmf_iso_notadjustedformetamorphosis, re_formula = NULL) %>%
  ungroup %>%
  select(-.row, -.chain, -.iteration) %>%
  mutate(.epred = (.epred*sd_conc) + mean_conc) %>% 
  pivot_wider(names_from = mean_n15, values_from = .epred) %>%
  mutate(slope = `1.4` - `0.4`) %>%  # multiplicative change per trophic level (same as this: 10^(`1.4` - `-2`))
  group_by(.draw, pfas_type) %>%
  reframe(slope = mean(slope))

tmf_slopes_bycompound = tmf_iso_slopes_bycompound %>% 
  group_by(pfas_type) %>% 
  median_qi(slope) %>% 
  left_join(pfas_names) %>% 
  bind_rows(tmf_sum_slope_summary %>% mutate(pfas_order = 0)) %>% 
  make_summary_table(, center = "slope") %>% 
  arrange(pfas_order)

write_csv(tmf_slopes_bycompound, file = "plots/ms_plots_tables/tmf_slopes_bycompound.csv")

tmf_iso_slopes_bycompound %>% 
  group_by(pfas_type) %>% 
  reframe(prob_greater = sum(slope > 0)/max(.draw))



# figX_TMF_combined -------------------------------------------------------
library(patchwork)
fig_tmf_combined = plot_tmf_sum/plot_tmf_overall_filtered + plot_layout(heights = c(1, 0.4))

ggsave(fig_tmf_combined, file = "plots/ms_plots_tables/fig_tmf_combine.jpg", width = 6.5, height = 8)

# TMF R2 ----------------------------------------------------------------------
brm_tmf_iso_notadjustedformetamorphosis = readRDS("models/brm_tmf_iso_notadjustedformetamorphosis.rds")

overall_r2 = bayes_R2(brm_tmf_iso_notadjustedformetamorphosis)

preds_brm = brm_tmf_iso_notadjustedformetamorphosis$data %>% 
  filter(pfas_type %in% c("PFOA", "PFNA", "PFUnA", "PFOS")) %>% 
  add_epred_draws(brm_tmf_iso_notadjustedformetamorphosis)

tmf_r2 = preds_brm %>%
  group_by(.draw, pfas_type) %>%
  reframe(
    var_fit = var(.epred),
    var_resid = var(log10_median_conc_s - .epred),
    r2 = var_fit / (var_fit + var_resid)) %>%
  group_by(pfas_type) %>%
  median_qi(r2) %>% 
  add_row(pfas_type = "Overall Model",
          r2 = overall_r2[1],
          .lower = overall_r2[3],
          .upper = overall_r2[4],
          .width = 0.95,
          .point = "median", 
          .interval = "qi")

write_csv(tmf_r2, file = "plots/ms_plots_tables/tmf_r2.csv")

# TMF slopes ----------------------------------------------------------------------
tmf_mod = readRDS("models/brm_tmf_iso_notadjustedformetamorphosis.rds")

tmf_slopes_perpfas = tmf_mod$data2 %>% 
  distinct(site, pfas_type) %>% 
  expand_grid(mean_n15 = c(0, 1)) %>% 
  filter(pfas_type %in% c("PFOA", "PFNA", "PFUnA", "PFOS")) %>% 
  add_epred_draws(tmf_mod) %>% 
  ungroup %>% 
  mutate(center_sum = attributes(tmf_mod$data2$log10_median_conc_s)$`scaled:center`,
         scale_sum = attributes(tmf_mod$data2$log10_median_conc_s)$`scaled:scale`) %>% 
  group_by(site, .draw, pfas_type) %>% 
  mutate(.epred_log10 = (.epred*scale_sum) + center_sum) %>%  
  select(site, pfas_type, mean_n15, .draw, .epred_log10) %>% 
  pivot_wider(names_from = mean_n15, values_from = .epred_log10) %>% 
  mutate(slope = `1` - `0`) %>% 
  group_by(pfas_type, .draw) %>% 
  reframe(slope = mean(slope)) 

tmf_slopes_perpfas_summary = tmf_slopes_perpfas %>% 
  group_by(pfas_type) %>% 
  median_qi(slope)

tmf_slopes_perpfas %>% 
  group_by(pfas_type) %>% 
  reframe(prob_positive = sum(slope > 0)/max(.draw))

tmf_slope_overall = tmf_slopes_perpfas %>% 
  group_by(.draw) %>% 
  reframe(slope = mean(slope)) %>% 
  median_qi(slope) %>% 
  mutate(pfas_type = "Overall Model")

write_csv(bind_rows(tmf_slopes_perpfas_summary,
                    tmf_slope_overall),
          file = "plots/ms_plots_tables/tmf_slopes.csv")

# Other Tables ------------------------------------------------------------------

insect_mass_per_taxon = readRDS(file = "posteriors/insect_mass_averaged_over_sites.rds") %>% 
  clean_names() %>% 
  mutate(units = "mass_per_insect_grams") %>% 
  select(-type) %>% 
  separate(type_taxon, into = c("stage", "taxon")) %>% 
  make_summary_table(center = "mean_gww", lower = "lower", upper = "upper", digits = 2) %>% 
  unite("stage_units", c(stage,units)) %>% 
  select(stage_units, taxon, center_interval) %>% 
  pivot_wider(names_from = stage_units, values_from = center_interval) %>% 
  clean_names()

sum_mass_mef_ratios = read_csv(file = "plots/ms_plots_tables/sum_partitioning_taxon.csv") %>% 
  left_join(insect_mass_per_taxon) %>% 
  left_join(read_csv(file = "plots/ms_plots_tables/ratio_in_sum_body_burden.csv") %>% 
              select(-metric) %>% rename(burden_ratio_adult_to_larval = center_interval)) %>% 
  left_join(read_csv(file = "plots/ms_plots_tables/change_in_sum_body_burden.csv")%>% 
              select(-metric) %>% rename(burden_diff_adult_to_larval = center_interval)) %>% 
  select(taxon, starts_with("larval"), starts_with("adult"), mef, starts_with("burden_diff"), everything())

write_csv(sum_mass_mef_ratios, file = "plots/ms_plots_tables/sum_mass_mef_ratios.csv")


# figx_isotopes -----------------------------------------------------------

brm_isotopes = readRDS(file = "models/brm_isotopes.rds")

library(tidybayes)
library(ggthemes)

iso_posts_raw = brm_isotopes$data2 %>% 
  distinct(site, sample_type, mean_15n, mean_13c) %>% 
  add_epred_draws(brm_isotopes) %>% 
  ungroup %>% 
  mutate(site = fct_relevel(site, "Burr Pond Brook", "Hop Brook", "Ratlum Brook", "Pequabuck River")) %>% 
  select(site, sample_type, .category, .draw, .epred, mean_15n, mean_13c) 
  
iso_posts = iso_posts_raw %>% 
  pivot_wider(names_from = .category, values_from = .epred) %>% 
  group_by(site, sample_type) %>% 
  reframe(d13c_median = median(d13c + mean_13c),
          d13c_lower = quantile(d13c + mean_13c, probs = 0.025),
          d13c_upper = quantile(d13c + mean_13c, probs = 0.975),
          d15n_median = median(d15n + mean_15n),
          d15n_lower = quantile(d15n + mean_15n, probs = 0.025),
          d15n_upper = quantile(d15n + mean_15n, probs = 0.975),
          d13c_median_c = median(d13c),
          d13c_lower_c = quantile(d13c, probs = 0.025),
          d13c_upper_c = quantile(d13c, probs = 0.975),
          d15n_median_c = median(d15n),
          d15n_lower_c = quantile(d15n, probs = 0.025),
          d15n_upper_c = quantile(d15n, probs = 0.975)) %>% 
  mutate(sample_type = str_replace(sample_type, "Emergent", "Adult")) %>% 
  mutate(sample_type = fct_relevel(sample_type, "Spider", "Adult Trichoptera", "Adult Plecoptera", "Adult Odonata", "Adult Ephemeroptera",
                                   "Adult Diptera"))

isotope_biplot = iso_posts %>% 
  ggplot(aes(x = d13c_median, y = d15n_median, color = sample_type, group = sample_type)) + 
  geom_point(aes(color = sample_type)) + 
  geom_errorbar(aes(ymin = d15n_lower, ymax = d15n_upper)) + 
  geom_errorbarh(aes(xmin = d13c_lower, xmax = d13c_upper)) +
  facet_wrap(~site, nrow = 1) + 
  guides(color = guide_legend(override.aes = list(shape = 19, alpha = 1))) +
  labs(color = "Organism",
       x = expression(paste(delta^{13}, "C")),
       y = expression(paste(delta^{15}, "N"))) +
  theme(legend.text = element_text(size = 8),
        legend.position = "top",
        legend.title = element_blank()) +
  scale_color_manual(values = c("black", "#B3DE69", "#FCCDE5","#8DA0CB", "#FDB462", "#8DD3C7")) + # mimic colors in Figure 3
  # geom_point(data = brm_isotopes$data) +
  NULL

ggsave(isotope_biplot, file = "plots/ms_plots_tables/isotope_biplot.jpg", dpi = 400, width = 6.5, height = 5)


isotope_biplot_singlepanel_each_site = iso_posts %>% 
  ggplot(aes(x = d13c_median_c, y = d15n_median_c, color = sample_type, group = sample_type)) + 
  geom_point(aes(color = sample_type, shape = site)) + 
  geom_errorbar(aes(ymin = d15n_lower_c, ymax = d15n_upper_c)) + 
  geom_errorbarh(aes(xmin = d13c_lower_c, xmax = d13c_upper_c)) +
  # facet_wrap(~site, nrow = 1) + 
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(color = "Organism",
       x = expression(paste(delta^{13}, "C")),
       y = expression(paste(delta^{15}, "N"))) +
  theme(legend.text = element_text(size = 8),
        legend.title = element_blank()) +
  scale_color_manual(values = c("black", "#B3DE69", "#FCCDE5","#8DA0CB", "#FDB462", "#8DD3C7")) + # mimic colors in Figure 3
  # geom_point(data = brm_isotopes$data) +
  NULL

ggsave(isotope_biplot_singlepanel_each_site,
       file = "plots/ms_plots_tables/isotope_biplot_singlepanel_each_site.jpg",
       dpi = 400, width = 6.5, height = 6)


isotope_biplot_summary = iso_posts_raw %>% 
  mutate(.epred = case_when(.category == "d13c" ~ .epred + mean_13c,
                            TRUE ~ .epred + mean_15n)) %>% 
  group_by(site, .category, sample_type) %>% 
  median_qi(.epred) %>% 
  make_summary_table() %>% 
  select(site, .category, center_interval, sample_type) %>% 
  pivot_wider(names_from = sample_type, values_from = center_interval)

write_csv(isotope_biplot_summary, file = "plots/ms_plots_tables/isotope_biplot_summary.csv")


# iso biplot averaged over sites
brm_isotopes_notfiltered = readRDS("models/brm_isotopes_notfiltered.rds")

iso_posts_siteaverage = brm_isotopes_notfiltered$data2 %>% 
  distinct(site, sample_type, mean_15n, mean_13c) %>% 
  add_epred_draws(brm_isotopes_notfiltered) %>% 
  ungroup %>% 
  mutate(site = fct_relevel(site, "Burr Pond Brook", "Hop Brook", "Ratlum Brook", "Pequabuck River")) %>% 
  select(sample_type, .category, .draw, .epred) %>% 
  group_by(sample_type, .category, .draw) %>% 
  reframe(.epred = mean(.epred)) %>% # average over sites
  pivot_wider(names_from = .category, values_from = .epred) %>% 
  group_by(sample_type) %>% 
  reframe(d13c_median_c = median(d13c),
          d13c_lower_c = quantile(d13c, probs = 0.025),
          d13c_upper_c = quantile(d13c, probs = 0.975),
          d15n_median_c = median(d15n),
          d15n_lower_c = quantile(d15n, probs = 0.025),
          d15n_upper_c = quantile(d15n, probs = 0.975)) %>% 
  mutate(sample_type = str_replace(sample_type, "Emergent", "Adult")) %>% 
  mutate(sample_type = fct_relevel(sample_type, "Spider", "Adult Trichoptera", "Adult Plecoptera", "Adult Odonata", "Adult Ephemeroptera",
                                   "Adult Diptera"))


isotope_biplot_siteaveraged = iso_posts_siteaverage %>% 
  ggplot(aes(x = d13c_median_c, y = d15n_median_c, color = sample_type, group = sample_type)) + 
  geom_point(aes(color = sample_type)) + 
  geom_errorbar(aes(ymin = d15n_lower_c, ymax = d15n_upper_c)) + 
  geom_errorbarh(aes(xmin = d13c_lower_c, xmax = d13c_upper_c)) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(color = "Organism",
       x = expression(paste(delta^{13}, "C")),
       y = expression(paste(delta^{15}, "N"))) +
  theme(legend.text = element_text(size = 8),
        legend.title = element_blank()) +
  scale_color_manual(values = c("black", "#B3DE69", "#FCCDE5","#8DA0CB", "#FDB462", "#8DD3C7")) + # mimic colors in Figure 3
  NULL


ggsave(isotope_biplot_siteaveraged,
       file = "plots/ms_plots_tables/isotope_biplot_siteaveraged.jpg",
       dpi = 400, width = 6.5, height = 6)


isotope_biplot_summary_siteaveraged = brm_isotopes$data2 %>% 
  distinct(site, sample_type, mean_15n, mean_13c) %>% 
  add_epred_draws(brm_isotopes) %>% 
  ungroup %>% 
  mutate(site = fct_relevel(site, "Burr Pond Brook", "Hop Brook", "Ratlum Brook", "Pequabuck River")) %>% 
  select(sample_type, .category, .draw, .epred) %>% 
  group_by(sample_type, .category, .draw) %>% 
  reframe(.epred = mean(.epred)) %>% 
  group_by(sample_type, .category) %>% 
  median_qi(.epred) %>% 
  make_summary_table() %>% 
  select(sample_type, .category, center_interval) %>% 
  pivot_wider(names_from = .category, values_from = center_interval)

write_csv(isotope_biplot_summary_siteaveraged, file = "plots/ms_plots_tables/isotope_biplot_summary_siteaveraged_centered.csv")  


brm_isotopes_notfiltered = readRDS("models/brm_isotopes_notfiltered.rds")

isotope_biplot_summary_siteaveraged_larvae = brm_isotopes_notfiltered$data2 %>% 
  distinct(site, sample_type, mean_15n, mean_13c) %>% 
  add_epred_draws(brm_isotopes_notfiltered) %>% 
  ungroup %>% 
  mutate(site = fct_relevel(site, "Burr Pond Brook", "Hop Brook", "Ratlum Brook", "Pequabuck River")) %>% 
  select(sample_type, .category, .draw, .epred) %>% 
  group_by(sample_type, .category, .draw) %>% 
  reframe(.epred = mean(.epred)) %>% 
  group_by(sample_type, .category) %>% 
  median_qi(.epred) %>% 
  make_summary_table() %>% 
  select(sample_type, .category, center_interval) %>% 
  pivot_wider(names_from = .category, values_from = center_interval) %>% 
  mutate(sample_type = case_when(sample_type == "Biofilm" ~ "Larval Biofilm",
                                 sample_type == "Spider" ~ "Larval Spider",
                                 TRUE ~ sample_type)) %>% 
  separate(sample_type, into = c("stage", "taxon")) %>% 
  pivot_longer(starts_with("d1")) %>% 
  unite('stage_name', c(stage, name)) %>% 
  pivot_wider(names_from = stage_name, values_from = value) %>% 
  filter(taxon != "Megaloptera") 

write_csv(isotope_biplot_summary_siteaveraged_larvae, 
          file = "plots/ms_plots_tables/isotope_biplot_summary_siteaveraged_centered_larvae.csv")  

isotope_biplot_summary_siteaveraged_uncentered = brm_isotopes$data2 %>% 
  distinct(site, sample_type, mean_15n, mean_13c) %>% 
  add_epred_draws(brm_isotopes) %>% 
  ungroup %>% 
  mutate(.epred = case_when(.category == "d13c" ~ .epred + mean_13c,
                            TRUE ~ .epred + mean_15n)) %>% 
  select(sample_type, .category, .draw, .epred) %>% 
  group_by(sample_type, .category, .draw) %>% 
  reframe(.epred = mean(.epred)) %>% 
  group_by(sample_type, .category) %>% 
  median_qi(.epred) %>% 
  make_summary_table() %>% 
  select(sample_type, .category, center_interval) %>% 
  pivot_wider(names_from = .category, values_from = center_interval)

write_csv(isotope_biplot_summary_siteaveraged_uncentered, 
          file = "plots/ms_plots_tables/isotope_biplot_summary_siteaveraged_uncentered.csv")  
