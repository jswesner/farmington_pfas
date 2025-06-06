library(tidybayes)
library(tidyverse)
library(brms)
library(janitor)
library(viridis)
library(ggh4x)
library(readxl)
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

# load models and data
insect_mass = read_excel("data/Farmington_InsectMasses.xlsx") %>% clean_names() %>% 
  mutate(units = "grams")

mean_insect_mass = insect_mass %>% 
  group_by(site, order, life_stage) %>% 
  reframe(mean_gdw = mean(0.2*composite_sample_mass_ww, na.rm = T), # 0.2 converts wet to dry mass
          sd_gdw = sd(0.2*composite_sample_mass_ww, na.rm = T)) %>% 
  rename(taxon = order) %>% 
  mutate(type = case_when(life_stage == "Adult" ~ "Emergent", T ~ life_stage),
         type_taxon = paste0(type, "_", taxon)) %>% 
  mutate(site = case_when(site == "Pequabuck Brook" ~ "Pequabuck River", TRUE ~ site))

hg4_taxon = readRDS(file = "models/hg4_taxon.rds")
merged_d2 = hg4_taxon$data2$merged_d2

max_conc = unique(merged_d2$max_conc)

# plot insect masses
mean_insect_mass %>% 
  ggplot(aes(x = taxon, color = life_stage,
             y = mean_gdw)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_gdw - sd_gdw, ymax = mean_gdw + sd_gdw)) +
  facet_wrap(~site)

# estimate body burden ----------------------------------------------------

insect_posts = merged_d2 %>% 
  distinct(site, type_taxon, taxon, type, max_conc, pfas_type) %>% 
  left_join(mean_insect_mass) %>% 
  mutate(type = as.factor(type),
         type = fct_relevel(type, "Larval", "Emergent")) %>% 
  add_epred_draws(hg4_taxon, ndraws = 500) %>% 
  mutate(.epred = .epred*max_conc) 

diff_summary = insect_posts %>% 
  filter(type %in% c("Emergent", "Larval")) %>% 
  mutate(ng_per_bug = (.epred/mean_gdw)) %>% 
  ungroup %>%
  select(site, taxon, type, ng_per_bug, .draw, pfas_type) %>%
  pivot_wider(names_from = type, values_from = ng_per_bug) %>% 
  group_by(pfas_type) %>% 
  mutate(diff = Emergent - Larval) %>% 
  group_by(pfas_type) %>% 
  mutate(diff = diff) %>% 
  filter(!is.na(diff)) %>% 
  # filter(pfas_type == "PFOS") %>%
  group_by(taxon, .draw, pfas_type) %>% 
  reframe(diff = mean(diff)) %>% 
  group_by(taxon, pfas_type) %>% 
  median_qi(diff, .width = 0.5) %>% 
  mutate(taxon = as.factor(taxon),
         taxon = fct_relevel(taxon, "Plecoptera", "Odonata", "Ephemeroptera", "Trichoptera", "Diptera"))


change_in_body_burden = diff_summary %>% 
  ggplot(aes(x = diff, y = pfas_type)) + 
  geom_point(aes(color = taxon), size = 0.8) + 
  geom_errorbarh(aes(xmin = .lower, xmax = .upper,
                     color = taxon), height = 0,
                 linewidth = 0.2) +
  # facet_grid2( pfas_type ~ .) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  xlim(-2000, 5000) +
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

ggsave(change_in_body_burden, file = "plots/change_in_body_burden.jpg", 
       width = 6.5, height = 4)



# estimate isotope and tmf ------------------------------------------------
iso_posts_notadjustedformetamorphosis = readRDS(file = "posteriors/iso_posts_notadjustedformetamorphosis.rds") 

pfas_posts = merged_d2 %>% 
  distinct(site, type_taxon, taxon, type, max_conc, pfas_type) %>% 
  filter(type != "Sediment") %>% 
  filter(type != "Detritus") %>%
  filter(type != "Water") %>% 
  filter(type != "Seston") %>% 
  mutate(taxon = case_when(is.na(taxon) ~ type,TRUE ~ taxon)) %>% 
  add_epred_draws(hg4_taxon, re_formula = NULL) %>% 
  mutate(.epred = .epred*max_conc) 

iso_post_summaries_notadjustedformetamorphosis = iso_posts_notadjustedformetamorphosis %>% 
  group_by(taxon, site, center_15n) %>% 
  reframe(mean_n15 = mean(.epred),
          median_n15 = median(.epred),
          sd_n15 = sd(.epred),
          lower_n15 = quantile(.epred, probs = 0.025),
          upper_n15 = quantile(.epred, probs = 0.975)) %>% 
  filter(site != "Russell Brook") # No spiders and only 2 dots for insects

pfas_conc_isotopes_notadjustedformetamorphosis = pfas_posts %>% 
  filter(site != "Russell Brook") %>% # only 2 data points. Can fit regression here and model below was having trouble so I removed.
  group_by(site, taxon, type, pfas_type)  %>% 
  reframe(mean_conc = mean(.epred),
          median_conc = median(.epred),
          sd_conc = sd(.epred),
          lower_conc = quantile(.epred, probs = 0.025),
          upper_conc = quantile(.epred, probs = 0.975)) %>% 
  left_join(iso_post_summaries_notadjustedformetamorphosis, relationship = "many-to-many") %>% 
  mutate(log10_mean_conc = log10(mean_conc),
         log10_sd_conc = sd(log10_mean_conc),
         log10_mean_conc_s = scale(log10_mean_conc),
         log10_median_conc = log10(median_conc),
         log10_median_conc_s = scale(log10_median_conc)) 

pfas_conc_isotopes_notadjustedformetamorphosis %>%
  # mutate(range = upper_conc - lower_conc) %>% arrange(-range) %>% select(range)
  # filter(site == "Burr Pond Brook") %>% 
  ggplot(aes(x = mean_n15, color = pfas_type, y = log10_median_conc)) + 
  geom_pointrange(aes(ymin = log10_median_conc + log10_sd_conc,
                      ymax = log10_median_conc - log10_sd_conc)) +
  # geom_text(aes(label = taxon)) +
  facet_grid2(site ~ pfas_type, scales = "free_y") +
  scale_color_custom() +
  # scale_y_log10() +
  # scale_x_log10() +
  guides(fill = "none",
         color = "none") +
  # geom_smooth(method = lm, se = T) +
  NULL

# fit model

brm_tmf_iso_notadjustedformetamorphosis = brm(log10_median_conc_s ~ mean_n15 + (1 + mean_n15|site) + (1 + mean_n15|pfas_type),
                  data = pfas_conc_isotopes_notadjustedformetamorphosis,
                  family = gaussian(),
                  prior = c(prior(normal(0, 1), class = Intercept),
                             prior(normal(0, 1), class = b),
                             prior(exponential(2), class = sd)),
                  iter = 2000, chains = 4)
# 
# saveRDS(brm_tmf_iso_notadjustedformetamorphosis, file = "models/brm_tmf_iso_notadjustedformetamorphosis.rds")


# # check model
pp_check(brm_tmf_iso_notadjustedformetamorphosis, type = "dens_overlay_grouped", group = "pfas_type")

brm_tmf_iso_notadjustedformetamorphosis = readRDS(file = "models/brm_tmf_iso_notadjustedformetamorphosis.rds")

posts_tmf_iso_notadjustedformetamorphosis = brm_tmf_iso_notadjustedformetamorphosis$data %>% 
  distinct(site, pfas_type) %>% 
  expand_grid(mean_n15 = seq(-2, 1, length.out = 30)) %>% 
  add_epred_draws(brm_tmf_iso_notadjustedformetamorphosis, re_formula = NULL)

mean_conc = attributes(pfas_conc_isotopes_notadjustedformetamorphosis$log10_median_conc_s)$`scaled:center`
sd_conc = attributes(pfas_conc_isotopes_notadjustedformetamorphosis$log10_median_conc_s)$`scaled:scale`

posts_tmf_iso_summary_notadjustedformetamorphosis = posts_tmf_iso_notadjustedformetamorphosis %>% 
  group_by(mean_n15, site, pfas_type) %>% 
  mutate(.epred_log10 = (.epred*sd_conc) + mean_conc) %>% 
  median_qi(.epred_log10)
  
brm_tmf_iso_raw_notadjustedformetamorphosis = brm_tmf_iso_notadjustedformetamorphosis$data %>%
  left_join(pfas_conc_isotopes_notadjustedformetamorphosis %>% ungroup %>% distinct(site, pfas_type)) %>% 
  mutate(.epred_log10 = (log10_median_conc_s*sd_conc) + mean_conc)

plot_tmf_isotopes_notadjustedformetamorphosis = posts_tmf_iso_summary_notadjustedformetamorphosis %>% 
  # filter(pfas_type == "6:2FTS") %>% 
  ggplot(aes(x = mean_n15, y = .epred_log10, fill = pfas_type)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line() +
  scale_fill_custom() + 
  scale_color_custom() +
  facet_grid2(site ~ pfas_type, scales = "free") +
  geom_point(data = brm_tmf_iso_raw_notadjustedformetamorphosis,
             aes(color = pfas_type),
             shape = 1,
             size = 0.5) +
  guides(fill = "none",
         color = "none") +
  scale_x_continuous(breaks = c(-1.5, 0, 1.5)) +
  theme(strip.text = element_text(size = 6),
        axis.text.x = element_text(size = 7)) +
  labs(x = expression(paste(delta^{15}, "N (centered)")),
       y = "log10(ppb)") +
  NULL

ggsave(plot_tmf_isotopes_notadjustedformetamorphosis, file = "plots/plot_tmf_isotopes_notadjustedformetamorphosis.jpg",
       width = 10, height = 4)


# slopes ------------------------------------------------------------------

plot_tmf_isotopes2_notadjustedformetamorphosis = posts_tmf_iso_summary_notadjustedformetamorphosis %>% 
  # filter(pfas_type == "6:2FTS") %>% 
  ggplot(aes(x = mean_n15, y = 10^.epred_log10, fill = pfas_type)) +
  geom_ribbon(aes(ymin = 10^.lower, 
                  ymax = 10^.upper), alpha = 0.3) +
  geom_line() +
  scale_fill_custom() + 
  scale_color_custom() +
  facet_grid2(site ~ pfas_type, scales = "free") +
  guides(fill = "none",
         color = "none") +
  geom_point(data = brm_tmf_iso_raw_notadjustedformetamorphosis,
             aes(color = pfas_type),
             shape = 1,
             size = 0.5) +
  scale_x_continuous(breaks = c(-1.5, 0, 1.5)) +
  theme(strip.text = element_text(size = 6),
        axis.text.x = element_text(size = 7)) +
  labs(x = expression(paste(delta^{15}, "N (centered)")),
       y = "PFAS Concentration (tissue ppb)") +
  NULL

ggsave(plot_tmf_isotopes2_notadjustedformetamorphosis, file = "plots/plot_tmf_isotopes2_notadjustedformetamorphosis.jpg",
       width = 10, height = 4)

tmf_iso_slopes_notadjustedformetamorphosis = brm_tmf_iso_notadjustedformetamorphosis$data %>% 
  distinct(site, pfas_type) %>% 
  # expand_grid(mean_n15 = seq(-2, 1.4, length.out = 2)) %>% # 3.4 per mill difference = 1 trophic level. Use -2 and 1.4 to keep within the values of centered n15
  expand_grid(mean_n15 = seq(0, 1, length.out = 2)) %>% 
  add_epred_draws(brm_tmf_iso_notadjustedformetamorphosis, re_formula = NULL) %>%
  ungroup %>% 
  select(-.row, -.chain, -.iteration) %>% 
  pivot_wider(names_from = mean_n15, values_from = .epred) %>% 
  mutate(slope = `1` - `0`,
         slope = 10^slope)

tmf_slope_table_notadjustedformetamorphosis = tmf_iso_slopes_notadjustedformetamorphosis %>% 
  group_by(site, pfas_type) %>% 
  median_qi(slope)

write_csv(tmf_slope_table_notadjustedformetamorphosis, file = "tables/tmf_slope_table_notadjustedformetamorphosis.csv")

library(ggridges)
tmf_slopes_notadjustedformetamorphosis = tmf_iso_slopes_notadjustedformetamorphosis %>% 
  ggplot(aes(x = slope, y = pfas_type, fill = site)) + 
  # stat_halfeye(aes(fill = site, color = site)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, scale = 1, 
                      alpha = 0.7) +
  scale_fill_brewer() +
  # scale_x_log10() +
  # ggthemes::scale_fill_colorblind() +
  # ggthemes::scale_color_colorblind() +
  labs(x = expression(paste("TMF (Fold increase in PFAS per unit increase in ", delta^{15}, "N)")),
       y = "",
       fill = "",
       color = "") +
  coord_cartesian(xlim = c(NA, 5))


ggsave(tmf_slopes_notadjustedformetamorphosis, file = "plots/tmf_slopes_notadjustedformetamorphosis.jpg", width = 6.5, height = 9)


# tmf slopes not site specific --------------------------------------------
tmf_iso_slopes_overall = brm_tmf_iso_notadjustedformetamorphosis$data %>% 
  distinct(pfas_type) %>% 
  # expand_grid(mean_n15 = seq(-2, 1.4, length.out = 2)) %>% # 3.4 per mill difference = 1 trophic level. Use -2 and 1.4 to keep within the values of centered n15
  expand_grid(mean_n15 = seq(0, 1, length.out = 2)) %>% 
  add_epred_draws(brm_tmf_iso_notadjustedformetamorphosis, 
                  re_formula = ~ (1 + mean_n15|pfas_type)) %>%
  ungroup %>% 
  select(-.row, -.chain, -.iteration) %>% 
  pivot_wider(names_from = mean_n15, values_from = .epred) %>% 
  mutate(slope = `1` - `0`,
         slope = 10^slope)


tmf_slopes_overall = tmf_iso_slopes_overall %>% 
  group_by(pfas_type) %>% 
  median_qi(slope)

write_csv(tmf_slopes_overall, file = "tables/tmf_slopes_overall.csv")

mean_conc = attributes(pfas_conc_isotopes_notadjustedformetamorphosis$log10_median_conc_s)$`scaled:center`
sd_conc = attributes(pfas_conc_isotopes_notadjustedformetamorphosis$log10_median_conc_s)$`scaled:scale`

posts_tmf_iso_overall = brm_tmf_iso_notadjustedformetamorphosis$data %>% 
  distinct(pfas_type) %>% 
  expand_grid(mean_n15 = seq(-2, 1, length.out = 30)) %>% 
  add_epred_draws(brm_tmf_iso_notadjustedformetamorphosis, re_formula = ~ (1 + mean_n15|pfas_type)) %>% 
  mutate(.epred_log10 = (.epred*sd_conc) + mean_conc) %>% 
  group_by(pfas_type, mean_n15) %>% 
  median_qi(.epred_log10)


plot_tmf_overall = posts_tmf_iso_overall %>% 
  # filter(pfas_type == "6:2FTS") %>% 
  ggplot(aes(x = mean_n15, y = .epred_log10, fill = pfas_type)) +
  geom_ribbon(aes(ymin = .lower, 
                  ymax = .upper), alpha = 0.3) +
  geom_line() +
  scale_fill_custom() + 
  scale_color_custom() +
  facet_wrap2(~ pfas_type, nrow = 1) +
  guides(fill = "none",
         color = "none") +
  geom_point(data = brm_tmf_iso_raw_notadjustedformetamorphosis,
             aes(color = pfas_type),
             shape = 1,
             size = 0.5) +
  scale_x_continuous(breaks = c(-1.5, 0, 1.5)) +
  theme(strip.text = element_text(size = 6),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 7)) +
  labs(x = expression(paste(delta^{15}, "N (centered)")),
       y = "PFAS Concentration\n(tissue ppb)") +
  NULL


ggsave(plot_tmf_overall, file = "plots/plot_tmf_overall.jpg",
       width = 10, height = 2)

# errors in variables ---------------------------------------
# this was not fitting efficiently - lots of high rhats even with tight priors.
# leave it for now
# insect_conc_isotopes %>%
#   filter(pfas_type == "PFOS") %>% 
#   ggplot(aes(x = mean_n15, color = pfas_type, y = log10_mean_conc)) + 
#   geom_pointrange(aes(xmin = mean_n15 - sd_n15, 
#                       xmax = mean_n15 + sd_n15,
#                       ymin = log10_mean_conc + log10_mean_conc_sd,
#                       ymax = log10_mean_conc - log10_mean_conc_sd)) +
#   # geom_text(aes(label = taxon)) +
#   facet_grid(site ~ pfas_type,
#              scales = "free_y") +
#   scale_color_custom() +
#   # scale_y_log10() +
#   # scale_x_log10() +
#   guides(fill = "none",
#          color = "none") +
#   geom_smooth(method = lm, se = T)
# 
# # errors in variables
# brm_tmf_iso_mi = brm(log10_median_conc | mi(log10_sd_conc) ~ mean_n15 +
#                        (1 + mean_n15|site) + (1 + mean_n15|pfas_type),
#                      data = pfas_conc_isotopes,
#                      family = gaussian(),
#                      prior = c(prior(normal(0, 1), class = "b"),
#                                prior(normal(5, 2), class = "Intercept"),
#                                prior(exponential(3), class = "sd"),
#                                prior(exponential(3), class = "sigma")),
#                      iter = 200, chains = 1, cores = 4)
# 
# saveRDS(brm_tmf_iso_mi, file = "models/brm_tmf_iso_mi.rds")
# 
# brm_tmf_iso_mi = update(brm_tmf_iso_mi, chains = 1, iter = 200, newdata = pfas_conc_isotopes)
# 
# # 
# 
# mi_data_grid_list = brm_tmf_iso_mi$data %>% 
#   group_by(site, pfas_type) %>% 
#   group_split()
# 
# sim_mi_list = list()
# 
# for(i in 1:length(mi_data_grid_list)){
#   sim_mi_list[[i]] = tibble(mean_n15 = seq(min(mi_data_grid_list[[i]]$mean_n15),
#                          max(mi_data_grid_list[[i]]$mean_n15),
#                          length.out = 30),
#                          site = unique(mi_data_grid_list[[i]]$site),
#                          pfas_type = unique(mi_data_grid_list[[i]]$pfas_type))
# }
# 
# 
# 
# posts_tmf_iso_mi = bind_rows(sim_mi_list) %>%
#   mutate(sd_n15 = mean(brm_tmf_iso_mi$data$sd_n15),
#          log10_mean_conc_sd = mean(brm_tmf_iso_mi$data$log10_mean_conc_sd)) %>% 
#   add_epred_draws(brm_tmf_iso_mi, re_formula = NULL)
# 
# posts_tmf_iso_summary_mi = posts_tmf_iso_mi %>% 
#   group_by(mean_n15, site, pfas_type) %>%  
#   median_qi(.epred) 
# 
# brm_tmf_iso_raw_mi = brm_tmf_iso_mi$data %>% 
#   filter(site != "Russell Brook") %>% 
#   left_join(insect_conc_isotopes %>% 
#               ungroup %>% 
#               distinct(n15_scale, n15_center, site, pfas_type)) 
# 
# 
# posts_tmf_iso_summary_mi %>% 
#   left_join(insect_conc_isotopes %>% 
#               ungroup %>% 
#               distinct(n15_scale, n15_center, site, pfas_type)) %>% 
#   ggplot(aes(x = (mean_n15 - n15_center)/n15_scale, y = .epred, fill = pfas_type)) +
#   geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
#   geom_line() +
#   scale_fill_custom() + 
#   facet_grid2(site ~ pfas_type) +
#   geom_point(data = brm_tmf_iso_raw_mi,
#              aes(color = pfas_type,
#                  y = log10_mean_conc),
#              shape = 1,
#              size = 0.5) +
#   theme(strip.text = element_text(size = 8)) +
#   labs(x = "\u2030N15") +
#   NULL
# 
