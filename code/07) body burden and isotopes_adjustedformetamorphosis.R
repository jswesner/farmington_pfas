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
  rename(taxon = order)

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

insect_posts = merged_d2  
  distinct(site, type_taxon, taxon, type, max_conc, pfas_type) %>% 
  left_join(mean_insect_mass, relationship = "many-to-many") %>% 
  add_epred_draws(hg4_taxon, ndraws = 500) %>% 
  mutate(.epred = .epred*max_conc) 

insect_posts %>% 
  mutate(ng_per_bug = .epred/mean_gdw) %>% 
  group_by(site, taxon, type_taxon, type, pfas_type) %>% 
  median_qi(ng_per_bug) %>% 
  ggplot(aes(x = type_taxon, y = ng_per_bug, ymin = .lower, ymax = .upper)) +
  geom_pointrange() +
  facet_wrap(~pfas_type) +
  # scale_y_log10() +
  NULL


insect_posts %>% 
  mutate(ng_per_bug = .epred/mean_gdw) %>% 
  ungroup %>%
  # filter(taxon == "Diptera") %>%
  # distinct(type) %>% 
  # filter(pfas_type == "PFOS") %>% 
  select(site, taxon, type, .epred, ng_per_bug) %>% 
  ggplot(aes(x = taxon, color = type, y = ng_per_bug)) + 
  # geom_point() +
  stat_pointinterval() +
  # geom_boxplot(aes(group = type)) +
  facet_wrap(~site) +
  # stat_pointinterval() +
  scale_y_log10() +
  NULL



# estimate isotope and tmf ------------------------------------------------
iso_posts = readRDS(file = "posteriors/iso_posts.rds") 

pfas_posts = merged_d2 %>% 
  distinct(site, type_taxon, taxon, type, max_conc, pfas_type) %>% 
  filter(type != "Sediment") %>% 
  filter(type != "Detritus") %>%
  filter(type != "Water") %>% 
  filter(type != "Seston") %>% 
  mutate(taxon = case_when(is.na(taxon) ~ type,TRUE ~ taxon)) %>% 
  add_epred_draws(hg4_taxon, re_formula = NULL) %>% 
  mutate(.epred = .epred*max_conc) 

iso_post_summaries = iso_posts %>% 
  group_by(taxon, site, center_15n) %>% 
  reframe(mean_n15 = mean(.epred),
          median_n15 = median(.epred),
          sd_n15 = sd(.epred),
          lower_n15 = quantile(.epred, probs = 0.025),
          upper_n15 = quantile(.epred, probs = 0.975)) %>% 
  filter(site != "Russell Brook") # No spiders and only 2 dots for insects

pfas_conc_isotopes = pfas_posts %>% 
  filter(site != "Russell Brook") %>% # only 2 data points. Can fit regression here and model below was having trouble so I removed.
  group_by(site, taxon, type, pfas_type)  %>% 
  reframe(mean_conc = mean(.epred),
          median_conc = median(.epred),
          sd_conc = sd(.epred),
          lower_conc = quantile(.epred, probs = 0.025),
          upper_conc = quantile(.epred, probs = 0.975)) %>% 
  left_join(iso_post_summaries, relationship = "many-to-many") %>% 
  mutate(log10_mean_conc = log10(mean_conc),
         log10_sd_conc = sd(log10_mean_conc),
         log10_mean_conc_s = scale(log10_mean_conc),
         log10_median_conc = log10(median_conc),
         log10_median_conc_s = scale(log10_median_conc)) 

pfas_conc_isotopes %>%
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

brm_tmf_iso = brm(log10_median_conc_s ~ mean_n15 + (1 + mean_n15|site) + (1 + mean_n15|pfas_type),
                  data = pfas_conc_isotopes,
                  family = gaussian(),
                  prior = c(prior(normal(0, 1), class = Intercept),
                             prior(normal(0, 1), class = b),
                             prior(exponential(2), class = sd)),
                  iter = 2000, chains = 4)
# 
# saveRDS(brm_tmf_iso, file = "models/brm_tmf_iso.rds")


# # check model
pp_check(brm_tmf_iso, type = "dens_overlay_grouped", group = "pfas_type")

brm_tmf_iso = readRDS(file = "models/brm_tmf_iso.rds")

posts_tmf_iso = brm_tmf_iso$data %>% 
  distinct(site, pfas_type) %>% 
  expand_grid(mean_n15 = seq(-2, 1, length.out = 30)) %>% 
  add_epred_draws(brm_tmf_iso, re_formula = NULL)

mean_conc = attributes(pfas_conc_isotopes$log10_median_conc_s)$`scaled:center`
sd_conc = attributes(pfas_conc_isotopes$log10_median_conc_s)$`scaled:scale`

posts_tmf_iso_summary = posts_tmf_iso %>% 
  group_by(mean_n15, site, pfas_type) %>% 
  mutate(.epred_log10 = (.epred*sd_conc) + mean_conc) %>% 
  median_qi(.epred_log10) %>% 
  left_join(pfas_conc_isotopes %>% ungroup %>% distinct(site, pfas_type)) 
  
brm_tmf_iso_raw = brm_tmf_iso$data %>%
  left_join(pfas_conc_isotopes %>% ungroup %>% distinct(site, pfas_type)) %>% 
  mutate(.epred_log10 = (log10_median_conc_s*sd_conc) + mean_conc)

plot_tmf_isotopes = posts_tmf_iso_summary %>% 
  # filter(pfas_type == "6:2FTS") %>% 
  ggplot(aes(x = mean_n15, y = .epred_log10, fill = pfas_type)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line() +
  scale_fill_custom() + 
  scale_color_custom() +
  facet_grid2(site ~ pfas_type, scales = "free") +
  geom_point(data = brm_tmf_iso_raw,
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

ggsave(plot_tmf_isotopes, file = "plots/plot_tmf_isotopes.jpg",
       width = 10, height = 4)


# slopes ------------------------------------------------------------------

plot_tmf_isotopes2 = posts_tmf_iso_summary %>% 
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
  geom_point(data = brm_tmf_iso_raw,
             aes(color = pfas_type),
             shape = 1,
             size = 0.5) +
  scale_x_continuous(breaks = c(-1.5, 0, 1.5)) +
  theme(strip.text = element_text(size = 6),
        axis.text.x = element_text(size = 7)) +
  labs(x = expression(paste(delta^{15}, "N (centered)")),
       y = "PFAS Concentration (tissue ppb)") +
  NULL

ggsave(plot_tmf_isotopes2, file = "plots/plot_tmf_isotopes2.jpg",
       width = 10, height = 4)

tmf_iso_slopes = brm_tmf_iso$data %>% 
  distinct(site, pfas_type) %>% 
  # expand_grid(mean_n15 = seq(-2, 1.4, length.out = 2)) %>% # 3.4 per mill difference = 1 trophic level. Use -2 and 1.4 to keep within the values of centered n15
  expand_grid(mean_n15 = seq(0, 1, length.out = 2)) %>% 
  add_epred_draws(brm_tmf_iso, re_formula = NULL) %>%
  ungroup %>% 
  select(-.row, -.chain, -.iteration) %>% 
  pivot_wider(names_from = mean_n15, values_from = .epred) %>% 
  mutate(slope = `1` - `0`,
         slope = 10^slope)

tmf_slope_table = tmf_iso_slopes %>% 
  group_by(site, pfas_type) %>% 
  median_qi(slope)

write_csv(tmf_slope_table, file = "tables/tmf_slope_table.csv")

library(ggridges)
tmf_slopes = tmf_iso_slopes %>% 
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


ggsave(tmf_slopes, file = "plots/tmf_slopes.jpg", width = 6.5, height = 9)




# compare slopes with and without metamorphosis correction ----------------

tmf_slope_table_notadjustedformetamorphosis = read_csv(file = "tables/tmf_slope_table_notadjustedformetamorphosis.csv") %>% 
  select(site, pfas_type, slope, .lower, .upper) %>% 
  rename(slope_notadjusted = slope, 
         .lower_notadjusted = .lower,
         .upper_notadjusted = .upper)

tmf_compare = tmf_slope_table %>% 
  select(site, pfas_type, slope, .lower, .upper) %>% 
  rename(slope_adjusted = slope, 
         .lower_adjusted = .lower,
         .upper_adjusted = .upper) %>% 
  left_join(tmf_slope_table_notadjustedformetamorphosis) %>% 
  ggplot(aes(x = slope_adjusted,
             y = slope_notadjusted)) +
  geom_point() +
  geom_errorbar(aes(ymin = .lower_notadjusted, ymax = .upper_notadjusted),
                linewidth = 0.1, alpha = 0.5) +
  geom_errorbarh(aes(xmin = .lower_adjusted, xmax = .upper_adjusted),
                 linewidth = 0.1, alpha = 0.5) +
  geom_abline() +
  labs(x = "TMF slope (adjusted for metamorphosis fractionation)",
       y = "TMF slope (not adjusted for metamorphosis fractionation)")


ggsave(tmf_compare, file = "plots/tmf_compare.jpg", width = 5, height = 5)

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
