library(readxl)
library(tidyverse)
library(janitor)
library(brms)
library(tidybayes)
library(ggforce)
library(ggthemes)

# load data
isotopes = readRDS(file = "data/isotopes.rds")

brm_isotopes = readRDS(file = "models/brm_isotopes.rds")


# data plots ----------------------------------------------------------

isotopes %>% 
  pivot_longer(cols = c(d15n, original_d15n)) %>%
  ggplot(aes(x = taxon, color = name, y = value)) + 
  geom_point(position = position_dodge(width = 0.2)) +
  facet_wrap(~site)

# isotope_quick_plot = isotopes %>% 
#   ggplot(aes(x = d13c, y = d15n)) + 
#   geom_point(aes(color = family)) +
#   geom_mark_ellipse() +
#   geom_point(data = isotopes %>% filter(family == "Spider"), color = "black", shape = 13) +
#   geom_mark_ellipse(data = isotopes %>% filter(family == "Spider"), fill = "grey", shape = 13) +
#   facet_wrap(~site) +
#   ylim(-4, 20) +
#   xlim(-40, -20) +
#   scale_color_colorblind()
# 
# ggsave(isotope_quick_plot, file = "plots/isotope_quick_plot.jpg", width = 6.5, height = 4)

raw_n15 = isotopes %>% 
  group_by(type_taxon) %>% 
  mutate(median = median(mean_centered_15n)) %>% 
  ggplot(aes(x = reorder(taxon, median),
             y = mean_centered_15n,
             color = taxon)) + 
  geom_point() +
  labs(y = "\u03b4 15N",
       x = "") +
  guides(color = "none") +
  facet_wrap(~site) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 0))

# ggsave(raw_n15, file = "plots/raw_n15.jpg", width = 6.5, height = 4)

# fit model ---------------------------------------------------------------

brm_isotopes_notadjustedformetamorphosis = brm(mean_centered_15n ~ taxon + (1 + taxon|site),
                   family = gaussian(),
                   prior = c(prior(normal(0, 1), class = Intercept),
                             prior(normal(0, 1), class = b),
                             prior(exponential(2), class = sd)),
                   data = isotopes)

saveRDS(brm_isotopes_notadjustedformetamorphosis, file = "models/brm_isotopes_notadjustedformetamorphosis.rds")
plot(conditional_effects(brm_isotopes_notadjustedformetamorphosis), points = T)

# posteriors --------------------------------------------------------------
brm_isotopes_notadjustedformetamorphosis = readRDS(file = "models/brm_isotopes_notadjustedformetamorphosis.rds")

baseline_n15 = readRDS(file = "posteriors/iso_posts_notadjustedformetamorphosis.rds") %>% 
  filter(taxon == "Biofilm") %>% 
  mutate(baseline_n15_raw = .epred + center_15n) %>% 
  ungroup %>% 
  select(site, .draw, baseline_n15_raw) 

iso_posts_notadjustedformetamorphosis = isotopes %>%
  distinct(taxon, site, mean_15n) %>% 
  rename(center_15n = mean_15n) %>% 
  add_epred_draws(brm_isotopes_notadjustedformetamorphosis, re_formula = ~ NULL) %>% 
  group_by(taxon) %>% 
  mutate(median = median(.epred)) %>% 
  left_join(baseline_n15) %>% 
  mutate(raw_n15 = .epred + center_15n,
         trophic_level = 1 + (raw_n15 - baseline_n15_raw)/3.4)

saveRDS(iso_posts_notadjustedformetamorphosis, file = "posteriors/iso_posts_notadjustedformetamorphosis.rds")

iso_posts_notadjustedformetamorphosis %>% 
  group_by(site) %>% 
  ggplot(aes(x = reorder(taxon, median),
             y = .epred)) + 
  stat_halfeye() +
  labs(y = "\u03b4 15N",
       x = "") +
  facet_wrap(~site) +
  geom_point(data = isotopes %>% left_join(iso_posts_notadjustedformetamorphosis %>% ungroup %>% distinct(taxon, median)), 
             aes(y = mean_centered_15n),
             fill = "white", shape = 21) +
  ylim(-4,4)


iso_posts_notadjustedformetamorphosis %>% 
  group_by(site, taxon) %>% 
  mutate(raw_n15 = .epred + center_15n) %>% 
  median_qi(raw_n15)
