library(readxl)
library(tidyverse)
library(janitor)
library(brms)
library(tidybayes)
library(ggforce)
library(ggthemes)

isotopes = read_excel("data/MergedDataFile_Farmington Stable Isotopes 2022-2023.xlsx", 
                       sheet = "Merged Data Sheet") %>% 
  clean_names() %>% 
  filter(media == "Invertebrate") %>% 
  group_by(site) %>%
  mutate(d15n_s = scale(d15n),
         mean_15n = mean(d15n, na.rm = T),
         sd_15n = sd(d15n, na.rm = T)) %>% 
  ungroup

brm_isotopes = readRDS(file = "models/brm_isotopes.rds")


# bivariate plot ----------------------------------------------------------

isotope_quick_plot = isotopes %>% 
  filter(family != "Spider") %>% 
  filter(lifestage != "Larval") %>% 
  ggplot(aes(x = d13c, y = d15n)) + 
  geom_point(aes(color = family)) +
  geom_mark_ellipse() +
  geom_point(data = isotopes %>% filter(family == "Spider"), color = "black", shape = 13) +
  geom_mark_ellipse(data = isotopes %>% filter(family == "Spider"), fill = "grey", shape = 13) +
  facet_wrap(~site) +
  ylim(-4, 20) +
  xlim(-40, -20) +
  scale_color_colorblind()


ggsave(isotope_quick_plot, file = "plots/isotope_quick_plot.jpg", width = 6.5, height = 4)

# fit model ---------------------------------------------------------------

# brm_isotopes = brm(d15n_s ~ family + (1|site) + (1|lifestage),
#                    family = gaussian(),
#                    prior = c(prior(normal(0, 1), class = Intercept),
#                              prior(normal(0, 1), class = b),
#                              prior(exponential(2), class = sd)),
#                    data = isotopes)

# saveRDS(brm_isotopes, file = "models/brm_isotopes.rds")
# brm_isotopes = readRDS(file = "models/brm_isotopes.rds")
# plot(conditional_effects(brm_isotopes), points = T)

# posteriors --------------------------------------------------------------

posts = isotopes %>%
  distinct(family, mean_15n, sd_15n, lifestage, site) %>% 
  add_epred_draws(brm_isotopes, re_formula = ~ (1|site) + (1|lifestage)) %>% 
  ungroup %>% 
  # mutate(.epred = (.epred*sd_15n) + mean_15n) %>% 
  group_by(family, site) %>% 
  mutate(median = median(.epred))


posts %>% 
  ggplot(aes(x = reorder(family, median),
             y = .epred,
             fill = lifestage)) + 
  stat_halfeye() +
  labs(y = "\u03b4 15N",
       x = "") +
  facet_wrap(~site) +
  geom_point(data = isotopes %>% left_join(posts %>% ungroup %>% distinct(family, median)), 
             aes(y = d15n_s, 
                 color = lifestage), 
             fill = "white", shape = 21)

iso_post_summaries = posts %>% 
  filter(lifestage == "Adult") %>% 
  group_by(family, site) %>% 
  median_qi(.epred) %>% 
  rename(delta_n15_s = .epred) %>% 
  mutate(taxon = family)


# regress with PFAS -------------------------------------------------------
# load PFAS sum posteriors 
posts_taxon = readRDS("posteriors/posts_taxon.rds")

posts_taxon %>% 
  filter(!is.na(taxon)) %>% 
  filter(.draw <= 100) %>% 
  group_by(type, site, taxon, .draw, order) %>% 
  reframe(.epred = sum(.epred)) %>% 
  left_join(iso_post_summaries) %>% 
  ggplot(aes(x = delta_n15_s, y = .epred)) + 
  # geom_point() +
  stat_halfeye(aes(group = delta_n15_s)) +
  scale_y_log10() +
  facet_wrap(~site)
