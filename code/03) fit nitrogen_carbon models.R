library(readxl)
library(tidyverse)
library(janitor)
library(brms)
library(tidybayes)

nitrogen_carbon = read_excel("data/MergedDataFile_Farmington Stable Isotopes 2022-2023.xlsx", 
                      sheet = "Merged Data Sheet") %>% 
  clean_names() %>% 
  filter(media == "Invertebrate") %>%
  pivot_longer(cols = c(starts_with("per"))) %>% 
  group_by(name) %>% 
  mutate(nutrient_s = scale(value),
         mean_nut = mean(value, na.rm = T),
         sd_nut = sd(value, na.rm = T))

brm_nitrogen = readRDS(file = "models/brm_nitrogen.rds")


# fit model ---------------------------------------------------------------
# brm_nitrogen = update(brm_isotopes, formula = nutrient_s ~ family + (1 | site) + (1 | lifestage) + (1|name), newdata = nitrogen_carbon)
# 
# saveRDS(brm_nitrogen, file = "models/brm_nitrogen.rds")
# 

# posteriors --------------------------------------------------------------

posts_nc = nitrogen_carbon %>% 
  distinct(mean_nut, sd_nut, name, family, site) %>% 
  add_epred_draws(brm_nitrogen, re_formula = ~ (1|name) + (1|site)) %>% 
  mutate(value = (.epred*sd_nut) + mean_nut)

posts_nc %>% 
  ggplot(aes(x = reorder(family, mean_nut), 
             y = value)) +
  facet_wrap(site~name, scales = "free_y") +
  stat_halfeye() +
  geom_point(data = nitrogen_carbon)
