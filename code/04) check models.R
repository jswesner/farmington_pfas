library(tidyverse)
library(brms)
library(tidybayes)

hg4_taxon = readRDS(file = "models/hg4_taxon.rds")

merged_d2 = readRDS("data/merged_d2.rds") %>% 
  mutate(type_taxon = sample_type) %>% 
  separate(sample_type, into = c('type', 'taxon')) %>% 
  group_by(type_taxon) %>% 
  mutate(max_conc = max(conc_ppb, na.rm = T)) %>% 
  mutate(conc_ppb_s = conc_ppb/max(conc_ppb))

pp_check(hg4_taxon) + scale_x_log10()

test = plot(conditional_effects(hg4_taxon), points = T)

test$type_taxon + scale_y_log10()

post = hg4_taxon$data %>% 
  select(-conc_ppb_s) %>% 
  distinct() %>% 
  left_join(merged_d2 %>% ungroup %>% distinct(type_taxon, max_conc)) %>%
  # mutate(max_conc = hg4_taxon$data2$max_conc_ppb) %>% 
  add_epred_draws(hg4_taxon) %>% 
  group_by(type_taxon, site) %>%
  mutate(.epred = .epred*max_conc) %>% 
  median_qi(.epred)

post %>% 
  ggplot(aes(x = type_taxon, y = .epred)) +
  geom_pointrange(aes(y = .epred, ymin = .lower, ymax = .upper)) +
  facet_wrap(~site) +
  scale_y_log10() +
  geom_point(data = merged_d2, aes(y = conc_ppb)) +
  theme(axis.text.x = element_text(angle = 90))
