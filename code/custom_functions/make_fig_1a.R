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

make_fig_1a = function(mod = NULL, mod_notaxa){
  mod1 = mod_notaxa
  mod1_taxa = mod
  
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
  
posts_sumpfas_summary %>%
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
    geom_line(data = line_posts %>% filter(type %in% c("Water", "Biofilm", "Larval", "Adult", "Tetragnathidae")),
              aes(group = group)) + 
    scale_fill_brewer(type = "qual", palette = 7) +
    scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
    NULL
}