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

make_fig_1 <- function(mod = NULL){
  
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
                                 center_interval = "median_cri", digits = 2) {
    df %>%
      mutate(
        .center_val = signif(.data[[center]], digits),
        .lower_val = signif(.data[[lower]], digits),
        .upper_val = signif(.data[[upper]], digits),
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
  merged_d2 = mod$data2
  
  max_conc = unique(merged_d2$max_conc)
  
  pfas_names = read_csv("data/pfas_names.csv")
  
  mod_dat = mod$data %>% 
    mutate(pfas_type = str_replace(pfas_type, "\\.", "\\:")) %>% 
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
    add_epred_draws(mod, re_formula = NULL, dpar = T) %>%
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
    add_epred_draws(seed = 20202, mod, re_formula = NULL, dpar = T) %>%
    mutate(.epred = .epred*unique(mod$data2$max_conc_ppb)) 
  
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
  
  mod_dat %>%
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
  
}
