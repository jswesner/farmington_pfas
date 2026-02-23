library(tidyverse)
library(brms)
library(readxl)

#~~~~This code refits the ppb model with just diptera and trichoptera. 
#~~~~It then compares it to the original model with all taxa
#~~~~The purpose is to compare the influence of larger taxa on the whole model in response to a reviewer's comments
# about the possible influence of large taxa that are not eaten by spiders

hg4_taxon = readRDS(file = "models/hg4_taxon.rds")    # pfas ppb model with all taxa
hg4_taxon_diptric = readRDS(file = "models/temporary/hg4_taxon_diptric.rds") # pfas ppb model with select taxa
# pfas ppb refit ----------------------------------------------------------------
# load data
merged_d2_diptric = readRDS("data/merged_d2.rds") %>% filter(is.na(taxon) |  taxon %in% c("Diptera", "Trichoptera"))

# fit model
# 
# hg4_taxon <- brm(bf(conc_ppb_s ~ type_taxon + (1 + type_taxon|site) + (1 + type_taxon|pfas_type),
#               hu ~ type_taxon + (1 + type_taxon|site) + (1 + type_taxon|pfas_type)),
#            data = merged_d2,
#            family = hurdle_gamma(link = "log"),
#            prior=c(prior(normal(-2,1), class = Intercept), # -5,2
#                    prior(normal(0,0.5), class = b), # was (0,2), then (0,0.1), then (0,0.5)
#                    prior(normal(0,1), class = b, dpar = hu),
#                    prior(normal(-1.5, 1), class = Intercept, dpar= hu)), # was (-2,1), then (-1.5,1)
#            iter = 2000 , chains = 4, cores = 4,
#            seed = 5,
#            # save_pars = save_pars(all = TRUE),
#            data2 = merged_d2))
# 
# saveRDS(hg4_taxon, file = "models/hg4_taxon.rds")

# hg4_taxon = readRDS(file = "models/hg4_taxon.rds")
# hg4_taxon_diptric = update(hg4_taxon, newdata = merged_d2_diptric, chains = 1, iter = 500, data2 = merged_d2_diptric)
# 
# saveRDS(hg4_taxon_diptric, file = "models/temporary/hg4_taxon_diptric.rds")




# sum pfas refit ----------------------------------------------------------
# Refit sum ppb for just diptera and trichoptera
# mod1_taxa = readRDS(file = "models/mod1_taxa.rds")    # sum pfas per taxon
# merged_d2_sum_taxa_diptric = readRDS(file = "data/merged_d2_sum_taxa_diptric.rds") 
# 
# mod1_taxa_diptric = update(mod1_taxa, newdata = merged_d2_sum_taxa_diptric, data2 = merged_d2_sum_taxa_diptric)
# saveRDS(mod1_taxa_diptric, file = "models/temporary/mod1_taxa_diptric.rds")
mod1_taxa_diptric = readRDS(file = "models/temporary/mod1_taxa_diptric.rds")

# Refit sum ppb for all trophic groups, but including only just diptera and trichoptera for emergers and larvae
merged_d2_sum_diptric = readRDS(file = "data/merged_d2_sum_diptric.rds")
# mod1 = readRDS(file = "models/mod1.rds")  
# mod1_diptric = update(mod1, newdata = merged_d2_sum_diptric, data2 = merged_d2_sum_diptric)
# saveRDS(mod1_diptric, file = "models/temporary/mod1_diptric.rds")
mod1_diptric = readRDS(file = "models/temporary/mod1_diptric.rds")


# plot --------------------------------------------------------------------
source("code/custom_functions/make_fig_1b.R") # load code to remake fig 1b
source("code/custom_functions/make_fig_1a.R") # load code to remake fig 1a

# fig1b_original = make_fig_1b(mod = hg4_taxon)
fig1b_diptric = make_fig_1(mod = hg4_taxon_diptric)

fig1a_original = make_fig_1a(mod_notaxa = mod1, mod = mod1_taxa)

diptric_lines = fig1a_diptric$layers[[3]]$data
diptric_points = fig1a_diptric$data
fig1a_diptric = fig1a_original +
  geom_pointrange(data = diptric_points,
                  aes(ymin = .lower, ymax = .upper),
                  shape = 24,
                  size = 0.6) +
  geom_line(data = diptric_lines, aes(x = type, y = sum_ppb,
                                      group = group))

ggsave(fig1a_diptric, file = "plots/temporary/fig1a_diptric.jpg", width = 6, height = 6)

# sum tmf -----------------------------------------------------------------
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

posts_sumpfas_diptric = mod1_diptric$data2 %>% 
  distinct(type, site, mean_sum_ppb) %>%
  mutate(order = case_when(type == "Water" ~ 1,
                           type == "Sediment" ~ 2,
                           type == "Biofilm" ~ 5,
                           type == "Detritus" ~ 3,
                           type == "Seston" ~ 4,
                           type == "Larval" ~ 6,
                           type == "Emergent" ~ 7,
                           type == "Tetragnathidae" ~ 8)) %>% 
  add_epred_draws(seed = 20202, mod1_diptric, re_formula = NULL) %>% 
  mutate(sum_ppb = (.epred - 0.0001)*mean_sum_ppb) %>% 
  mutate(sum_ppb = case_when(sum_ppb < 0.000101 ~ 0, TRUE ~ sum_ppb - 0.0001)) %>% 
  group_by(type, .draw, order) %>% 
  reframe(sum_ppb = mean(sum_ppb)) %>% 
  mutate(type = case_when(type == "Emergent" ~ "Adult", T ~ type))


post_sumpfas_ratios_diptric = posts_sumpfas_diptric %>% 
  select(-order) %>% 
  pivot_wider(names_from = type, values_from = sum_ppb) %>% 
  mutate(a_bw = Biofilm/Water,
         b_lb = Larval/Biofilm,
         c_al = Adult/Larval,
         d_ta = Tetragnathidae/Adult) %>% 
  select(.draw, a_bw, b_lb, c_al, d_ta) %>% 
  pivot_longer(cols = -.draw)


post_sumpfas_ratios_diptric %>% 
  group_by(name) %>% 
  median_qi(value)

post_sumpfas_ratios_diptric %>% 
  group_by(name) %>% 
  reframe(prob = sum(value>1)/max(.draw))
