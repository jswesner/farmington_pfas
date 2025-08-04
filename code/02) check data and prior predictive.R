library(brms)
library(tidyverse)
library(tidybayes)
library(janitor)

# proportion of zeros in modeled data
merged_d2 %>% filter(conc_ppb == 0) %>% nrow/nrow(merged_d2)

# proportion of zeros in full data (i.e., all PFAS compounds)
full_data = read_csv("data/Farmington_PFAS_FWstudy_Datarelease.csv") %>% 
  clean_names() %>% 
  rename(sample_id = field_id,
         pfas_type = compound) %>% 
  mutate(nondetects = case_when(conc == "ND" ~ "ND", # not detected
                                conc == "NM" ~ "NM", # not measured
                                conc == 'na' ~ 'na')) %>%  # not available
  mutate(conc_ppb = parse_number(conc)) %>% 
  mutate(conc_ppb = case_when(nondetects == "ND" ~ 0, TRUE ~ conc_ppb)) 

full_data %>% 
  filter(conc_ppb == 0) %>% nrow/nrow(full_data)

length(unique(full_data$pfas_type)) # number of pfas compounds tested for = 28

full_data %>% 
  group_by(pfas_type) %>% 
  add_tally(name = "total") %>% 
  filter(conc_ppb > 0) %>% 
  group_by(pfas_type, total) %>% 
  tally() %>% 
  arrange(-n) %>% 
  mutate(prop = n/total)
  

# proportion of zeros by group --------------------------------------------

full_data %>% 
  group_by(sample_type) %>% 
  mutate(is_zero = case_when(conc_ppb == 0 ~ 1,TRUE ~ 0)) %>% 
  add_tally() %>% 
  group_by(sample_type, n) %>% 
  reframe(zero = sum(is_zero)) %>% 
  mutate(prop = zero/n) %>% 
  arrange(-prop)

# prior simulation --------------------------------------------------------

#            prior=c(prior(normal(-2,1), class = Intercept), # -5,2
#                    prior(normal(0,0.5), class = b), # was (0,2), then (0,0.1), then (0,0.5)
#                    prior(normal(0,1), class = b, dpar = hu),
#                    prior(normal(-1.5, 1), class = Intercept, dpar= hu)), # was (-2,1), then (-1.5,1)
n = 1000
tibble(alpha_lin = rnorm(n, -2, 1),
       beta_lin = rnorm(n, 0, 0.5),
       alpha_logit = rnorm(n, -1.5, 1),
       beta_logit = rnorm(n, 0, 1)) %>% 
  mutate(intercept_linear = exp(alpha_lin + beta_lin*0),
         intercept_logit = inv_logit_scaled(alpha_logit + beta_logit*0),
         group1_linear1 = exp(alpha_lin + beta_lin*1),
         group1_logit1 = inv_logit_scaled(alpha_logit + beta_logit*1)) %>% 
  mutate(draw = 1:nrow(.)) %>% 
  select(starts_with("intercept"), starts_with("group"), draw) %>% 
  pivot_longer(cols = -draw) %>% 
  separate(name, into = c("response", "model")) %>% 
  ggplot(aes(x = response, y = value)) + 
  stat_halfeye() +
  facet_wrap(~model, scales = "free") +
  scale_y_log10() +
  NULL
  
