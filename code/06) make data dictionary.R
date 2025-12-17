
# 1) get readme files
readme_files <- list.files(path = "data", 
                           pattern = "README.csv$", full.names = TRUE)

read_and_clean <- function(file) {
  df <- read_csv(file)
  # Remove columns that are entirely NA or blank
  df[, colSums(!is.na(df) & df != "") > 0]
}

readme_data <- lapply(readme_files, read_and_clean)

# 2) get data files
data_files <- 
  list.files(
    path = "data",
    pattern = "(Datarelease.csv$|InsectMasses.csv$)", # gets files ending in Datarelease OR InsectMasses
    full.names = TRUE
  )

raw_data <- lapply(data_files, read_and_clean)

# 3) extract data types (integer, numeric, etc.)
data_types_list = NULL
  
for(i in 1:length(raw_data)){
  data_types_list[[i]] = tibble(
  col= names(raw_data[[i]]),
  data_type = sapply(raw_data[[i]], function(col) class(col)[1])  # Extract the primary class of each column
) %>% rename(`Column name` = col)
}

data_types = purrr::map2(readme_data, data_types_list, ~ left_join(.x, .y)) 


# 4) extract data ranges

data_ranges = NULL
for(i in 1:length(raw_data)){
  data_ranges[[i]] = raw_data[[i]] %>%
    mutate(across(everything(), as.character)) %>% 
    pivot_longer(cols = everything()) %>% 
    arrange(name) %>% 
    group_by(name) %>%
    rename(`Column name` = name) %>% 
    summarise(value = paste(unique(value), collapse = ", "))
}

# 5) make data dictionary
# first for the four data files
data_dictionary = NULL  

for(i in 1:length(readme_data)){
  data_dictionary[[i]] = readme_data[[i]] %>% 
  drop_na("Column name") %>% 
  left_join(data_ranges[[i]]%>% 
              drop_na("Column name"), by = "Column name") %>% 
  left_join(data_types[[i]] %>% select(-Description) %>% 
              drop_na("Column name"), by = "Column name") %>% 
  mutate(dataset = basename(data_files[i])) %>% 
  select(dataset, `Column name`, data_type, everything())
}

# then add the site location dictionary
site_locations = read_csv("data/SiteLocations.csv") %>%
  mutate(across(everything(), as.character)) %>%  
  pivot_longer(cols = everything()) %>% 
  group_by(name) %>% 
  rename(`Column name` = name) %>% 
  summarize(value = paste(unique(value), collapse = ", ")) %>% 
  mutate(dataset = "SiteLocations.csv",
         data_type = c("numeric", "numeric", "character", "character", "character"))

data_dictionary[[5]] = site_locations

write_csv(bind_rows(data_dictionary), file = "data/data_dictionary.csv")

                  