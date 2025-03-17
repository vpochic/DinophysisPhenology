### SCOPUS Dinophysis ###
# V. POCHIC, 2025-03-17

# A little script to plot data from SCOPUS on Dinophysis and other 
# dinoflagellate genera in the scientific literature.

# When you compare the trend of documents in SCOPUS mentioning Dinophysis, 
# it looks like it stands out from other, non-toxic dinoflagellates 
# (i.e. Ornithocercus, Protoperidinium in our case study).

# All 3 genera were mentioned in SCOPUS documents (albeit barely) BEFORE 1980, 
# the year Dinophysis fortii was revealed to be toxic by Yasumoto and his 
# colleagues.

# packages ####
library(tidyverse)

## import and merge data ####

# SCOPUS results ####

# for search on:

# "Dinophysis"
SCOPUS_dinophysis <- read.csv(file = 'SCOPUS/data/Scopus-10-Analyze-Year_all-fields_Dinophysis.csv',
                                  header = FALSE) %>%
  # Remove first 4 rows that are useless
  slice_tail(n=-4) %>%
  # Create a variable for the search term
  mutate(term = 'Dinophysis')

# add column names
colnames(SCOPUS_dinophysis) <- c('Year', 'n', 'term')

# "Ornithocercus"
SCOPUS_ornitho <- read.csv(file = 'SCOPUS/data/Scopus-10-Analyze-Year_all-fields_Ornithocercus.csv',
                       header = FALSE) %>%
  # Remove first 4 rows that are useless
  slice_tail(n=-4) %>%
  # Create a variable for the search term
  mutate(term = 'Ornithocercus')

# add column names
colnames(SCOPUS_ornitho) <- c('Year', 'n', 'term')

# "Protoperidinium"
SCOPUS_protop <- read.csv(file = 'SCOPUS/data/Scopus-10-Analyze-Year_all-fields_Protoperidinium.csv',
                           header = FALSE) %>%
  # Remove first 4 rows that are useless
  slice_tail(n=-4) %>%
  # Create a variable for the search term
  mutate(term = 'Protoperidinium')

# add column names
colnames(SCOPUS_protop) <- c('Year', 'n', 'term')

# vector of years ####

# Create a vector of years to introduce zeros in the SCOPUS datasets (for years
# with 0 articles), for each of the 4 terms ('dinoflagellate' + 3 genus names)

# The 'Year' vector starts in 1937 because that is the first year a document in
# SCOPUS mentions either of the 3 genera
Years_dinophysis <- as_tibble(seq(from = 1937, to = 2025, by = 1)) %>%
  mutate(term = 'Dinophysis')
Years_protop <- as_tibble(seq(from = 1937, to = 2025, by = 1)) %>%
  mutate(term = 'Protoperidinium')
Years_ornitho <- as_tibble(seq(from = 1937, to = 2025, by = 1)) %>%
  mutate(term = 'Ornithocercus')

# Merge all 4
Years = bind_rows(Years_dinophysis, Years_ornitho, Years_protop)

# Appropriate column names
colnames(Years) <- c('Year', 'term')

# Mutate 'Year' as character to fit SCOPUS datasets
Years <- Years %>%
  mutate(Year = as.character(Year))

# merge ####

# Merge the SCOPUS datasets and the Years tibble
SCOPUS_dino <- right_join(
  # All the SCOPUS datasets in 1 tibble
  bind_rows(SCOPUS_dinophysis, SCOPUS_ornitho, SCOPUS_protop), 
  # The Years tibble
                         Years,
  # By 'Year' and 'term'
  by = c('Year', 'term')) %>%
  
  # mutate 'Year' as a number
  mutate(Year = as.integer(Year)) %>%
  
  # Replace NAs with 0s
  mutate(n = ifelse(is.na(n) == TRUE, 0, n)) %>%
  
  # Term as factor so we order it as we want
  mutate(term = as_factor(term)) %>%
  mutate(term = fct_relevel(term, c('Dinophysis', 'Protoperidinium', 
                                    'Ornithocercus'))) %>%
  arrange(desc(Year))

### plot ####

# 3 step color scale
palette_3_steps <- c('#711412', '#FC4D6B', '#2E6CD9')

# Plot all years prior to 2025
ggplot(data = subset(SCOPUS_dino, Year < 2025)) +
  # Plotting data
  # geom_point(aes(x = Year, y = n, color = term), size = 4) +
  geom_path(aes(x = Year, y = n, 
                color = term), linewidth = 1, alpha = .9) +
  # x-axis breaks and limits
  scale_x_continuous(limits = c(1930,2024),
                     breaks = c(1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980,
                                1990, 2000, 2010, 2020)) +
  # color scale and legend title
  scale_color_discrete(type = palette_3_steps,
                       guide = guide_legend(title = 'Search term:')) +
  
  # Plot and axis titles
  labs(title = 'Documents in SCOPUS mentioning 
different dinoflagellate genera', 
       y = 'Number of articles', x = 'Year'
  ) +
  #Theme
  theme_classic() +
  theme(legend.position = 'bottom')

# save the plot
ggsave('SCOPUS_articles.tiff', height = 130, width = 164, units = 'mm',
       dpi = 300, compression = 'lzw')
