### Rrs spectra analysis from Sentinel-2 images of Mesodinium blooms ###

# V. POCHIC
# 2025/11/18

## Packages
library(tidyverse)
library(ggnewscale)

# Import data ####

# From the Loire bloom
spectra_loire <- read.csv2('Data/Satellite/Rrs_spectra/Spectral_data_20170412_Loire_Meso.csv', header = TRUE,
                            fileEncoding = 'ISO-8859-1') %>%
  # group
  group_by(Category, Replicate) %>%
  # Create an 'Individual' variable for plotting
  mutate(Individual = paste(Site, Category, Replicate, collapse = '_')) %>%
  # Factorise and reorder Category
  mutate(Category = as_factor(Category)) %>%
  mutate(Category = fct_relevel(Category, c('Ocean', 'Mesodinium')))

# From the Seine bloom
spectra_seine <- read.csv2('Data/Satellite/Rrs_spectra/Spectral_data_20170820_Seine_Meso.csv', header = TRUE,
                           fileEncoding = 'ISO-8859-1') %>%
  # group
  group_by(Category, Replicate) %>%
  # Create an 'Individual' variable for plotting
  mutate(Individual = paste(Site, Category, Replicate, collapse = '_')) %>%
  # Factorise and reorder Category
  mutate(Category = as_factor(Category)) %>%
  mutate(Category = fct_relevel(Category, c('Ocean', 'Mesodinium')))

# Summarise data ####

## Loire bloom
spectra_loire_sum <- spectra_loire %>%
  group_by(Category, Wavelength, Site) %>%
  select(-c('Band','Replicate')) %>%
  summarise(Rrs_mean = mean(Rrs), Rrs_sd = sd(Rrs), .groups = 'keep') %>%
  # Factorise and reorder Category
  mutate(Category = as_factor(Category)) %>%
  mutate(Category = fct_relevel(Category, c('Ocean', 'Mesodinium')))

## Seine bloom
spectra_seine_sum <- spectra_seine %>%
  group_by(Category, Wavelength, Site) %>%
  select(-c('Band','Replicate')) %>%
  summarise(Rrs_mean = mean(Rrs), Rrs_sd = sd(Rrs), .groups = 'keep') %>%
  # Factorise and reorder Category
  mutate(Category = as_factor(Category)) %>%
  mutate(Category = fct_relevel(Category, c('Ocean', 'Mesodinium')))

# Join both datasets
spectra_join <- bind_rows(spectra_loire, spectra_seine)
spectra_join_sum <- bind_rows(spectra_loire_sum, spectra_seine_sum)

# Plot data ####

# Color palette
palette_2 <- c('#000E53','#711412')

# Plot
ggplot() +
  ## Rrs spectra
  # Ribbon for sd
  geom_ribbon(data = spectra_join_sum, aes(x = Wavelength,
                                              ymin = Rrs_mean - Rrs_sd,
                                              ymax = Rrs_mean + Rrs_sd,
                                              fill = Category), alpha = .15) +
  # Lines for individual replicates
  geom_path(data = spectra_join, aes(x = Wavelength, y = Rrs, 
                                        color = Category, group = Individual),
            alpha = .35, linewidth = .4) +
  # Lines for means
  geom_path(data = spectra_join_sum, aes(x = Wavelength, y = Rrs_mean, 
                                            color = Category, group = Category),
            alpha = .7, linewidth = .8) +
  # Points for means
  geom_point(data = spectra_join_sum, aes(x = Wavelength, y = Rrs_mean, 
                                             color = Category),
             shape = 19, size = 3) +
  # color scale
  scale_color_discrete(type = palette_2, labels = c(
  'Clear
waters', 
  'Mesodinium
bloom')) +
  scale_fill_discrete(type = palette_2, guide = 'none') +
  
  # Facets by site
  facet_wrap(facets = 'Site', nrow = 2) +
  
  # X-axis scaling
  scale_x_continuous(limits = c(400,900)) +
  # Labels
  labs(title = NULL, y = 'Rrs (sr-1)',
       x = 'Wavelength (nm)', color = 'Category:') +
  # Theme
  theme_classic() +
  theme(legend.position = 'bottom')

# Nice!
# Save that
# ggsave('Plots/Satellite/Rrs_spectra_meso_blooms.tiff', height = 95, width = 82,
#        units = 'mm', dpi = 300, compression = 'lzw')
