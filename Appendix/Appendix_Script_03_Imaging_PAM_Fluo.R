#### Additional script: Imaging PAM fluorometer data ###
# Author: V. POCHIC
# Last modif: 2026/07/06

## General information ####

## Description: 
# In this script we plot the Fv/Fm data acquired on Dinophysis cultures with
# an Imaging PAM fluorometer during our starvation experiment in 2024

## Files required:
# Data/Starvation_Experiment/2024/PAM_fluo/[...].csv
## --> PAM fluorometer data for 2 cultures (Fed and Unfed) on day 44 of the
# experiment

## /!\ Please note that these are data acquired with an instrument and a method
# we were not completely familiar with. Hence, these are preliminary data and
# the files processed here do not reflect the entire set of
# measurements, which is why we chose not to publish them in our article.
# Use and/or interpret these data at your own peril! /!\

## Outputs:
## Figures:
# Figures A05 - Plot of the PSII Yield as a function of irradiance for Fed and 
# Unfed Dinophysis cultures on day 44 of our starvation experiment.

#### Required packages ####
library(tidyverse)

####-----------------------------------------------------------------------####
#### Import and curate data ####

Fed_data <- read.csv2('Appendix/Data/PAM_Fluo/Fed1_20240828.csv',
                      header = TRUE, fileEncoding = 'ISO-8859-1') %>%
  mutate(Condition = 'Fed')
  

Unfed_data <- read.csv2('Appendix/Data/PAM_Fluo/Unfed2_20240828.csv',
                      header = TRUE, fileEncoding = 'ISO-8859-1') %>%
  mutate(Condition = 'Unfed')

# Merge datasets
PAM_data <- bind_rows(Fed_data, Unfed_data)

# Compute some stats
PAM_data_stats <- PAM_data %>%
  group_by(Condition, PAR) %>%
  summarise(Y_PSII.mean = mean(Y_PSII),
            Y_PSII.sd = sd(Y_PSII),
            .groups = 'keep')

#### Plot ####

# We will have a rather simple approach here: we'll plot the PSII Yield data
# for the 2 conditions, with mean +- sd

palette_2 <- c('grey10', 'red2')

ggplot() +
  # First, individual points
  geom_point(data = PAM_data, aes(x = PAR, y = Y_PSII, color = Condition),
           size = 2, alpha = .35, stroke = 1, shape = 21) +
  # Then, mean + sd
  geom_errorbar(data = PAM_data_stats, aes(x = PAR, 
                                           ymin = Y_PSII.mean - Y_PSII.sd,
                                           ymax = Y_PSII.mean + Y_PSII.sd),
                width = 10, color = 'black') +
  geom_point(data = PAM_data_stats, aes(x = PAR, y = Y_PSII.mean,
                                        color = Condition), 
             size = 5, stroke = 1.25, shape = 21) +
  
  # Color palette
  scale_color_discrete(palette = palette_2) +
  
  # Labels
  labs(x = 'Light (ppfd)', y = 'PSII Quantum Yield',
       color = 'Condition: ') +
  
  # Theme
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.background = element_rect(fill = '#F8F7F5', color = 'grey10',
                                         linewidth = .25),
        legend.frame = element_rect(fill = 'transparent', color = 'grey10',
                                    linewidth = .25))

# I've seen better plots, but it's not awful. That will do.
# ggsave('Appendix/Plots/Imaging_PAM_Fluo/FigA05_Imaging_PAM_Fluo_YPSII_day44.tiff',
#        height = 100, width = 150, units = 'mm', dpi = 600, compression = 'lzw')
  
