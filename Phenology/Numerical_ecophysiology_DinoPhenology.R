#### Numerical ecophysiology :) ##
## V. POCHIC
# 2024-04-24

# /!\ This script requires data tables generated with the 'Deriv_GAMs_DinoPheno-
# logy' and 'Dino phenology analysis' scripts /!\

### What is the relationship between Dinophysis accumulation/loss rate and
# different parameters ?

#### Packages ####
library(tidyverse)

#### Import data ####

# Table with temperature
Table_temp <- read.csv2('Table_hydro_daily_20240424.csv', header = TRUE,
                        fileEncoding = 'ISO-8859-1')
# Table with Dinophysis GAM derivatives
Deriv_dino <- read.csv2('Deriv_dino_20240424.csv', header = TRUE,
                        fileEncoding = 'ISO-8859-1')

# Modelled Mesodinium counts
response_pred_plot_Meso <- read.csv2('response_pred_plot_MESO_20240424.csv', 
                                     header = TRUE,
                                     fileEncoding = 'ISO-8859-1')

# Modelled Dinophysis counts
response_pred_plot_Dino <- read.csv2('response_pred_plot_DINO_20240424.csv', 
                                     header = TRUE,
                                     fileEncoding = 'ISO-8859-1')

#### How does Dinophysis accumulation rates vary with temperature ? ####

# Join the two tables
Deriv_dino_temp <- left_join(Table_temp, Deriv_dino,
                             by = c('Code_point_Libelle', 'Day'))
# Just need to plot this

# Some aesthetics

# Relevel factor to order sampling sites
Deriv_dino_temp <- Deriv_dino_temp %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

# a nice color palette
pheno_palette12 <- c('red3', 'orangered1', 'dodgerblue4', 'dodgerblue1', 
                     'chartreuse4', 'chartreuse2','goldenrod3', 'darkgoldenrod1',
                     'darkorchid4', 'darkorchid1', 'firebrick1', 'deeppink2'
)

# And let's plot
# A nice plot
ggplot(Deriv_dino_temp, aes(x = TEMP.med, y = derivative, 
                       color = Code_point_Libelle)) +
  # Derivative as points (function of median temperature)
  geom_point(aes(x = TEMP.med, y = derivative), 
             size = 1.5, alpha = .5) +
  geom_errorbar(aes(x = TEMP.med, ymin = lower, ymax = upper,
                    color = Code_point_Libelle)) +
  # Draw a line at 0 to separate accumulation from loss
  geom_line(aes(x = TEMP.med, y = 0), color = 'grey10', linewidth = .7) +
  labs(y = c(expression(paste("1st derivative of Dinophysis GAM (d"^"1",")"))),
       x = "Sea surface temperature (median of fortnight, °C)",
       title = "Dinophysis accumulation versus SST"
  ) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  scale_fill_discrete(type = pheno_palette12, guide = 'none') +
  theme_classic()

# It's kinda ugly but at least we got something
# Save it
# ggsave('DinoDeriv_temperature_12sites.tiff', dpi = 300, height = 225, width = 300,
#        units = 'mm', compression = 'lzw')

#### How do Dinophysis accumulation rates vary with Mesodinium cell counts ? ####

# Join the Dinophysis derivative and the Mesodinium model prediction
Deriv_dino_Meso <- left_join(response_pred_plot_Meso, Deriv_dino,
                             by = c('Code_point_Libelle', 'Day'))
# Just need to plot this

# Some aesthetics

# Relevel factor to order sampling sites
Deriv_dino_Meso <- Deriv_dino_Meso %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Bouzigues (a)', 'Diana centre'))

# a nice color palette
pheno_palette10 <- c('red3', 'orangered1', 'dodgerblue4', 'dodgerblue1', 
                     'chartreuse4', 'chartreuse2','goldenrod3', 'darkgoldenrod1',
                     'darkorchid1', 'deeppink2'
)

# And let's plot
# A nice plot
ggplot(Deriv_dino_Meso, aes(x = median.fit, y = derivative, 
                            color = Code_point_Libelle)) +
  # Derivative as points (function of median temperature)
  geom_point(aes(x = median.fit, y = derivative), 
             size = 1.5, alpha = .5) +
  # We add vertical AND horizontal error bars to figure uncertainty in both the
  # Dinophysis and Mesodinium models
  geom_errorbar(aes(xmin = lwrS, xmax = uprS,
                    ymin = lower, ymax = upper,
                    color = Code_point_Libelle)) +
  # Draw a line at 0 to separate accumulation from loss
  geom_line(aes(x = median.fit, y = 0), color = 'grey10', linewidth = .7) +
  labs(y = c(expression(paste("1st derivative of Dinophysis GAM (d"^"1",")"))),
       x = "Modelled Mesodinium cell count",
       title = "Dinophysis accumulation versus Mesodinium cell count"
  ) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette10, guide = 'none') +
  scale_fill_discrete(type = pheno_palette10, guide = 'none') +
  theme_classic()

# It's quite strange
# Save it
# ggsave('DinoDeriv_Mesodinium_10sites.tiff', dpi = 300, height = 225, width = 300,
#        units = 'mm', compression = 'lzw')

#### How do Dinophysis cell counts vary with Mesodinium cell counts ? ####

# Join the tables of Mesodinium and Dinophysis modelled counts
Modelled_Dino_Meso <- left_join(response_pred_plot_Meso, response_pred_plot_Dino,
                             by = c('Code_point_Libelle', 'Day'),
                             suffix = c('.meso','.dino'))
# Just need to plot this

# Some aesthetics

# Relevel factor to order sampling sites
Modelled_Dino_Meso <- Modelled_Dino_Meso %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Bouzigues (a)', 'Diana centre'))

# a nice color palette (only 10 sites because Mesodinium)
pheno_palette10 <- c('red3', 'orangered1', 'dodgerblue4', 'dodgerblue1', 
                     'chartreuse4', 'chartreuse2','goldenrod3', 'darkgoldenrod1',
                     'darkorchid1', 'deeppink2'
)

# And let's plot
# A nice plot
ggplot(Modelled_Dino_Meso, aes(x = median.fit.meso, y = median.fit.dino, 
                            color = Code_point_Libelle)) +
  # Derivative as points (function of median temperature)
  geom_point(aes(x = median.fit.meso, y = median.fit.dino), 
             size = 1.5, alpha = .5) +
  # We add vertical AND horizontal error bars to figure uncertainty in both the
  # Dinophysis and Mesodinium models
  geom_errorbar(aes(xmin = lwrS.meso, xmax = uprS.meso,
                    ymin = lwrS.dino, ymax = uprS.dino,
                    color = Code_point_Libelle)) +
  # Labels
  labs(y = "Modelled Dinophysis cell count",
       x = "Modelled Mesodinium cell count",
       title = "Dinophysis versus Mesodinium cell counts (modelled)"
  ) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette10, guide = 'none') +
  scale_fill_discrete(type = pheno_palette10, guide = 'none') +
  theme_classic()

# It's quite nice! Reminds me of Lotka-VOlterra predator/prey dynamics
# Save it
# ggsave('Dino_vs_Meso_modelled_10sites.tiff', dpi = 300, height = 225, width = 300,
#        units = 'mm', compression = 'lzw')
