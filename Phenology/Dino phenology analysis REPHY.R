#### Phenology analysis - Dinophysis REPHY
## V. POCHIC
# 2024-04-03

#### Packages ####
library(tidyverse)
library(ggplot2)
library(viridis)

##### Import data ####

Table_Dino_zeros <- read.csv2('Table_Dino_zeros.csv', header = TRUE, 
                              fileEncoding = "ISO-8859-1")


#### Seasonality ####

### Create a seasonality dataset (with the calendar day/julian day variable)
Season_Dino <- Table_Dino_zeros %>%
  # Only FLORTOT
  filter(Code.parametre == 'FLORTOT') %>%
  # create a calendar day variable
  mutate(Day = as.numeric(yday(Date))) %>%
  mutate(Date = ymd(Date)) %>%
  # create the log of abundance + 1
  mutate(log_c = log10(Dinophysis_genus+1)) %>%
  # create a 'true count' variable
  # this variable corresponds to the number of cells that were actually
  # counted by the operator (in 10 mL). This variable will follow a Poisson
  # distribution, contrary to Dinophysis_genus because the conversion from
  # 10 mL to 1L (*100) prevents some intermediate values (e.g., 150 cells.L-1)
  #### Note that for the 2 sites in Arcachon, counts are *probably* done in
  # 100 mL, so we need to specify that
  mutate(true_count = ifelse(
    # Si
    Code_point_Libelle == 'Arcachon - Bouée 7' |
      # Ou
      Code_point_Libelle == 'Teychan bis',
    # Alors
    Dinophysis_genus/10, 
    # Sinon
    Dinophysis_genus/100)) %>%
  # Filter out exceptionnaly high counts (>500 cells observed in 10 mL)
  # this represents 3 events (2 in Antifer and 1 in Cabourg)
  filter(true_count < 500) %>%
  # Filter some counts that do not correspond to the protocol (in Arcachon)
  filter(Dinophysis_genus %% 10 == 0) %>%
  # converting the site to factor for the model
  mutate(Code_point_Libelle = as.factor(Code_point_Libelle))


# Checking the homogeneity of the sampling between years
# We make a histogram of samples for the 17 years (2006-2022)

hist(Season_Dino$Year, breaks = 17)
# Approximately homogeneous around 25 per year, a bit more sampling in 2006
# (almost 35), a bit less in 2020 (15-20) due to Covid.

# Homogeneity across the year, for the 26 fortnights in a year 
# (fortnightly sampling)
hist(Season_Dino$Day, breaks = 26)
# Approximately OK, except for the last fortnight of the year 
# (Christmas + NY's eve)

# This seems quite fine!!!

#### Max Dino abundance ####

Max_Dino <- Season_Dino %>%
  group_by(Code_point_Libelle, Year) %>%
  # Select only the max in each year
  filter(Dinophysis_genus == max(Dinophysis_genus)) %>%
  filter(Year <= 2022) %>%
  # Filter out 'max' values that are 0s because no Dinophysis was observed 
  # that year
  filter(Dinophysis_genus != 0) %>%
  ungroup() %>%
  group_by(Code_point_Libelle)

# Plot aesthetics
pheno_palette13 <- c('red3', 'orangered1', 'chocolate4', 'chocolate',
                     'dodgerblue4', 'dodgerblue1', 'chartreuse4', 'chartreuse2',
                     'goldenrod3', 'darkgoldenrod1',
                     'darkorchid4', 'darkorchid1', 'deeppink2')

# Re-level factors
Max_Dino <- Max_Dino %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Diana centre'))

### Violin plots of Dinophysis maxima per year (weighted by 'true_count')
ggplot(data = Max_Dino, aes(x = Code_point_Libelle, y = Day)) +
  # geom_violin(aes(x = Code_point_Libelle, y = Day, 
  #                           fill = Code_point_Libelle,
  #                           color = Code_point_Libelle,
  #                 weight = true_count), alpha = .4, lwd = .75) +
  geom_point(aes(x = Code_point_Libelle, y = Day, 
                 fill = Code_point_Libelle,
                 color = Code_point_Libelle), alpha = .4, size = 4,
             position = position_jitter(width = .05)) +
  scale_color_discrete(type = pheno_palette13, guide = 'none') +
  scale_fill_discrete(type = pheno_palette13, guide = 'none') +
  theme_classic()

#### Seasonality of Dino presence ####

Season_plot <- Season_Dino %>%
  group_by(Code_point_Libelle, Year) %>%
  filter(Year <= 2022) %>%
  # Filter out zero values (the violin plot doesn't need them)
  filter(Dinophysis_genus != 0) %>%
  ungroup() %>%
  group_by(Code_point_Libelle)

# Plot aesthetics
pheno_palette13 <- c('red3', 'orangered1', 'chocolate4', 'chocolate',
                     'dodgerblue4', 'dodgerblue1', 'chartreuse4', 'chartreuse2',
                     'goldenrod3', 'darkgoldenrod1',
                     'darkorchid4', 'darkorchid1', 'deeppink2')

# Re-level factors
Season_plot <- Season_plot %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Diana centre'))

### Violin plot of Dinophysis observations (weighted by 'true_count'), plus
# all Dinophysis observations as points as points
ggplot(data = Season_plot, aes(x = Code_point_Libelle, y = Day)) +
  geom_violin(aes(x = Code_point_Libelle, y = Day,
                  fill = Code_point_Libelle,
                  color = Code_point_Libelle,
                  weight = true_count),
              alpha = .4, lwd = .75) +
  geom_point(data = Season_plot, aes(x = Code_point_Libelle, y = Day, 
                                  fill = Code_point_Libelle), alpha = .4, size = 2.2,
             shape = 21,
             color = 'black',
             position = position_jitter(width = .05)) +
  scale_color_discrete(type = pheno_palette13, guide = 'none') +
  scale_fill_discrete(type = pheno_palette13, guide = 'none') +
  theme_classic()


#### Seasonality of Dino presence + maxima ####

Season_plot <- Season_Dino %>%
  group_by(Code_point_Libelle, Year) %>%
  filter(Year <= 2022) %>%
  # Filter out zero values (the violin plot doesn't need them)
  filter(Dinophysis_genus != 0) %>%
  ungroup() %>%
  group_by(Code_point_Libelle)

# Plot aesthetics
pheno_palette13 <- c('red3', 'orangered1', 'chocolate4', 'chocolate',
                     'dodgerblue4', 'dodgerblue1', 'chartreuse4', 'chartreuse2',
                     'goldenrod3', 'darkgoldenrod1',
                     'darkorchid4', 'darkorchid1', 'deeppink2')

# Re-level factors
Season_plot <- Season_plot %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Diana centre'))

### Violin plot of Dinophysis observations (unweighted), plus
# maxima as points
ggplot(data = Season_plot, aes(x = Code_point_Libelle, y = Day)) +
  geom_violin(aes(x = Code_point_Libelle, y = Day,
                  fill = Code_point_Libelle,
                  color = Code_point_Libelle),
                  alpha = .4, lwd = .75) +
  geom_point(data = Max_Dino, aes(x = Code_point_Libelle, y = Day, 
                 fill = Code_point_Libelle), alpha = .4, size = 3,
             shape = 21,
             color = 'black',
             position = position_jitter(width = .05)) +
  scale_color_discrete(type = pheno_palette13, guide = 'none') +
  scale_fill_discrete(type = pheno_palette13, guide = 'none') +
  theme_bw()

# ggsave('Dinophysis_phenology_violin.tiff', width = 300, height = 150, 
#        units = 'mm', compression = 'lzw')

#### Test : Fragilidium ####

# Import complete taxon data
# Table with phytoplankton counts
Table_phyto_taxon <- read.csv2('Table1_phyto_taxon.csv', fileEncoding = "ISO-8859-1")

# Select only stations of interest from 2006 --> 2022
Table_phyto_OL <- filter(Table_phyto_taxon, Code_point_Libelle %in% 
                           c(# Nord-Pas de Calais
                             'Point 1 Boulogne', 'At so', 
                             # Baie de Seine
                             'Antifer ponton pétrolier', 'Cabourg',
                             # Bretagne Sud
                             'Men er Roue', 'Ouest Loscolo',
                             # Pertuis charentais
                             'Auger', 'Le Cornard',
                             # Arcachon
                             'Arcachon - Bouée 7', 'Teychan bis',
                             # Mediterranée
                             'Parc Leucate 2', 'Bouzigues (a)', 'Diana centre')
) %>%
  filter(Year >= 2006 & Year <= 2022)

test_fragi <- Table_phyto_OL %>%
pivot_wider(names_from = Taxon, values_from = Comptage, values_fill = 0) %>%
  select(contains('Fragilidium') | 
           contains(c('Code.Region', 'Code_point_Libelle', 'Year', 'Month',
                      'Date', 'Code.parametre', 'SALI', 'TEMP'))) %>%
  filter(Fragilidium != 0)

# Not many observations of Fragilidium