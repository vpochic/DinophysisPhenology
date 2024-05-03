#### Phenology analysis - Dinophytes and diatoms REPHY
## V. POCHIC
# 2024-05-03

#### Packages ####
library(tidyverse)
library(ggplot2)
library(viridis)
library(see)
library(ggnewscale)
library(ghibli)
library(cmocean)

#### Import data ####

# Table with phytoplankton counts
Table_phyto_taxon <- read.csv2('Table1_phyto_taxon.csv', fileEncoding = "ISO-8859-1")

# Select only stations of interest from 2007 --> 2022
Table_phyto_select <- filter(Table_phyto_taxon, Code_point_Libelle %in% 
                           c(# Baie de Seine
                             'Antifer ponton pétrolier', 'Cabourg',
                             # Bretagne Sud
                             'Men er Roue', 'Ouest Loscolo',
                             # Pertuis charentais
                             'Auger', 'Le Cornard',
                             # Arcachon
                             'Arcachon - Bouée 7', 'Teychan bis',
                             # Mediterranée
                             'Parc Leucate 2', 'Bouzigues (a)', 'Sète mer',
                             'Diana centre')
) %>%
  filter(Year >= 2007 & Year <= 2022)

# To regroup all taxa in the group 'Dinophyceae', we need the magic list
Phylum_Class <- read.csv2('Liste_phylum.classe_REPHY_accents.csv', header = TRUE,
                          fileEncoding = 'ISO-8859-1')

# And left join
Table_phyto_select <- left_join(Table_phyto_select, Phylum_Class, by = 'Taxon',
                                suffix = c('',''))

#### More zeros ####

### The goal here is to create a table where all absences of detection are 
# listed as zeros
Table_phylum_zeros <- Table_phyto_select %>%
  # Eliminate the Taxon column
  select(-c('Taxon')) %>%
  # Only FLORTOT
  filter(Code.parametre == 'FLORTOT') %>%
  # Group by the variables defining the sample and the phylum
  group_by(Code_point_Libelle, Year, Month, Date, ID.interne.passage, 
           Phylum.Classe) %>%
  # Add all different values from a same phylum (for example, 200 Chaetoceros +
  # 100 Rhizosolenia = 300 Bacillariophyceae)
  summarise(Comptage = sum(Comptage), .groups = 'keep') %>%
  # And pivot wider with the Phyla
  pivot_wider(names_from = Phylum.Classe, values_from = Comptage, values_fill = 0)

# Nickel chrome

# Save Table_phylum_zeros so we don't have to re-run the model to plot again
# write.csv2(Table_phylum_zeros, 'Table_phylum_zeros.csv', row.names = FALSE,
#            fileEncoding = "ISO-8859-1")


### Clean the environment
rm(Table_phyto_taxon)
rm(Table_phyto_select)


#### Seasonality ####

### Create a seasonality dataset (with the calendar day/julian day variable)
Season_phylum <- Table_phylum_zeros %>%
  # create a calendar day variable
  mutate(Day = as.numeric(yday(Date))) %>%
  mutate(Date = ymd(Date)) %>%
  # Computing a 'fortnight' variable that matches the sampling frequency
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2)) %>%
  #### Note that for the 2 sites in Arcachon, counts are sometimes done in
  # 100 mL, but we can't really know which ones. So we stick with the cells.L-1
  # values.
  # converting the site to factor for the model
  mutate(Code_point_Libelle = as.factor(Code_point_Libelle))


# Checking the homogeneity of the sampling between years
# We make a histogram of samples for the 17 years (2006-2022)

hist(Season_phylum$Year, breaks = 16)
# Approximately homogeneous around 300 per year (all sites together), a shifty
# bit in the first year (weird plotting bug, not found in the raw data), a bit 
# less in 2020 (15-20) due to Covid.

# Homogeneity across the year, for the 26 fortnights in a year 
# (fortnightly sampling)
hist(Season_phylum$Day, breaks = 26)
# Approximately OK, except for the last fortnight of the year 
# (Christmas + NY's eve)

# This seems quite fine!!!

#### Max Dino abundance ####

Max_Dinophytes <- Season_phylum %>%
  group_by(Code_point_Libelle, Year) %>%
  # Select only the max in each year
  filter(Dinophyceae == max(Dinophyceae)) %>%
  filter(Year <= 2022) %>%
  # Filter out 'max' values that are 0s because no Dinophyte was observed 
  # that year (unlikely)
  filter(Dinophyceae != 0) %>%
  ungroup() %>%
  group_by(Code_point_Libelle)

# Plot aesthetics
pheno_palette12 <- c('red3', 'orangered1', 'dodgerblue4', 'dodgerblue1', 
                     'chartreuse4', 'chartreuse2','goldenrod3', 'darkgoldenrod1',
                     'darkorchid4', 'darkorchid1', 'firebrick1', 'deeppink2'
                     )

# Re-level factors
Max_Dinophytes <- Max_Dinophytes %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

### Violin plots of Dinophysis maxima per year (weighted by 'true_count')
ggplot(data = Max_Dinophytes, aes(x = Code_point_Libelle, y = Day)) +
  # geom_violin(aes(x = Code_point_Libelle, y = Day, 
  #                           fill = Code_point_Libelle,
  #                           color = Code_point_Libelle,
  #                 weight = true_count), alpha = .4, lwd = .75) +
  geom_point(aes(x = Code_point_Libelle, y = Day, 
                 fill = Code_point_Libelle,
                 color = Code_point_Libelle,
                 size = log10(Dinophyceae)), alpha = .4,
                 position = position_jitter(width = .05)) +
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  scale_fill_discrete(type = pheno_palette12, guide = 'none') +
  scale_size_continuous(guide = 'none') +
  theme_classic()

#### Max Diatom abundance ####

Max_Diatoms <- Season_phylum %>%
  group_by(Code_point_Libelle, Year) %>%
  # Select only the max in each year
  filter(Bacillariophyceae == max(Bacillariophyceae)) %>%
  filter(Year <= 2022) %>%
  # Filter out 'max' values that are 0s because no diatom was observed 
  # that year (unlikely)
  filter(Bacillariophyceae != 0) %>%
  ungroup() %>%
  group_by(Code_point_Libelle)

# Plot aesthetics
pheno_palette12 <- c('red3', 'orangered1', 'dodgerblue4', 'dodgerblue1', 
                     'chartreuse4', 'chartreuse2','goldenrod3', 'darkgoldenrod1',
                     'darkorchid4', 'darkorchid1', 'firebrick1', 'deeppink2'
)

# Re-level factors
Max_Diatoms <- Max_Diatoms %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

### Violin plots of Dinophysis maxima per year (weighted by 'true_count')
ggplot(data = Max_Diatoms, aes(x = Code_point_Libelle, y = Day)) +
  # geom_violin(aes(x = Code_point_Libelle, y = Day, 
  #                           fill = Code_point_Libelle,
  #                           color = Code_point_Libelle,
  #                 weight = true_count), alpha = .4, lwd = .75) +
  geom_point(aes(x = Code_point_Libelle, y = Day, 
                 fill = Code_point_Libelle,
                 color = Code_point_Libelle,
                 size = log10(Bacillariophyceae)), alpha = .4,
             position = position_jitter(width = .05)) +
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  scale_fill_discrete(type = pheno_palette12, guide = 'none') +
  scale_size_continuous(guide = 'none') +
  theme_classic()

#### Seasonality of Dinophyte presence ####

Season_plot <- Season_phylum %>%
  group_by(Code_point_Libelle, Year) %>%
  filter(Year <= 2022) %>%
  # Filter out zero values (the violin plot doesn't need them)
  filter(Dinophyceae != 0) %>%
  ungroup() %>%
  group_by(Code_point_Libelle)

# Plot aesthetics
pheno_palette12 <- c('red3', 'orangered1', 'dodgerblue4', 'dodgerblue1', 
                     'chartreuse4', 'chartreuse2','goldenrod3', 'darkgoldenrod1',
                     'darkorchid4', 'darkorchid1', 'firebrick1', 'deeppink2'
                     )

# Re-level factors
Season_plot <- Season_plot %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

### Violin plot of Dinophyceae observations (weighted by 'true_count'), plus
# all Dinophyceae observations as points
ggplot(data = Season_plot, aes(x = Code_point_Libelle, y = Day)) +
  geom_violin(aes(x = Code_point_Libelle, y = Day,
                  fill = Code_point_Libelle,
                  color = Code_point_Libelle,
                  weight = log10(Dinophyceae)),
              alpha = .4, lwd = .75) +
  geom_point(data = Season_plot, aes(x = Code_point_Libelle, y = Day, 
                                  fill = Code_point_Libelle, 
                                  size = log10(Dinophyceae)), 
             alpha = .4,
             shape = 21,
             color = 'black',
             position = position_jitter(width = .05)) +
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  scale_fill_discrete(type = pheno_palette12, guide = 'none') +
  scale_size_continuous(guide = 'none') +
  theme_classic()


#### Seasonality of Dinophyte presence + maxima ####

Season_plot <- Season_phylum %>%
  group_by(Code_point_Libelle, Year) %>%
  filter(Year <= 2022) %>%
  # Filter out zero values (the violin plot doesn't need them)
  filter(Dinophyceae != 0) %>%
  ungroup() %>%
  group_by(Code_point_Libelle)

# Plot aesthetics
pheno_palette12 <- c('red3', 'orangered1', 'dodgerblue4', 'dodgerblue1', 
                     'chartreuse4', 'chartreuse2','goldenrod3', 'darkgoldenrod1',
                     'darkorchid4', 'darkorchid1', 'firebrick1', 'deeppink2'
)

# Re-level factors
Season_plot <- Season_plot %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

### Violin plot of Dinophyceae observations (weighted by log10(count)), plus
# maxima as points
ggplot(data = Season_plot, aes(x = Code_point_Libelle, y = Day)) +
  geom_violinhalf(aes(x = Code_point_Libelle, y = Day,
                  fill = Code_point_Libelle,
                  color = Code_point_Libelle
                  ,weight = log10(Dinophyceae)
                  ),
                  alpha = .4, lwd = .75) +
  geom_point(data = Max_Dinophytes, aes(x = Code_point_Libelle, y = Day, 
                 fill = Code_point_Libelle, size = log10(Dinophyceae)), 
             alpha = .4,
             shape = 21,
             color = 'black',
             position = position_jitter(width = .05)) +
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  scale_fill_discrete(type = pheno_palette12, guide = 'none') +
  scale_size_continuous(guide = 'none') +
  # Labels
  labs(title = "Phenology of Dinophyceae observations",
       y = "Day of the year",
       x = "Sampling site" 
       ) +
  # this flips the y and x axes so that violin plots are horizontal
  coord_flip() +
  # and this ensures that the first (= northernmost) sites are at the top of
  # the plot
  scale_x_discrete(limits = rev) +
  theme_bw()

# ggsave('Dinophytes_phenology_half-violin_weighted.tiff', width = 150, height = 300,
#        units = 'mm', compression = 'lzw')


#### Seasonality of diatom presence + maxima ####

# We can do exactly the same with Bacillariophyceae
Season_plot <- Season_phylum %>%
  group_by(Code_point_Libelle, Year) %>%
  filter(Year <= 2022) %>%
  # Filter out zero values (the violin plot doesn't need them)
  filter(Bacillariophyceae != 0) %>%
  ungroup() %>%
  group_by(Code_point_Libelle)

# Plot aesthetics
pheno_palette12 <- c('red3', 'orangered1', 'dodgerblue4', 'dodgerblue1', 
                     'chartreuse4', 'chartreuse2','goldenrod3', 'darkgoldenrod1',
                     'darkorchid4', 'darkorchid1', 'firebrick1', 'deeppink2'
)

# Re-level factors
Season_plot <- Season_plot %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

### Violin plot of Bacillariophyceae observations (weighted by log10(count)), plus
# maxima as points
ggplot(data = Season_plot, aes(x = Code_point_Libelle, y = Day)) +
  geom_violinhalf(aes(x = Code_point_Libelle, y = Day,
                      fill = Code_point_Libelle,
                      color = Code_point_Libelle
                      ,weight = log10(Bacillariophyceae)
  ),
  alpha = .4, lwd = .75) +
  geom_point(data = Max_Diatoms, aes(x = Code_point_Libelle, y = Day, 
                                        fill = Code_point_Libelle, 
                                     size = log10(Bacillariophyceae)), 
             alpha = .4,
             shape = 21,
             color = 'black',
             position = position_jitter(width = .05)) +
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  scale_fill_discrete(type = pheno_palette12, guide = 'none') +
  scale_size_continuous(guide = 'none') +
  # Labels
  labs(title = "Phenology of Bacillariophyceae observations",
       y = "Day of the year",
       x = "Sampling site" 
  ) +
  # this flips the y and x axes so that violin plots are horizontal
  coord_flip() +
  # and this ensures that the first (= northernmost) sites are at the top of
  # the plot
  scale_x_discrete(limits = rev) +
  theme_bw()

# ggsave('Diatoms_phenology_half-violin_weighted.tiff', width = 150, height = 300,
#        units = 'mm', compression = 'lzw')

