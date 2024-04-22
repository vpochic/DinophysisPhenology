#### Phenology analysis - Dinophysis REPHY
## V. POCHIC
# 2024-04-19

#### Packages ####
library(tidyverse)
library(ggplot2)
library(viridis)
library(see)
library(ggnewscale)

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

#### More zeros ####

### The goal here is to create a table where all absences of detection are 
# listed as zeros
Table_Dino_zeros <- Table_phyto_select %>%
  pivot_wider(names_from = Taxon, values_from = Comptage, values_fill = 0) %>%
  select(starts_with('Dinophysis') | 
           contains(c('Code.Region', 'Code_point_Libelle', 'Year', 'Month',
                      'Date', 'Code.parametre', 'SALI', 'TEMP')))

Table_Dino_zeros <- Table_Dino_zeros %>%
  # Getting rid of some unwelcome guests (who got onboard because the taxon
  # name contains 'sali')
  select(-(contains('Diplo')))

# Save Table_Dino_zeros so we don't have to re-run the model to plot again
# write.csv2(Table_Dino_zeros, 'Table_Dino_zeros2.csv', row.names = FALSE,
#            fileEncoding = "ISO-8859-1")


### Clean the environment
rm(Table_phyto_taxon)
rm(Table_phyto_select)


#### Seasonality ####

### Create a seasonality dataset (with the calendar day/julian day variable)
Season_Dino <- Table_Dino_zeros %>%
  # Only FLORTOT
  filter(Code.parametre == 'FLORTOT') %>%
  # create a calendar day variable
  mutate(Day = as.numeric(yday(Date))) %>%
  mutate(Date = ymd(Date)) %>%
  # Computing a 'fortnight' variable that matches the sampling frequency
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2)) %>%
  #### Note that for the 2 sites in Arcachon, counts are sometimes done in
  # 100 mL, so we need to filter out those counts
  # For this, we go and look for any count values IN ANY DINOPHYSIS TAXON that
  # cannot correspond to 10 mL counts (i.e., not multiples of 100)
  ### We replace these values with 0. Note that for a given date, we can have
  # one Dinophysis species counted in 10 mL and another in 100 mL
  ## That's why we have to proceed taxon by taxon
  mutate(across(.cols = c('Dinophysis', 'Dinophysis + phalacroma', 'Dinophysis acuta',
                          'Dinophysis acuminata', 'Dinophysis caudata', 'Dinophysis tripos',
                          'Dinophysis sacculus', 'Dinophysis fortii', 
                          'Dinophysis hastata + odiosa'), .fns = ~ ifelse(. %% 100 != 0, 0, .))) %>%
  # And creating some count variables for Dinophysis as a genus
  mutate(Dinophysis_genus = rowSums(across(contains('Dinophysis')))) %>%
  # create the log of abundance + 1
  mutate(log_c = log10(Dinophysis_genus+1)) %>%
  # create a 'true count' variable
  # this variable corresponds to the number of cells that were actually
  # counted by the operator (in 10 mL). This variable will follow a Poisson
  # distribution, contrary to Dinophysis_genus because the conversion from
  # 10 mL to 1L (*100) prevents some intermediate values (e.g., 150 cells.L-1)
  mutate(true_count = Dinophysis_genus/100) %>%
  # Filter out exceptionnaly high counts (>500 cells observed in 10 mL)
  # this represents 3 events (2 in Antifer and 1 in Cabourg)
  filter(true_count < 500) %>%
  # converting the site to factor for the model
  mutate(Code_point_Libelle = as.factor(Code_point_Libelle))


# Checking the homogeneity of the sampling between years
# We make a histogram of samples for the 17 years (2006-2022)

hist(Season_Dino$Year, breaks = 16)
# Approximately homogeneous around 300 per year (all sites together), a shifty
# bit in the first year (weird plotting bug, not found in the raw data), a bit 
# less in 2020 (15-20) due to Covid.

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
pheno_palette12 <- c('red3', 'orangered1', 'dodgerblue4', 'dodgerblue1', 
                     'chartreuse4', 'chartreuse2','goldenrod3', 'darkgoldenrod1',
                     'darkorchid4', 'darkorchid1', 'firebrick1', 'deeppink2'
                     )

# Re-level factors
Max_Dino <- Max_Dino %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

### Violin plots of Dinophysis maxima per year (weighted by 'true_count')
ggplot(data = Max_Dino, aes(x = Code_point_Libelle, y = Day)) +
  # geom_violin(aes(x = Code_point_Libelle, y = Day, 
  #                           fill = Code_point_Libelle,
  #                           color = Code_point_Libelle,
  #                 weight = true_count), alpha = .4, lwd = .75) +
  geom_point(aes(x = Code_point_Libelle, y = Day, 
                 fill = Code_point_Libelle,
                 color = Code_point_Libelle,
                 size = log_c), alpha = .4,
                 position = position_jitter(width = .05)) +
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  scale_fill_discrete(type = pheno_palette12, guide = 'none') +
  scale_size_continuous(guide = 'none') +
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

### Violin plot of Dinophysis observations (weighted by 'true_count'), plus
# all Dinophysis observations as points
ggplot(data = Season_plot, aes(x = Code_point_Libelle, y = Day)) +
  geom_violin(aes(x = Code_point_Libelle, y = Day,
                  fill = Code_point_Libelle,
                  color = Code_point_Libelle,
                  weight = true_count),
              alpha = .4, lwd = .75) +
  geom_point(data = Season_plot, aes(x = Code_point_Libelle, y = Day, 
                                  fill = Code_point_Libelle, size = log_c), 
             alpha = .4,
             shape = 21,
             color = 'black',
             position = position_jitter(width = .05)) +
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  scale_fill_discrete(type = pheno_palette12, guide = 'none') +
  scale_size_continuous(guide = 'none') +
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

### Violin plot of Dinophysis observations (weighted by true_count), plus
# maxima as points
ggplot(data = Season_plot, aes(x = Code_point_Libelle, y = Day)) +
  geom_violinhalf(aes(x = Code_point_Libelle, y = Day,
                  fill = Code_point_Libelle,
                  color = Code_point_Libelle
                  ,weight = true_count
                  ),
                  alpha = .4, lwd = .75) +
  geom_point(data = Max_Dino, aes(x = Code_point_Libelle, y = Day, 
                 fill = Code_point_Libelle, size = log_c), 
             alpha = .4,
             shape = 21,
             color = 'black',
             position = position_jitter(width = .05)) +
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  scale_fill_discrete(type = pheno_palette12, guide = 'none') +
  scale_size_continuous(guide = 'none') +
  # Labels
  labs(title = "Phenology of Dinophysis observations",
       y = "Day of the year",
       x = "Sampling site" 
       ) +
  # this flips the y and x axes so that violin plots are horizontal
  coord_flip() +
  # and this ensures that the first (= northernmost) sites are at the top of
  # the plot
  scale_x_discrete(limits = rev) +
  theme_bw()

# ggsave('Dinophysis_phenology_half-violin_unweighted.tiff', width = 150, height = 300,
#        units = 'mm', compression = 'lzw')


#### Adding a heatmap for temperature ####

Season_plot <- Season_Dino %>%
  group_by(Code_point_Libelle, Year) %>%
  # Filter out zero values (the violin plot doesn't need them)
  filter(Dinophysis_genus != 0) %>%
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

# Import the hydrological data
Table_hydro_fortnightly = read.csv2('Table_hydro_fortnightly.csv', header = TRUE,
                                    fileEncoding = 'ISO-8859-1')

# Create a daily table for plotting the heatmap
# We create a vector for 'Day of the year'
Daily_basis <- expand_grid(Day = seq(1,365))
# And one for 'Fortnight' that matches
Daily_basis$Fortnight = as.vector(c(rep(1:26, each = 14), 26))
# We add the site data
Daily_basis <- expand_grid(Daily_basis,
                           Code_point_Libelle = unique(Table_hydro_fortnightly$Code_point_Libelle))

# We arrange Daily_basis by Site
Daily_basis <- Daily_basis %>%
  group_by(Code_point_Libelle, Day) %>%
  arrange(Code_point_Libelle)

# Nickel chrome!

Table_hydro_daily <- left_join(Daily_basis, Table_hydro_fortnightly,
                               by = c('Code_point_Libelle', 'Fortnight'),
                               suffix  = c('','')) %>%
  filter(is.na(TEMP.med) == FALSE) %>%
  group_by(Day, Code_point_Libelle)

### Violin plot of Dinophysis observations (weighted by true_count), plus
# maxima as points, and with the heatmap for temperature
ggplot(data = Season_plot, aes(x = Code_point_Libelle, y = Day)) +
  # First part of the plot : the heatmap of temperature
  geom_tile(data = Table_hydro_daily,
            aes(x = Code_point_Libelle, y = Day, fill = TEMP.med
                                   )) +
  scale_fill_distiller(palette = 'RdBu') +
  # Labels
  labs(title = "Phenology of Dinophysis observations",
       y = "Day of the year",
       x = 'Sampling site',
       fill = 'Temperature (°C)'
  ) +
  # This resets the color scale for fill
  new_scale_fill() +
  # Second part of the plot : the violin plot
  geom_violinhalf(aes(x = Code_point_Libelle, y = Day,
                      fill = Code_point_Libelle,
                      color = Code_point_Libelle
                      ,weight = true_count
  ),
  alpha = .4, lwd = .75) +
  geom_point(data = Max_Dino, aes(x = Code_point_Libelle, y = Day, 
                                  fill = Code_point_Libelle, size = log_c), 
             alpha = .4,
             shape = 21,
             color = 'black',
             position = position_jitter(width = .05)) +
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  scale_fill_discrete(type = pheno_palette12, guide = 'none') +
  scale_size_continuous(guide = 'none') +
  # this flips the y and x axes so that violin plots are horizontal
  coord_flip() +
  # and this ensures that the first (= northernmost) sites are at the top of
  # the plot
  scale_x_discrete(limits = rev) +
  theme(
    panel.background = element_rect(fill = 'transparent', 
                                   linewidth = .3, 
                                   color = 'grey10'),
    legend.background = element_rect(linewidth = .5, color = 'grey10'),
    legend.title = element_text(size = 11, color = 'grey5'),
    legend.frame = element_rect(linewidth = .2, color = 'grey25'),
    legend.ticks = element_line(linewidth = .2, color = 'grey25'),
    legend.position = 'bottom'
  )

# Save that beautiful shit!!!!
ggsave('Dinophysis_phenology_half-violin_weighted_heatmap.tiff', width = 160, height = 320,
       units = 'mm', compression = 'lzw', dpi = 300)

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