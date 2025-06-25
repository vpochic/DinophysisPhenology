###### Generalised Linear Models of Dinophysis change rate ###
## V. POCHIC
# 2025-06-23

# The goal here is to examine the relationship between Dinophysis change rate
# (the derivative of the phenology GAM) and several selected environmental
# parameters

# Based on the tutorial on GLMs by Bede Davies
# https://bedeffinianrowedavies.com/statisticstutorials/gaussianglms

# Careful! We need data files computed with the following scripts:
# GAM_Dino_all_sites-mgcv.R (for the GAM models)
# Deriv_GAMs_DinoPhenology_REPHY.R (for the derivatives = acc. rates)
# Dino_phenology_heatmaps.R (for the environmental variables)

### Packages ####
library(tidymodels)
library(tidyverse)
library(vip)
library(ranger)
library(DescTools)
library(viridis)
library(ggpointdensity)
library(RColorBrewer)
library(ggnewscale)

### Import data ####
## Part 1: Dinophysis sampling data ####
# That's the file Season_Dino.csv with all sampling events in all sites for
# 2007-2022
Season_Dino <- read.csv2('Data/REPHY_outputs/Season_Dino_20250604.csv', header = TRUE, 
                         fileEncoding = 'ISO-8859-1')


Season_Dino_12sites <- filter(Season_Dino, Code_point_Libelle %in% 
                                c(# Pas de Calais
                                  'Point 1 Boulogne', 'At so',
                                  # Baie de Seine
                                  'Antifer ponton pétrolier', 'Cabourg',
                                  # Bretagne Nord
                                  'les Hébihens', 'Loguivy',
                                  # Bretagne Sud
                                  'Men er Roue', 'Ouest Loscolo',
                                  # Pertuis charentais
                                  'Auger', 'Le Cornard',
                                  # Arcachon
                                  'Arcachon - Bouée 7', 'Teychan bis')) %>%
  # Get rid of the one case in which Heure is ""
  # (Antifer, 2020/05/27, 0 Dinophysis)
  filter(Heure != "") %>%
  # Get a date-time variable by pasting together date and time
  mutate(DateTime = paste0(Date, ' ', Heure)) %>%
  # Get the Target date in the right format for corresponding with the stratif
  # dataset
  mutate(Target_Date = ymd_hms(DateTime)) %>%
  # And get the Date variable in the right format
  mutate(Date = ymd(Date))

## Part 2: Mesodinium sampling data ####
# That's the file Season_Meso.csv with all sampling events in all sites for
# 2007-2022
Season_Meso <- read.csv2('Data/REPHY_outputs/Season_Meso_20250604.csv', 
                         header = TRUE, fileEncoding = 'ISO-8859-1')


Season_Meso_12sites <- filter(Season_Meso, Code_point_Libelle %in% 
                                c(# Pas de Calais
                                  'Point 1 Boulogne', 'At so',
                                  # Baie de Seine
                                  'Antifer ponton pétrolier', 'Cabourg',
                                  # Bretagne Nord
                                  'les Hébihens', 'Loguivy',
                                  # Bretagne Sud
                                  'Men er Roue', 'Ouest Loscolo',
                                  # Pertuis charentais
                                  'Auger', 'Le Cornard',
                                  # Arcachon
                                  'Arcachon - Bouée 7', 'Teychan bis')) %>%
  # Get rid of the potential cases in which Heure is ""
  filter(Heure != "") %>%
  # Get a date-time variable by pasting together date and time
  mutate(DateTime = paste0(Date, ' ', Heure)) %>%
  # Get the Target date in the right format for corresponding with the stratif
  # dataset
  mutate(Target_Date = ymd_hms(DateTime)) %>%
  # And get the Date variable in the right format
  mutate(Date = ymd(Date))

# Same number of events as the Dinophysis tibble: seems ok

## Part 3: Pigment data ####
# We want to get the HPLC pigment data to look at the alloxanthin

Table_pigments <- read.csv2('Data/REPHY_outputs/Table_pigments_2007-2022.csv', 
                            header = TRUE, fileEncoding = 'ISO-8859-1')

# Only keep the 4 sites we will use for our second RF model
Table_pigments_select <- filter(Table_pigments, Code_point_Libelle %in% 
                                  c(# Baie de Seine
                                    'Antifer ponton pétrolier', 'Cabourg',
                                    # Bretagne Sud
                                    'Men er Roue', 'Ouest Loscolo')) %>%
  # And get the Date variable in the right format
  mutate(Date = ymd(Date))

# This dataset is a bit different from the others.

## Part 4: gam derivative data ####

# We import the dataframe of the derivative of the GAM.
# From the script Deriv_GAMs_DinoPhenology...
gam_Dino.d <- read.csv2('Data/GAM_outputs/Derivatives_GAM_Dino_12sites_20241203.csv', 
                        header = TRUE, fileEncoding = 'ISO-8859-1')

# Or with the multiyear data
gam_Dino.d_multiyear <- read.csv2('Data/GAM_outputs/gam_Dino_multiyear_deriv_20241204.csv', 
                                  header = TRUE, fileEncoding = 'ISO-8859-1') %>%
  # Put Date in date format
  mutate(Date = ymd(Date))

## Part 5: environmental data ####
### Hydrological parameters form the REPHY dataset
Table_hydro <- read.csv2('Data/REPHY_outputs/Table1_hydro_select.csv', header = TRUE, 
                         fileEncoding = "ISO-8859-1")

# Select only stations of interest from 2007 --> 2022
Table_hydro_select <- filter(Table_hydro, Code_point_Libelle %in% 
                               c(# Pas de Calais
                                 'Point 1 Boulogne', 'At so',
                                 # Baie de Seine
                                 'Antifer ponton pétrolier', 'Cabourg',
                                 # Bretagne Nord
                                 'les Hébihens', 'Loguivy',
                                 # Bretagne Sud
                                 'Men er Roue', 'Ouest Loscolo',
                                 # Pertuis charentais
                                 'Auger', 'Le Cornard',
                                 # Arcachon
                                 'Arcachon - Bouée 7', 'Teychan bis')
) %>%
  # Only after 2007
  filter(Year >= 2007) %>%
  # Only events for which Temperature and Salinity were properly recorded
  filter(is.na(TEMP) == FALSE) %>%
  filter(is.na(SALI) == FALSE) %>%
  # Change to date format
  mutate(Date = ymd(Date)) %>%
  mutate(Month = month(Date)) %>%
  mutate(Day = as.numeric(yday(Date))) %>%
  # Create a 'fortnight' variable to match the sampling frequency
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2))

# Checking for NAs
# Temperature
isNA_TEMP <- filter(Table_hydro_select, is.na(TEMP) == TRUE)
# Clean
# Salinity
isNA_SALI <- filter(Table_hydro_select, is.na(SALI) == TRUE)
# Clean
# Chlorophyll a
isNA_CHLOROA <- filter(Table_hydro_select, is.na(CHLOROA) == TRUE)
# Not clean. 733 obs with CHLOROA = NA, from 2007 to 2022, in several sites

## Part 6: Model data ####

### Stratification index from the GAMAR model ##
Table_stratif <- read.csv2('Data/Models/GAMAR/Outputs/Stratif_index_GAMAR_12sites.csv',
                           header = TRUE, fileEncoding = 'ISO-8859-1') %>%
  # Getting it in date format
  mutate(Target_Date = ymd_hms(Target_Date))

# The stratif file is much heavier because we have ALL THE VALUES FOR ALL THE
# DATES ON ALL SITES
# By joining our files and ditching what's superfluous, we should manage

### Atmospheric data from the ERA5 model ##

# These data are in 3 period formats (daily, weekly, fortnightly) and in 2 stats
# format (median over the period or mean over the period)

## Daily
# Mean
Table_era5_daily_mean <- read.csv2('Data/Models/ERA5/Outputs/era5_dataset_daily_mean.csv', 
                                   header = TRUE, fileEncoding = 'ISO-8859-1')
# Median
Table_era5_daily_median <- read.csv2('Data/Models/ERA5/Outputs/era5_dataset_daily_mean.csv', 
                                     header = TRUE, fileEncoding = 'ISO-8859-1')

## Weekly
# Mean
Table_era5_weekly_mean <- read.csv2('Data/Models/ERA5/Outputs/era5_dataset_weekly_mean.csv', 
                                    header = TRUE, fileEncoding = 'ISO-8859-1')
# Median
Table_era5_weekly_median <- read.csv2('Data/Models/ERA5/Outputs/era5_dataset_weekly_mean.csv', 
                                      header = TRUE, fileEncoding = 'ISO-8859-1')

## Fortnightly
# Mean
Table_era5_fortnight_mean <- read.csv2('Data/Models/ERA5/Outputs/era5_dataset_fortnight_mean.csv', 
                                       header = TRUE, fileEncoding = 'ISO-8859-1')
# Median
Table_era5_fortnight_median <- read.csv2('Data/Models/ERA5/Outputs/era5_dataset_fortnight_mean.csv', 
                                         header = TRUE, fileEncoding = 'ISO-8859-1')

#### Joining datasets ####

## First, let's join Season_Dino and stratif to see if this goes according to plan
# (1 sampling event = 1 row in the stratif dataset)
Table_data_RF <- left_join(Season_Dino_12sites, Table_stratif, 
                           by = c('Target_Date', 'Code_point_Libelle'),
                           suffix = c('',''))

# It works. It just works.

# Now let's join the hydrology table

Table_data_RF <- left_join(Table_data_RF, Table_hydro_select, 
                           by = c('Code_point_Libelle', 'Date'),
                           suffix = c('',''))

# Fantastic

## Third, the ERA5 data (only daily mean for now)
# events are associated by day and year (and sampling site of course)
Table_data_RF <- left_join(Table_data_RF, Table_era5_daily_mean, 
                           by = c('Year', 'Day', 'Code_point_Libelle'),
                           suffix = c('',''))

# Good.

# And now the tricky part: we project the derivative on every date based on the
# calendar day.
## Note: it would be even better to integrate the random effect of the year in
# the future (done)

# First, we remove duplicates in the GAM derivatives dataset (not the one with 
# every year) because Year isn't included
gam_Dino.d_distinct <- distinct(gam_Dino.d)

Table_data_RF_uniyear <- left_join(Table_data_RF, gam_Dino.d_distinct, 
                                   by = c('Code_point_Libelle', 'Day'), suffix = c('',''))

# For the one with 
Table_data_RF_multiyear <- left_join(Table_data_RF, gam_Dino.d_multiyear, 
                                     by = c('Code_point_Libelle', 'Date'), 
                                     suffix = c('','')) %>%
  # change the name of the derivative for consistency
  mutate(.derivative = deriv) %>%
  # Sampling site as factor and reorder
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis'))

# Let's see how much data we have in each site
ggplot(Table_data_RF_multiyear) +
  geom_histogram(aes(x = Code_point_Libelle, fill = Code_point_Libelle),
                 stat = 'count') +
  theme_classic()
# That's a nice rainbow. It seems we have between 250 (Antifer, At so) and 400
# (Auger, Le Cornard) data points. This seems ok.

### WARNING ### 2025/03/17
# There seems to be a big problem with TURB values (only 505 are not NA). It
# probably has to do with FNU to NTU conversion that was not done in the data
# curating scripts. We should correct that later.
#### Model building - 12 sites ####

# We will first build models based on the 12 sites at our disposal in the
# Atlantic-Channel region. Later, we will focus on only 4 sites for which we 
# have additional data. (But let's go 1 step at a time)

## Temperature ####

# The first effect we want to test is that of sea surface temperature.

# First, let's plot the data
pheno_palette12 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646')

ggplot(data = Table_data_RF_multiyear,
       aes(x = TEMP, y = .derivative, color = Code_point_Libelle)) +
  # Points
  geom_point(size = 2, alpha = .7) +
  # Add 2 "ghost points" to force the minimum y-axis range to -0.1;0.1
  geom_point(aes(x = 10, y = -.1), color = 'transparent', fill = 'transparent') +
  geom_point(aes(x = 10, y = .1), color = 'transparent', fill = 'transparent') +
  # Color scale (sampling sites)
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  # New color scale for temperature rug plot
  new_scale_color() +
  geom_rug(data = Table_data_RF_multiyear, aes(x = TEMP, y = NULL,
                                               color = TEMP)) +
  scale_color_distiller(palette = 'RdBu', direction = -1, guide = 'none') +
  # facets
  facet_wrap(facets = 'Code_point_Libelle', scales = 'free') +
  # Labels
  labs(x = 'Sea Surface Temperature (°C)', y = 'GAM derivative',) +
  # Theme
  theme_classic()

# This seems pretty fckin difficult
