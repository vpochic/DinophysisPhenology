#### Script: analysing Continuous Plankton Recorder (CPR) data ###
# Author: V. POCHIC 
# Last modif: 2026/02/20

## General information ####

## Description: 
# This script takes a specific CPR dataset (found here: https://doi.mba.ac.uk/data/3428) 
# and uses it to analyse the phenology of Dinophysis in Europe since the 1950s

## Files required:
# Data/CPR/CPR_DREAMproject_Data_Dinophysis_03022025.csv
# (CPR data on Dinophysis)
# Data/CPR/CPR_regions_ids.csv
# Data/CPR/CPR_regions_polygon.csv
# (Small tables used for defining regions in the analysis)
# Data/REPHY_outputs/Season_Dino.csv
# (Only used to get the lat/lon coordinates of REPHY sampling sites, for plotting)

## Outputs:
## Figures:
# Figure 1 A and B - Seasonality of Dinophysis in European waters
# Figure S4 A and B - Sampling effort in the CPR dataset

### Required packages ####

library(tidyverse)
library(viridis)
library(ggnewscale)
library(deeptime)
library(mapdata)
library(spatstat.univar)

####-----------------------------------------------------------------------####

#### Import data ####

# CPR data
CPR_data <- read.csv('Data/CPR/CPR_DREAMproject_Data_Dinophysis_03022025.csv',
                      header = TRUE, fileEncoding = 'ISO-8859-1')

# A world map
Worldmap <- map_data('worldHires')

#### 1st step: Plot the data on a map ####

# Plotting a map with Dinophysis absence data in black and presence in red

ggplot() +
  # Plotting the land
  geom_polygon(data = Worldmap, aes(x = long, y = lat, group = group), 
                        fill = "grey80", color = 'gray10', linewidth = .25)+
  coord_fixed(xlim=c(-25,
                     max(CPR_data$Longitude)), 
              ylim=c(min(CPR_data$Latitude),
                     max(CPR_data$Latitude)), ratio=1)+
  # Plotting the samples without Dinophysis (many samples have a weird value of
  # 1e-10 for Dinophysis, we will consider it as an absence)
  geom_point(data = subset(CPR_data, Dinophysis_Total < 1), 
             aes(x = Longitude, y = Latitude),
             color = 'black', size = .5, alpha = .8) +
  # Plotting the samples with Dinophysis
  geom_point(data = subset(CPR_data, Dinophysis_Total > 1), 
             aes(x = Longitude, y = Latitude),
             color = 'red', size = .5, alpha = .8) +
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)',
       title = "CPR data in the NE Atlantic")+
  theme_classic()

# Nice one!

#### Bin the data by coordinate tiles ####

# Our objective here is to bin our data in "tiles" of 0.5*0.5 degrees of lon and
# lat.
# For this, we will create a lat.bin and a lon.bin variables

CPR_data_binned <- CPR_data %>%
  mutate(lat.bin = 
           # The reasoning here is: if the first digit of the latitude is
           # between 0 and 4, we will round it to .25 ; if it's between 5 and 9,
           #we'll round it to .75 IN POSITIVE OR NEGATIVE!!!
           ifelse(Latitude >= 0,
                  # if positive
             ifelse(
             firstdigit(Latitude-trunc(Latitude, digits = 0)) <= 4,
                  trunc(Latitude, digits = 0) + .25,
                  trunc(Latitude, digits = 0) + .75),
             # if negative
             ifelse(
               firstdigit(Latitude-trunc(Latitude, digits = 0)) <= 4,
               trunc(Latitude, digits = 0) - .25,
               trunc(Latitude, digits = 0) - .75))) %>%
  
  # Same thing for longitude
  mutate(lon.bin = 
           ifelse(Longitude >= 0,
                  # if positive
                  ifelse(
                    firstdigit(Longitude-trunc(Longitude, digits = 0)) <= 4,
                    trunc(Longitude, digits = 0) + .25,
                    trunc(Longitude, digits = 0) + .75),
                  # if negative
                  ifelse(
                    firstdigit(Longitude-trunc(Longitude, digits = 0)) <= 4,
                    trunc(Longitude, digits = 0) - .25,
                    trunc(Longitude, digits = 0) - .75)))

# Now, we can summarise the dataset by grouping it thanks to our new bins
CPR_data_sum <- CPR_data_binned %>%
  # Create a dummy variable that identifies Dinophysis positive counts
  # (0 = no Dinophysis in sample, 1 = 1 or more Dinophysis in sample)
  mutate(Dinodummy = ifelse(Dinophysis_Total > 1, 1, 0)) %>%
  # Grouping by the bin variables
  group_by(lat.bin, lon.bin) %>%
  # Summarise the dataset as we wish to
  summarise(
    # Total sum of Dinophysis counted in a bin
    sumDino = sum(Dinophysis_Total),
    # Number of samples with Dinophysis in a bin
    nsamples_Dino = sum(Dinodummy),
    # Total number of samples in a bin
    nsamples = n(),
    .groups = 'keep')

# Now, for each bin, we can calculate a proportion of samples containing 
# Dinophysis, and a mean number of Dinophysis per sample
CPR_data_sum <- CPR_data_sum %>%
  mutate(prop_pos = nsamples_Dino/nsamples) %>%
  mutate(mean_Dino = sumDino/nsamples)



#### Plotting a heatmap of Dinophysis presence ####

# Now the map will display the proportion of "Dinophysis-positive" samples with
# a color scale

# First, the map
ggplot() +
  # Plotting the tiles of longitude and latitude depending on Dinophysis presence
  geom_rect(data = subset(CPR_data_sum, nsamples>=3),
            aes(xmin = lon.bin-.25, xmax = lon.bin+.25,
                ymin = lat.bin-.25, ymax = lat.bin+.25,
                fill = log(prop_pos+1))) +
  scale_fill_viridis(
    limits = c(0, log(0.5+1)),
    breaks = c(log(0+1), log(0.1+1), log(0.25+1), (log(0.5+1))),
    labels = c('0', '0.1', '0.25', '0.5')) +
  # Plotting the land
  geom_polygon(data = Worldmap, aes(x = long, y = lat, group = group), 
               fill = "grey80", color = 'gray10', linewidth = .25)+
  
  # Limits of the map
  coord_fixed(xlim=c(-18,
                     18), 
              ylim=c(35,
                     65),
              ratio = 1.15)+
  # Labels
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)',
       title = 'Dinophysis in the NE Atlantic (CPR data)',
       fill = 
         'Dinophysis presence
(proportion of samples,
log scale)'
  )+
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.frame = element_rect(color = 'black'),
        legend.ticks = element_line(color = 'black'),
        legend.background = element_rect(fill = 'white', color = 'black'))

# Some areas seem to be visited more often than others by that good ol' Dinophysis.
# May this just be due to differences in sampling effort?

#### Plotting a heatmap of sampling effort ####

# Same map but with the sampling effort by tile (number of samples)

ggplot() +
  # Plotting the tiles of longitude and latitude depending on sampling effort
  geom_rect(data = CPR_data_sum,
            aes(xmin = lon.bin-.25, xmax = lon.bin+.25,
                ymin = lat.bin-.25, ymax = lat.bin+.25, # width = .5, height = .5,
                fill = log10(nsamples))) +
  scale_fill_viridis(option = 'plasma') +
  # Plotting the land
  geom_polygon(data = Worldmap, aes(x = long, y = lat, group = group), 
               fill = "grey80", color = 'gray10', linewidth = .25)+
  coord_fixed(xlim=c(-25,
                     max(CPR_data$Longitude)), 
              ylim=c(min(CPR_data$Latitude),
                     max(CPR_data$Latitude)))+
  # Labels
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)',
       title = "CPR sampling in the NE Atlantic", 
       fill = 'Number of samples (log10)') +
  theme_bw() +
  theme(legend.position = 'bottom')

# Sampling effort probably has something to do with the Dino "hotspots", but not
# entirely for sure.

#### Regionalising ####

# Based on the map of Dinophysis presence, we identify 4 hotspots of
# Dinophysis in the NE Atlantic: the Iberian coast (Galicia and Portugal),
# the Celtic seas, the Eastern North Sea (Scotland) and the Western North Sea
# (Germany, Denmark, Norway)

# We will attribute region numbers to these regions to study more specifically
# Dinophysis in each of them.

CPR_data_regions <- CPR_data_binned %>%
  # Define regions
  mutate(region = ifelse(
    # Iberian Coast (#1)
    lat.bin <= 44.5 & lat.bin >= 36.5 & lon.bin <= -8.5 & lon.bin >= -11,
    1,
    ifelse(
      # Celtic seas (#2)
      lat.bin <= 52 & lat.bin >= 48 & lon.bin <= -2.5 & lon.bin >= -9,
      2,
      ifelse(
        # Western North Sea (#3)
        lat.bin <= 60 & lat.bin >= 52.5 & lon.bin <= 1.5 & lon.bin >= -3.5,
        3,
        ifelse(
          # Eastern North Sea (#4)
          lat.bin <= 60 & lat.bin >= 52.5 & lon.bin <= 14.5 & lon.bin >= 3.5,
          4,
          # Else, #0
          0))))) %>%
  # Create an appropriate date variable
  mutate(Date = paste(Year, Month, Day, sep = '-')) %>%
  mutate(Date = ymd(Date)) %>%
  # Get a day of the year variable
  mutate(DOY = yday(Date)) %>%
  # Creating "week" and "fortnight" variables
  mutate(Week = week(Date)) %>%
  # 'ceiling' takes the upper integer of a decimal number
  mutate(Fortnight = ceiling(Week/2))

# Do the regions correspond to what we want?
ggplot() +
  # Plotting the tiles of longitude and latitude depending on the region
  geom_rect(data = CPR_data_regions,
            aes(xmin = lon.bin-.25, xmax = lon.bin+.25,
                ymin = lat.bin-.25, ymax = lat.bin+.25, # width = .5, height = .5,
                fill = as_factor(region))) +
  scale_fill_discrete() +
  # Plotting the land
  geom_polygon(data = Worldmap, aes(x = long, y = lat, group = group), 
               fill = "grey80", color = 'gray10', linewidth = .25)+
  coord_fixed(xlim=c(-25,
                     max(CPR_data$Longitude)), 
              ylim=c(min(CPR_data$Latitude),
                     max(CPR_data$Latitude)))+
  # Labels
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)',
       title = "Regionalization of the NE Atlantic")+
  theme_bw() +
  theme(legend.position = 'bottom')

# Ugly as hell, but it shows our regionalisation works.

# We will want to plot the regions as polygons on top of the Dinophysis presence
# map

####-----------------------------------------------------------------------####
#### Maps of CPR data with regions ####

#--- Import data ---#

# Get the data for REPHY sampling sites coordinates
Season_Dino <- read.csv2('Data/REPHY_outputs/Season_Dino.csv', header = TRUE, 
                         fileEncoding = 'ISO-8859-1')

# Summarise coordinates by site
Table_sites <- Season_Dino %>%
  group_by(Code_point_Libelle) %>%
  summarise(# Longitude
    Longitude.mean = mean(Longitude),
    # Latitude
    Latitude.mean = mean(Latitude),
    .groups = 'keep') %>%
  # Relevel sites to make them appear in the desired order
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

# Get data from 2 little tibbles for the defined regions, to know which info
# to pass onto geom_polygon()
# Polygon coordinates
regions_polygon <- read.csv2('Data/CPR/CPR_regions_polygon.csv', header = TRUE,
                             fileEncoding = 'ISO-8859-1')
# And a georeferenced tibble to plot the id of each region
regions_id <- read.csv2('Data/CPR/CPR_regions_ids.csv', header = TRUE,
                        fileEncoding = 'ISO-8859-1')

#--- Plotting time ---#

# Define nice color palettes

# Color palette for the REPHY sampling sites
pheno_palette16 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646',
                     '#642C3A', '#DEB1CC', '#FC4D6B', '#791D40')

# Color palette for the CPR-defined regions
palette_regions4 <- c('orange1', 'royalblue2', 'red4', 'forestgreen')

# Plot the damn map

# The same bit of script can be used to generate 2 maps: the one for Dinophysis
# presence, and one for sampling effort of the CPR. To generate the second one,
# replace in the ggplot command the parts where alternative options are indicated
# as "#comments"
ggplot() +
  # Plotting the tiles of longitude and latitude depending on Dinophysis presence
  geom_rect(data = subset(CPR_data_sum, nsamples>=3),
            aes(xmin = lon.bin-.25, xmax = lon.bin+.25,
                ymin = lat.bin-.25, ymax = lat.bin+.25,
                fill = 
                  log(prop_pos+1) # Dinophysis presence version
                  # log10(nsamples) # Sampling effort version
                  )) +
  # adjust the limits of the color scale so it is compatible with the second
  # plot. No legend here.
  scale_fill_viridis(
    # option = 'plasma' # for sampling effort version
        # All 3 options below are for Dinophysis presence version
    limits = c(0, log(0.5+1)), #
    breaks = c(log(0+1), log(0.1+1), log(0.25+1), (log(0.5+1))), #
    labels = c('0', '0.1', '0.25', '0.5') #
    ) + 
  # Plotting the land
  geom_polygon(data = Worldmap, aes(x = long, y = lat, group = group), 
               fill = "grey80", color = 'gray10', linewidth = .25)+
  # Plotting the limits of the regions
  geom_polygon(data = regions_polygon, aes(x = lon, y = lat, group = region,
                                           color = as_factor(region)), 
               fill = "transparent", linewidth = 1) +
  # Plotting the id of each region
  geom_text(data = regions_id, aes(x = lon, y = lat, color = as_factor(region),
                                   label = region), size = 6)+
  scale_color_discrete(type = palette_regions4, guide = 'none') +
  new_scale_color() +
  
  ### Plotting REPHY sampling sites
  geom_point(data = Table_sites, aes(x = Longitude.mean, y = Latitude.mean, 
                                     color = Code_point_Libelle), size = 2.5)+
  scale_color_discrete(type = pheno_palette16, guide = 'none') +
  
  # Limits of the map
  coord_fixed(xlim=c(-18,
                     18), 
              ylim=c(35,
                     65),
              ratio = 1.15)+
  # Labels
  labs(title = 'A',
       y = 'Latitude (degrees)', x = 'Longitude (degrees)',
       fill = 
         'Dinophysis presence
(proportion of samples,
log scale)'
#        'Sampling effort
# (log10[number of samples])'
  )+
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.frame = element_rect(color = 'black'),
        legend.ticks = element_line(color = 'black'),
        legend.background = element_rect(fill = 'white', color = 'black'))

## Save the map (Dinophysis presence version)
# ggsave('Plots/CPR/Fig1A_CPR_Dinophysis_map.tiff', height = 190, width = 160,
#        units = 'mm', dpi = 300, compression = 'lzw')

## Save the map (Sampling effort version)
# ggsave('Plots/CPR/FigS4A_CPR_sampling_effort_map.tiff', height = 190,
# width = 160, units = 'mm', dpi = 300, compression = 'lzw')


### Hovmoller diagrams ####

# We'll try to plot the seasonality of Dinophysis for each region across time,
# with a Hovmoller diagram

CPR_data_hovmoller <- CPR_data_regions %>%
  # We convert the infinitesimal values into zeros
  mutate(Dinophysis_Total = ifelse(Dinophysis_Total < 1, 0,
                                   Dinophysis_Total)) %>%
  # We then group by Year and Month and region
  group_by(region, Year, Month) %>%
  # Process like before
  # Create a dummy variable that identifies Dinophysis positive counts
  # (0 = no Dinophysis in sample, 1 = 1 or more Dinophysis in sample)
  mutate(Dinodummy = ifelse(Dinophysis_Total > 1, 1, 0)) %>%
  # Summarise the dataset as we wish to
  summarise(
    # mean of Dinophysis counted
    mean_Dino = mean(Dinophysis_Total),
    # Total sum of Dinophysis counted
    sumDino = sum(Dinophysis_Total),
    # Number of samples with Dinophysis
    nsamples_Dino = sum(Dinodummy),
    # Total number of samples
    nsamples = n(),
    .groups = 'keep')

# Create a variable of the proportion of Dinophysis "positive" samples
CPR_data_hovmoller <- CPR_data_hovmoller %>%
  mutate(prop_pos = nsamples_Dino/nsamples)

# Ok. Now we plot.

#--- Plotting Howmoller diagrams ---#

# Let's add color to the background of the facets
# For this, we need a dataframe that links names of regions with specific colors
# The name of the columns in the dataframe is super important here!
df_color <- data.frame(name = 1:4, 
                       color = c('orange1', 'royalblue2', 'red4', 'forestgreen'))
# Ok nice

# One command, 2 versions of the plot (adjust comments accordingly)
ggplot(data = subset(CPR_data_hovmoller, region > 0)) +
  geom_tile(aes(x = Year, y = Month, fill = 
                  # log(prop_pos+1) # Dinophysis presence version
                  log10(nsamples) # Sampling effort version
                )) +
  scale_fill_viridis(
    option = 'plasma' # For sampling effort version
    ## All 3 options below are for the Dinophysis presence version.
    # limits = c(0, log(2)), #
    # breaks = c(log(0+1), log(0.1+1), log(0.25+1), log(0.5+1), log(2)), #
    # labels = c('0', '0.1', '0.25', '0.5', '1') #
  ) +
  facet_wrap_color(facets = 'region', colors = df_color) +
  scale_y_continuous(limits = c(0, 13), 
                     breaks = c(1,2,3,4,5,6,7,8,9,10,11,12)) +
  labs(title = 'B',
       fill = 
                'Sampling effort
         (log10[number of samples])'
#          'Dinophysis presence
# (proportion of samples,
# log scale)'
  ) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.frame = element_rect(color = 'black'),
        legend.ticks = element_line(color = 'black'),
        legend.background = element_rect(fill = 'white', color = 'black'))

# Perfect! We see that the sampling effort in region 1 is too scarce to
# infer anything about Dinophysis phenology in that region. 
## Save that!
# Dinophysis presence version:
# ggsave('Plots/CPR/Fig1B_CPR_Hovmoller_Dinophysis.tiff', height = 95, width = 160, units = 'mm',
#        dpi = 300, compression = 'lzw')
# CPR sampling effort version:
# ggsave('Plots/CPR/FigS4B_CPR_Hovmoller_sampling_effort.tiff', height = 95, width = 160, units = 'mm',
#        dpi = 300, compression = 'lzw')

####-------------------------- End of script ------------------------------####