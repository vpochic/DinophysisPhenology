### CPR data for Dinophysis phenology ###

# V. POCHIC, 2025-04-15

# This script is made for plotting and analysing data from the CPR.

## packages ####

library(tidyverse)
library(viridis)
library(ggmap)
library(mapdata)
library(maps)
library(spatstat.univar)

## import data ####

# CPR data
CPR_data <- read.csv('Data/CPR/CPR_DREAMproject_Data_Dinophysis_03022025.csv',
                      header = TRUE, fileEncoding = 'ISO-8859-1')

# A world map
Worldmap <- map_data('worldHires')

## Plot the data on a map ####

# plotting the map

ggplot() +
  # Plotting the land
  geom_polygon(data = Worldmap, aes(x = long, y = lat, group = group), 
                        fill = "white", color = 'gray10', linewidth = .25)+
  coord_fixed(xlim=c(-25,
                     max(CPR_data$Longitude)), 
              ylim=c(min(CPR_data$Latitude),
                     max(CPR_data$Latitude)), ratio=1)+
  # Plotting the samples without Dinophysis
  geom_point(data = subset(CPR_data, Dinophysis_Total == 0), 
             aes(x = Longitude, y = Latitude),
             color = 'black', size = .5, alpha = .8) +
  # Plotting the samples with Dinophysis
  geom_point(data = subset(CPR_data, Dinophysis_Total > 0), 
             aes(x = Longitude, y = Latitude),
             color = 'red', size = .5, alpha = .8) +
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)',
       title = "CPR data in the NE Atlantic")+
  theme_classic()

# Nice one!

# ggsave('CPR_map_dotplot.tiff', width = 164, height = 180, units = 'mm',
#        dpi = 300, compression = 'lzw')

## Bin the data by coordinate tiles ####

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
  mutate(Dinodummy = ifelse(Dinophysis_Total > 0, 1, 0)) %>%
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



## Plotting a heatmap of Dinophysis presence ####

# plotting the map

ggplot() +
  # Plotting the tiles of longitude and latitude depending on Dinophysis 
  # presence
  geom_rect(data = CPR_data_sum,
             aes(xmin = lon.bin-.25, xmax = lon.bin+.25,
                 ymin = lat.bin-.25, ymax = lat.bin+.25, # width = .5, height = .5,
                 fill = prop_pos*log10(nsamples+1))) +
  scale_fill_viridis() +
  # Plotting the land
  geom_polygon(data = Worldmap, aes(x = long, y = lat, group = group), 
               fill = "grey80", color = 'gray10', linewidth = .25)+
  coord_fixed(xlim=c(-25,
                     max(CPR_data$Longitude)), 
              ylim=c(min(CPR_data$Latitude),
                     max(CPR_data$Latitude)))+
  # Labels
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)',
       title = "Dinophysis in the NE Atlantic (CPR data)")+
  theme_bw() +
  theme(legend.position = 'bottom')

# This seems ok

# ggsave('CPR_heatmap_Dino_presence.tiff', width = 164, height = 180, units = 'mm',
#        dpi = 300, compression = 'lzw')

## Plotting a heatmap of sampling effort ####

# plotting the map

ggplot() +
  # Plotting the tiles of longitude and latitude depending on Dinophysis 
  # presence
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
       title = "CPR sampling in the NE Atlantic")+
  theme_bw() +
  theme(legend.position = 'bottom')

# This seems ok

# ggsave('CPR_sampling_effort.tiff', width = 164, height = 180, units = 'mm',
#        dpi = 300, compression = 'lzw')
