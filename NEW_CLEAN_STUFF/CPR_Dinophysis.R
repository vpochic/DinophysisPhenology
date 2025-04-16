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
                 fill = log10(mean_Dino+1)*prop_pos*log10(nsamples+1))) +
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

# ggsave('CPR_heatmap_Dino_presence2.tiff', width = 164, height = 180, units = 'mm',
#        dpi = 300, compression = 'lzw')

## Plotting a heatmap of sampling effort ####

# plotting the map

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
       title = "CPR sampling in the NE Atlantic")+
  theme_bw() +
  theme(legend.position = 'bottom')

# This seems ok

# ggsave('CPR_sampling_effort.tiff', width = 164, height = 180, units = 'mm',
#        dpi = 300, compression = 'lzw')

## Exploring the seasonality ####

# First, let's just plot the month of maximum Dinophysis abundance recorded

# Now, we can summarise the dataset by grouping it thanks to our new bins
CPR_data_max <- CPR_data_binned %>%
  # Grouping by the bin variables
  group_by(lat.bin, lon.bin) %>%
  # Get only rows with maximum abundance
  top_n(1, Dinophysis_Total) %>%
  # Filter 0 abundances
  filter(Dinophysis_Total > 0) %>%
  # Create a date variable in date format
  mutate(Date = paste(Year, Month, Day, sep = '-')) %>%
  mutate(Date = ymd(Date)) %>%
  # Get a julian day variable
  mutate(Julian_day = yday(Date))

## Plot the month of the recorded maximum for each bin where Dinophysis was
# observed

ggplot() +
  # Plotting the tiles of longitude and latitude depending on Dinophysis 
  # presence
  geom_rect(data = CPR_data_max,
            aes(xmin = lon.bin-.25, xmax = lon.bin+.25,
                ymin = lat.bin-.25, ymax = lat.bin+.25, # width = .5, height = .5,
                fill = Month)) +
  scale_fill_viridis(option = 'viridis') +
  # Plotting the land
  geom_polygon(data = Worldmap, aes(x = long, y = lat, group = group), 
               fill = "grey80", color = 'gray10', linewidth = .25)+
  coord_fixed(xlim=c(min(CPR_data$Longitude),
                     max(CPR_data$Longitude)), 
              ylim=c(min(CPR_data$Latitude),
                     max(CPR_data$Latitude)))+
  # Labels
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)',
       title = "Dinophysis seasonality the NE Atlantic (CPR data)")+
  theme_bw() +
  theme(legend.position = 'bottom')

# Doesn't seem really decisive

# Let's try something else

# We'll try to get the mode month (the month with the most observations) for
# each bin

# For this we define a function to calculate the mode
find_mode <- function(x) {
  u <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}

CPR_data_season <- CPR_data_binned %>%
  # Grouping by the bin variables
  group_by(lat.bin, lon.bin) %>%
  # Filter only positive Dinophysis observations
  filter(Dinophysis_Total > 0) %>%
  # Get the mode month
  reframe(mode_month = find_mode(Month),
          # We also get the sum of Dinophysis observed to decide which to keep
          # when there is more than 1 mode
          sum_Dino = sum(Dinophysis_Total)) %>%
  # Now we keep only the ones with the max abundance of Dino
  group_by(lat.bin, lon.bin) %>%
  top_n(1, sum_Dino)

## Let's look at the map

ggplot() +
  # Plotting the tiles of longitude and latitude depending on Dinophysis 
  # presence
  geom_rect(data = CPR_data_season,
            aes(xmin = lon.bin-.25, xmax = lon.bin+.25,
                ymin = lat.bin-.25, ymax = lat.bin+.25, # width = .5, height = .5,
                fill = mode_month)) +
  scale_fill_viridis(option = 'viridis') +
  # Plotting the land
  geom_polygon(data = Worldmap, aes(x = long, y = lat, group = group), 
               fill = "grey80", color = 'gray10', linewidth = .25)+
  coord_fixed(xlim=c(min(CPR_data$Longitude),
                     max(CPR_data$Longitude)), 
              ylim=c(min(CPR_data$Latitude),
                     max(CPR_data$Latitude)))+
  # Labels
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)',
       title = "Dinophysis seasonality in the NE Atlantic (CPR data)")+
  theme_bw() +
  theme(legend.position = 'bottom')

# But there may be a bias in the monthly sampling. Let's have a look

# Overall dataset
ggplot(data = CPR_data_binned) + 
  geom_bar(aes(x = as_factor(Month))) +
  scale_x_discrete() +
  theme_classic()

# Seems approximately even. Is this the case across the map?

CPR_data_sampling <- CPR_data_binned %>%
  # Grouping by the bin variables
  group_by(lat.bin, lon.bin) %>%
  # Get the mode month
  reframe(mode_month = find_mode(Month),
          # We also get the sum of Dinophysis observed to decide which to keep
          # when there is more than 1 mode (arbitrary)
          sum_Dino = sum(Dinophysis_Total)) %>%
  # Now we keep only the ones with the max abundance of Dino
  group_by(lat.bin, lon.bin) %>%
  top_n(1, sum_Dino)


ggplot() +
  # Plotting the tiles of longitude and latitude depending on Dinophysis 
  # presence
  geom_rect(data = CPR_data_sampling,
            aes(xmin = lon.bin-.25, xmax = lon.bin+.25,
                ymin = lat.bin-.25, ymax = lat.bin+.25, # width = .5, height = .5,
                fill = as_factor(mode_month))) +
  scale_fill_viridis(option = 'viridis', discrete = TRUE) +
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

### Regionalising ####

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
      lat.bin <= 53 & lat.bin >= 45 & lon.bin <= -1 & lon.bin >= -10.5,
      2,
      ifelse(
        # Western North Sea (#3)
        lat.bin <= 60 & lat.bin >= 52.5 & lon.bin < 3.5 & lon.bin >= -3.5,
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
  # Get a julian day variable
  mutate(Julian_day = yday(Date))

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

# We can also plot the regions as polygons on top of the Dinophysis presence map

# For this, we import a tibble containing the necessary information to provide to
# geom_polygon

regions_polygon <- read.csv2('Data/CPR/CPR_regions_polygon.csv', header = TRUE,
                             fileEncoding = 'ISO-8859-1')

ggplot() +
  # Plotting the tiles of longitude and latitude depending on Dinophysis presence
  geom_rect(data = CPR_data_sum,
            aes(xmin = lon.bin-.25, xmax = lon.bin+.25,
                ymin = lat.bin-.25, ymax = lat.bin+.25, # width = .5, height = .5,
                fill = log10(mean_Dino+1)*prop_pos*log10(nsamples+1))) +
  scale_fill_viridis() +
  # Plotting the land
  geom_polygon(data = Worldmap, aes(x = long, y = lat, group = group), 
               fill = "grey80", color = 'gray10', linewidth = .25)+
  # Plotting the limits of the regions
  geom_polygon(data = regions_polygon, aes(x = lon, y = lat, group = region,
                                           color = as_factor(region)), 
               fill = "transparent", linewidth = 1)+
  scale_color_discrete() +
  # Limits of the map
  coord_fixed(xlim=c(-25,
                     max(CPR_data$Longitude)), 
              ylim=c(min(CPR_data$Latitude),
                     max(CPR_data$Latitude)))+
  # Labels
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)',
       title = "Dinophysis in the NE Atlantic (CPR data)")+
  theme_bw() +
  theme(legend.position = 'bottom')

# Plot the seasonality of Dinophysis depending on the region
ggplot(data = CPR_data_regions) +
  geom_point(aes(x = Julian_day, y = Dinophysis_Total,
                 color = as_factor(region))) +
  facet_wrap(facets = 'region', scales = 'free_y') +
  scale_color_discrete(guide = 'none') +
  theme_classic()
