### CPR data for Dinophysis phenology ###

# V. POCHIC, 2025-11-19

# This script is made for plotting and analysing data from the CPR.

## packages ####

library(tidyverse)
library(viridis)
library(cmocean)
library(ggnewscale)
library(ggmap)
library(mapdata)
library(maps)
library(spatstat.univar)
library(ggpubr)

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

# ggsave('CPR_map_dotplot2.tiff', width = 164, height = 180, units = 'mm',
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
       title = "Dinophysis in the NE Atlantic (CPR data)", 
       fill = 'Dinophysis presence index') +
  theme_bw() +
  theme(legend.position = 'bottom')

# This seems ok

# ggsave('CPR_heatmap_Dino_presence3.tiff', width = 164, height = 180, units = 'mm',
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
       title = "CPR sampling in the NE Atlantic", 
       fill = 'Number of samples (log10)') +
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
       title = "Dinophysis seasonality in the NE Atlantic (CPR data)")+
  theme_bw() +
  theme(legend.position = 'bottom')

# Doesn't seem really decisive

# But maybe we can look at it another way?

ggplot(data = CPR_data_max) +
  geom_point(aes(x = Month, y = lat.bin, color = as_factor(Month)), size = 2.5,
             alpha = .8) +
  scale_color_viridis(discrete = TRUE, guide = 'none') +
  labs(y = 'Latitude (degrees)', x = 'Month with maximum Dinophysis count',
       title = "Month vs latitude (Dino seasonality)")+
  theme_classic()

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
  coord_fixed(xlim=c(-25,
                     max(CPR_data$Longitude)), 
              ylim=c(min(CPR_data$Latitude),
                     max(CPR_data$Latitude)))+
  # Labels
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)',
       title = "Dinophysis seasonality in the NE Atlantic (CPR data)")+
  theme_bw() +
  theme(legend.position = 'bottom')

# Ok

# The dot plot
ggplot(data = CPR_data_season) +
  geom_point(aes(x = mode_month, y = lat.bin, color = as_factor(mode_month)), size = 2.5,
             alpha = .8) +
  scale_color_viridis(discrete = TRUE, guide = 'none') +
  labs(y = 'Latitude (degrees)', x = 'Month with most Dinophysis count',
       title = "Month vs latitude (Dino seasonality)")+
  theme_classic()

# This is almost exactly the same figure than the one with the maximum Dinophysis
# count.

### Dino across latitudes with all the data

CPR_data_Dino <- CPR_data %>%
  # Get the date in the right format
  # Create a date variable in date format
  mutate(Date = paste(Year, Month, Day, sep = '-')) %>%
  mutate(Date = ymd(Date)) %>%
  # Get a julian day variable
  mutate(Julian_day = yday(Date))

# Plot it as a dotplot
ggplot() +
  # Observations without Dinophysis
  geom_point(data = subset(CPR_data_Dino, Dinophysis_Total == 0),
             aes(x = Julian_day, y = Latitude), color = 'gray70', alpha = .95) +
  # Observations with Dinophysis
  geom_point(data = subset(CPR_data_Dino, Dinophysis_Total>0),
             aes(x = Julian_day, y = Latitude, color = log10(Dinophysis_Total)), 
             size = 2.5,
             alpha = .8) +
  scale_color_viridis() +
  labs(y = 'Latitude (degrees)', x = 'Julian day of observation',
       title = "Dino observations (day vs latitude)")+
  theme_classic() +
  theme(legend.position = 'bottom')
  

# Is there a bias in the sampling itself? Let's have a look

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
  # Plotting the tiles of longitude and latitude depending on the most sampled
  # month
  geom_rect(data = CPR_data_sampling,
            aes(xmin = lon.bin-.25, xmax = lon.bin+.25,
                ymin = lat.bin-.25, ymax = lat.bin+.25, # width = .5, height = .5,
                fill = as_factor(mode_month))) +
  scale_fill_viridis(option = 'cividis', discrete = TRUE) +
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

# No real pattern here.

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
  # Get a julian day variable
  mutate(Julian_day = yday(Date)) %>%
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

# We can also plot the regions as polygons on top of the Dinophysis presence map

# For this, we import a tibble containing the necessary information to provide to
# geom_polygon

regions_polygon <- read.csv2('Data/CPR/CPR_regions_polygon.csv', header = TRUE,
                             fileEncoding = 'ISO-8859-1')

# And a georeferenced tibble to plot the id of each region
regions_id <- read.csv2('Data/CPR/CPR_regions_ids.csv', header = TRUE,
                        fileEncoding = 'ISO-8859-1')

ggplot() +
  # Plotting the tiles of longitude and latitude depending on Dinophysis presence
  geom_rect(data = subset(CPR_data_sum, prop_pos<1),
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
  # Plotting the id of each region
  geom_text(data = regions_id, aes(x = lon, y = lat, color = as_factor(region),
                                   label = region), size = 6)+
  scale_color_discrete(guide = 'none') +
  # Limits of the map
  coord_fixed(xlim=c(-25,
                     max(CPR_data$Longitude)), 
              ylim=c(min(CPR_data$Latitude),
                     max(CPR_data$Latitude)))+
  # Labels
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)',
       # title = "Dinophysis in the NE Atlantic (CPR data)", 
       fill = 'Dinophysis presence index')+
  theme_bw() +
  theme(legend.position = 'bottom')

## Save plot
# ggsave('CPR_Dinophysis_presence_regions2.tiff', height = 180, width = 164,
#        units = 'mm', dpi = 300, compression = 'lzw')

### The map for the common plot with REPHY data ####

# Get the data for REPHY sampling sites
# Import data for site coordinates
Season_Dino <- read.csv2('Data/REPHY_outputs/Season_Dino_20250604.csv', header = TRUE, 
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

# Color palette for the REPHY sampling sites
pheno_palette16 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646',
                     '#642C3A', '#DEB1CC', '#FC4D6B', '#791D40')

# COlor palette for the CPR-defined regions
palette_regions4 <- c('orange1', 'royalblue2', 'red4', 'forestgreen')

ggplot() +
  # Plotting the tiles of longitude and latitude depending on Dinophysis presence
  geom_rect(data = subset(CPR_data_sum, nsamples>=3),
            aes(xmin = lon.bin-.25, xmax = lon.bin+.25,
                ymin = lat.bin-.25, ymax = lat.bin+.25, # width = .5, height = .5,
                fill = log10(mean_Dino+1)*prop_pos*log10(nsamples+1))) +
  # adjust the limits of the color scale because of some extremes (prop_pos = 1)
  # that dwarf the rest
  scale_fill_viridis() +
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
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)',
       # title = "Dinophysis in the NE Atlantic (CPR data)", 
       fill = 'Dinophysis presence index')+
  theme_bw() +
  theme(legend.position = 'bottom')

## Save plot
# ggsave('Plots/CPR/CPR_REPHY_Dinophysis_map_composite.tiff', height = 200, width = 164,
#        units = 'mm', dpi = 300, compression = 'lzw')

# Plot the seasonality of Dinophysis depending on the region
ggplot(data = subset(CPR_data_regions, region > 0)) +
  # Here we try to divide by the lowest number of Dinophysis recorded, to get
  # a somewhat continuous distribution of our data.
  geom_point(aes(x = Julian_day, y = (Dinophysis_Total/15000),
                 color = as_factor(region))) +
  facet_wrap(facets = 'region', scales = 'free_y') +
  scale_color_discrete(guide = 'none') +
  theme_classic()

## Save plot
# ggsave('CPR_Dinophysis_seasonality_regions.tiff', height = 164, width = 164,
#        units = 'mm', dpi = 300, compression = 'lzw')

# Looking at timing of maximum
CPR_data_max <- CPR_data_regions %>%
  # Grouping by the bin variables
  group_by(lat.bin, lon.bin) %>%
  # Get only rows with maximum abundance
  top_n(1, Dinophysis_Total) %>%
  # Filter 0 or close to 0 abundances
  filter(Dinophysis_Total > 1)

# Plot the seasonality of Dinophysis depending on the region
ggplot(data = subset(CPR_data_max, region > 0)) +
  geom_point(aes(x = Julian_day, y = Latitude,
                 color = log10(Dinophysis_Total))) +
  geom_smooth(aes(x = Julian_day, y = Latitude),
              method = 'lm', alpha = .3) +
  facet_wrap(facets = 'region', scales= 'free_y') +
  scale_color_viridis() +
  theme_classic() +
  theme(legend.position = 'bottom')

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

### First, let's create a Hovmoller diagram of sampling effort in the different 
# regions

ggplot(data = subset(CPR_data_hovmoller, region > 0)) +
  geom_tile(aes(x = Year, y = Month, fill = log10(nsamples+1))) +
  scale_fill_viridis(option = 'plasma') +
  facet_wrap(facets = 'region') +
  scale_y_continuous(limits = c(0, 13), 
                     breaks = c(1,2,3,4,5,6,7,8,9,10,11,12)) +
  labs(title = 'Sampling effort over the study period,
by region',
       fill = 'Number of samples (log10)') +
  theme_classic() +
  theme(legend.position = 'bottom')

# Save that!
# ggsave('Plots/CPR/CPR_Hovmoller_sampling.tiff', height = 150, width = 164, units = 'mm',
#        dpi = 300, compression = 'lzw')

### Then, Hovmoller diagram of Dinophysis presence

ggplot(data = subset(CPR_data_hovmoller, region > 0)) +
  geom_tile(aes(x = Year, y = Month, fill = log(prop_pos+1))) +
  scale_fill_viridis(breaks = c(log(0+1), log(0.1+1), log(0.25+1), (log(0.5+1))),
                     labels = c('0', '0.1', '0.25', '0.5')) +
  facet_wrap(facets = 'region') +
  scale_y_continuous(limits = c(0, 13), 
                     breaks = c(1,2,3,4,5,6,7,8,9,10,11,12)) +
  labs(title = 'Dinophysis seasonality by region, CPR data',
       fill = 'Dinophysis presence 
(proportion of samples,
log scale)') +
  theme_classic() +
  theme(legend.position = 'bottom')

# Save that!
# ggsave('Plots/CPR/CPR_Hovmoller_proportion_4regions.tiff', height = 125, width = 164, units = 'mm',
#        dpi = 300, compression = 'lzw')

### Then, Hovmoller diagram of Dinophysis abundance

ggplot(data = subset(CPR_data_hovmoller, region > 0)) +
  geom_tile(aes(x = Year, y = Month, 
                fill = log10(mean_Dino+1))) +
  scale_fill_viridis() +
  facet_wrap(facets = 'region') +
  scale_y_continuous(limits = c(0, 13), 
                     breaks = c(1,2,3,4,5,6,7,8,9,10,11,12)) +
  labs(title = 'Dinophysis seasonality over the study period, 
by region',
       fill = 'Mean Dinophysis abundance (log10)') +
  theme_classic() +
  theme(legend.position = 'bottom')

# Save that!
# ggsave('Plots/CPR/CPR_Hovmoller_abundance.tiff', height = 150, width = 164, units = 'mm',
#        dpi = 300, compression = 'lzw')

### Plots for publication ####
# Let's combine 2 plots: the map and the hovmoller of Dino presence.
# ggarrange won't do the trick here, so we'll build 2 compatible plots that
# we will then arrange on ppt

# First, the map
ggplot() +
  # Plotting the tiles of longitude and latitude depending on Dinophysis presence
  geom_rect(data = subset(CPR_data_sum, nsamples>=3),
            aes(xmin = lon.bin-.25, xmax = lon.bin+.25,
                ymin = lat.bin-.25, ymax = lat.bin+.25, # width = .5, height = .5,
                fill = log10(nsamples))) + #log(prop_pos+1)
  # adjust the limits of the color scale so it is compatible with the second
  # plot. No legend here.
  scale_fill_viridis(option = 'plasma'
    # limits = c(0, 0.694),
    # breaks = c(log(0+1), log(0.1+1), log(0.25+1), (log(0.5+1))),
    # labels = c('0', '0.1', '0.25', '0.5'), guide = 'none'
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
       # title = "Dinophysis in the NE Atlantic (CPR data)"
       )+
  theme_bw() +
  theme(legend.position = 'bottom')

# Save the map
# ggsave('Plots/CPR/CPR_sampling_effort_map_composite_pub.tiff', height = 200, width = 164,
#        units = 'mm', dpi = 300, compression = 'lzw')

# Then, the Hovmoller diagram
ggplot(data = subset(CPR_data_hovmoller, region > 0)) +
  geom_tile(aes(x = Year, y = Month, fill = log10(nsamples))) + # log(prop_pos+1)
  scale_fill_viridis(option = 'plasma'
    # limits = c(0,0.694),
    # breaks = c(log(0+1), log(0.1+1), log(0.25+1), log(0.5+1), log(2)),
    #                  labels = c('0', '0.1', '0.25', '0.5', '1')
    ) +
  facet_wrap(facets = 'region') +
  scale_y_continuous(limits = c(0, 13), 
                     breaks = c(1,2,3,4,5,6,7,8,9,10,11,12)) +
  labs(title = 'B',
       fill = 'Sampling effort 
(number of samples, log10)'
#        'Dinophysis presence 
# (proportion of samples,
# log scale)'
       ) +
  theme_classic() +
  theme(legend.position = 'bottom')

# Save that!
# ggsave('Plots/CPR/CPR_Hovmoller_sampling_effort_4regions_pub.tiff', height = 125, width = 164, units = 'mm',
#        dpi = 300, compression = 'lzw')

### Hovmoller diagrams by fortnight (does not work well) ####

# Is it better to visualise by fortnight rather than month?
  
# Let's see

CPR_data_hovmoller_fortnight <- CPR_data_regions %>%
  # We convert the infinitesimal values into zeros
  mutate(Dinophysis_Total = ifelse(Dinophysis_Total < 1, 0,
                                   Dinophysis_Total)) %>%
  # We then group by Year and Fortnight and region
  group_by(region, Year, Fortnight) %>%
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

### First, let's create a Hovmoller diagram of sampling effort in the different 
# regions

ggplot(data = subset(CPR_data_hovmoller_fortnight, region > 0)) +
  geom_tile(aes(x = Year, y = Fortnight, fill = log10(nsamples+1))) +
  scale_fill_viridis(option = 'plasma') +
  facet_wrap(facets = 'region') +
  scale_y_continuous(limits = c(0, 27), 
                     breaks = c(1,5,10,15,26)) +
  labs(title = 'Sampling effort over the study period,
by region',
       fill = 'Number of samples (log10)') +
  theme_classic() +
  theme(legend.position = 'bottom')

# No, it does not improve visualisation. Many fortnights were left unsampled 
# over the study period

### Then, Hovmoller diagram of mean Dinophysis abundance

ggplot(data = subset(CPR_data_hovmoller_fortnight, region > 0)) +
  geom_tile(aes(x = Year, y = Fortnight, fill = log10(mean_Dino+1))) +
  scale_fill_viridis() +
  facet_wrap(facets = 'region') +
  scale_y_continuous(limits = c(0, 27), 
                     breaks = c(1,5,10,15,26)) +
  labs(title = 'Dinophysis seasonality over the study period, 
by region',
       fill = 'Mean abundance of Dinophysis (log10)') +
  theme_classic() +
  theme(legend.position = 'bottom')

# Same thing
