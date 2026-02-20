#### Script: Map of REPHY sampling sites ###
# Author: V. POCHIC, with help from code written by JY Dias :)
# Last modif: 2026/02/20


## General information ####

## Description: 
# This script makes a map of the REPHY sampling sites used in the study, with
# rivers.

## Files required:
# Data/REPHY_outputs/Season_Dino.csv
# To get the coordonates of the sampling sites.
# Data/Geography/shapefile_rivers/...
# A shapefile for rivers.
## Most stuff about rivers is adapted from a blog post by Milos Popovic
## (https://milospopovic.net/map-rivers-with-sf-and-ggplot2-in-r/)
## Thank you very much Milos, and big thanks to Simon Oiry (https://oirysimon.com/)
## for helping with shapefiles

## Outputs:
# Figure 2 - Map of REPHY sampling sites (only the map, the rest of the figure
# was drawn on powerpoint)

#### Required packages ####

library(tidyverse)
library(mapdata)
library(sf)

####-----------------------------------------------------------------------####
#### Import data ####
# Import data for site coordinates
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


# Import a worldmap
Worldmap <- map_data('worldHires')

## River data
sf_rivers_eu <- read_sf('Data/Geography/shapefile_rivers/HydroRIVERS_v10_eu.shp') %>%
  filter(ORD_FLOW <= 5)

#### Plotting the map ####

# Color palettes
pheno_palette16 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646',
                     '#642C3A', '#DEB1CC', '#FC4D6B', '#791D40')

# Plotting the map (takes a long time with the big shapefile version)
ggplot() + geom_polygon(data = Worldmap, aes(x = long, y = lat, group = group), 
                        fill = "white", color = 'gray10', linewidth = .25)+
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)',
       title = NULL)+
  
  # Rivers
  geom_sf(data = sf_rivers_eu, aes(linewidth = factor(ORD_FLOW), 
                                   alpha = factor(ORD_FLOW)),
          color = 'dodgerblue4') +
  # Alpha and linewidth scales
  scale_alpha_manual(values = c("3" = 1, "4" = .8, "5" = .3), 
                     guide = 'none') + #, "6" = .3, "7" = .2
  scale_linewidth_manual(values = c("3" = .4, "4" = .3, "5" = .2), 
                         guide = 'none') + # , "6" = .1
  # Bound coordinates
  coord_sf(
    crs = NULL,
    xlim = c(-5.5,9.5),
    ylim = c(41,51.5)
  ) +
  
  # points for sampling sites
  geom_point(data = Table_sites, aes(x = Longitude.mean, y = Latitude.mean, 
                              fill = Code_point_Libelle), size =4.5,
             shape = 21, alpha = .75, stroke = .1, color = 'gray10')+
  # fill scale
  scale_fill_discrete(type = pheno_palette16, guide = 'none') +

# Theme
  theme(panel.grid.major = element_line(color = 'gray35', linewidth = .2),
        panel.background = element_rect(fill = '#BBD4F2'),panel.border = element_blank(),
        legend.background = element_rect(fill = "white"),legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 8,color = "black"),
        legend.title = element_text(size = 8, color = "black",face="bold"),
        plot.background = element_rect(color = "white", linewidth = 0.001),
        legend.position = "bottom",plot.title = element_text(size=10,color ="black"),
        axis.text = element_text(size=8),axis.title.x = element_text(size=9),
        axis.title.y = element_text(size=9))

# Save the map
# ggsave('Plots/REPHY/Fig2_map_DinoPhenology_rivers.tiff', dpi = 300,
# width = 80, height = 80,
# units = 'mm', compression = 'lzw')

####-------------------------- End of script ------------------------------####