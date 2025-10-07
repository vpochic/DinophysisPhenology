###################### A map of REPHY sites for Dinophysis phenology ###
### V. POCHIC, 2025/09/02
##(Thank you J-YDi :)
## Some stuff about rivers is heavily adapted from a blog post by Milos Popovic
## (https://milospopovic.net/map-rivers-with-sf-and-ggplot2-in-r/)
## Thank you very much Milos, and big thanks to Simon Oiry (https://oirysimon.com/)
# for helping with shapefiles

#### Required packages ####

library(tidyverse)
library(ggmap)
library(mapdata)
library(maps)
library(ggthemes)
library(raster)
library(rasterVis)
library(sf)
library(httr)

#### Import data ####
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

# Exclude Mediterranean sites
Table_sites_Atl <- Table_sites %>%
  filter(Code_point_Libelle %in% c('Point 1 Boulogne', 'At so',
                                   'Antifer ponton pétrolier', 'Cabourg',
                                   'les Hébihens', 'Loguivy',
                                   'Men er Roue', 'Ouest Loscolo',
                                   'Le Cornard', 'Auger',
                                   'Arcachon - Bouée 7', 'Teychan bis'))

# Only Seine bay and Southern Brittany
Table_sites_Bio <- Table_sites %>%
  filter(Code_point_Libelle %in% c('Antifer ponton pétrolier', 'Cabourg',
                                   'Men er Roue', 'Ouest Loscolo'))


# Import a worldmap
Worldmap <- map_data('worldHires')

## River data
sf_rivers <- read_sf('Data/Geography/shapefile_rivers/Riviere_inf_4.shp') %>%
  mutate(Width = ifelse(ORD_FLOW == 3, 1, .75))
sf_rivers <- st_cast(sf_rivers, "MULTILINESTRING")

#### Make the map ####

# Color palettes
pheno_palette16 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646',
                     '#642C3A', '#DEB1CC', '#FC4D6B', '#791D40')

pheno_palette12 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646')

pheno_palette4 <- c('red3', 'orangered', 
                     '#2156A1', '#5995E3')

# plotting the map
ggplot() + geom_polygon(data = Worldmap, aes(x = long, y = lat, group = group), 
                        fill = "white", color = 'gray10', size = .25)+
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)',
       title = NULL)+
  
  # Rivers
  geom_sf(data = sf_rivers, color = 'dodgerblue4',
          linewidth = .4) +
  # Color scale for rivers
  scale_color_discrete(type = palette_rivers) +
  # Bound coordinates
  coord_sf(
    crs = NULL,
    xlim = c(-5.5,9.5),
    ylim = c(41,51.5)
  ) +
  
  geom_point(data = Table_sites, aes(x = Longitude.mean, y = Latitude.mean, 
                              fill = Code_point_Libelle), size =4.5,
             shape = 21, alpha = .75, stroke = .1, color = 'gray10')+
  
  
  
  # Pas de Calais
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Point 1 Boulogne"),
            aes(x = Longitude.mean + 1.7, y = Latitude.mean,label = "Point 1 Boulogne" ), stat = 'unique',
            size = 2,color="black",fontface = "bold") +
  geom_text(data = filter(Table_sites,Code_point_Libelle == "At so"),
            aes(x = Longitude.mean + .7, y = Latitude.mean,label = "At so" ), stat = 'unique',
            size = 2,color="black",fontface = "bold") +
  # Baie de Seine
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Antifer ponton pétrolier"),
            aes(x = Longitude.mean + 2, y = Latitude.mean,label = "Antifer ponton pétrolier" ), stat = 'unique',
            size = 2,color="black",fontface = "bold") +
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Cabourg"),
            aes(x = Longitude.mean + 1, y = Latitude.mean,label = "Cabourg" ), stat = 'unique', 
            size = 2,color="black",fontface = "bold") +
  # Bretagne Nord
  geom_text(data = filter(Table_sites,Code_point_Libelle == "les Hébihens"),
            aes(x = Longitude.mean + 1.5, y = Latitude.mean,label = "les Hébihens" ), stat = 'unique',
            size = 2,color="black",fontface = "bold") +
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Loguivy"),
            aes(x = Longitude.mean - .9, y = Latitude.mean +0.2,label = "Loguivy" ),
            stat = 'unique', size = 2,color="black",fontface = "bold") +
  # Bretagne Sud
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Men er Roue"),
            aes(x = Longitude.mean - 1.3, y = Latitude.mean,label = "Men er Roue" ), stat = 'unique',
            size = 2,color="black",fontface = "bold") +
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Ouest Loscolo"),
            aes(x = Longitude.mean + 1.3, y = Latitude.mean,label = "Ouest Loscolo" ), stat = 'unique',
            size = 2,color="black",fontface = "bold") +
  # Pertuis
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Le Cornard"),
            aes(x = Longitude.mean + 1.1, y = Latitude.mean,label = "Le Cornard" ), stat = 'unique',
            size = 2,color="black",fontface = "bold") +
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Auger"),
            aes(x = Longitude.mean + .8, y = Latitude.mean,label = "Auger" ), stat = 'unique',
            size = 2,color="black",fontface = "bold") +
  # Arcachon
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Teychan bis"),
            aes(x = Longitude.mean + 1.17, y = Latitude.mean,label = "Teychan bis" ), stat = 'unique',
            size = 2,color="black",fontface = "bold") +
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Arcachon - Bouée 7"),
            aes(x = Longitude.mean + 2, y = Latitude.mean,label = "Arcachon - Bouée 7" ),
            stat = 'unique', size = 2,color="black",fontface = "bold") +
  # Méditerranée
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Parc Leucate 2")[1,],
            aes(x = Longitude.mean - 1.5, y = Latitude.mean -0.11,label = "Parc Leucate 2" ), stat = 'unique',
            size = 2,color="black",fontface = "bold") +
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Sète mer"),
            aes(x = Longitude.mean - 1.1, y = Latitude.mean -0.05,label = "Sète mer" ), stat = 'unique',
            size = 2,color="black",fontface = "bold") +
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Bouzigues (a)")[1,],
            aes(x = Longitude.mean - 1.5, y = Latitude.mean +0.15,label = "Bouzigues (a)" ), stat = 'unique',
            size = 2,color="black",fontface = "bold") +
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Diana centre")[1,],
            aes(x = Longitude.mean - 1.4, y = Latitude.mean,label = "Diana centre" ), stat = 'unique',
            size = 2,color="black",fontface = "bold") +
  
  # geom_text(aes(x = 2.75, y = 48,label = "FRANCE" ), stat = 'unique', 
  #           size = 7,color="burlywood4",fontface = "bold") +
  
  geom_text(aes(x = -3.5, y = 46, label = 'Atlantic Ocean'), stat = 'unique', 
            size = 3,color = "lightblue4",fontface="italic",angle=-45) +
  geom_text(aes(x = -3.2, y = 50, label = 'English Channel'), stat = 'unique', 
            size = 3,color = "lightblue4",fontface="italic",angle=28) +
  geom_text(aes(x = 5, y = 42, label = 'Mediterranean Sea'), stat = 'unique', 
            size = 3,color = "lightblue4",fontface="italic",angle=8) +
  
  scale_fill_discrete(type = pheno_palette16, guide = 'none') +
  theme(panel.grid.major = element_line(color = 'gray10', size = .25),
        panel.background = element_rect(fill = '#BBD4F2'),panel.border = element_blank(),
        legend.background = element_rect(fill = "white"),legend.key = element_rect(fill = "white"),
        legend.text = element_text(size=8,color = "black"),
        legend.title = element_text(size = 8, color = "black",face="bold"),
        plot.background = element_rect(color = "white", linewidth = 0.001),
        legend.position = "bottom",plot.title = element_text(size=10,color ="black"),
        axis.text = element_text(size=10),axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11))

# Save the map
# ggsave('Plots/REPHY/MAP_DinoPhenology_rivers.tiff', dpi = 300, width = 100, height = 100,
# units = 'mm', compression = 'lzw')

# plotting the map (without Mediterranean sites)
ggplot() + geom_polygon(data = Worldmap, aes(x = long, y = lat, group = group), 
                        fill = "white", color = 'gray10', size = .25)+
  coord_fixed(xlim=c(-5.5,9.5), ylim=c(41,51.5), ratio=1.4)+
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)',
       title = "Maps of the REPHY sampling stations")+
  theme_gdocs()+
  geom_point(data = Table_sites_Atl, aes(x = Longitude.mean, y = Latitude.mean, 
                                     fill = Code_point_Libelle), size =8,
             shape = 21, alpha = .75, stroke = .1, color = 'gray10')+
  # Pas de Calais
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Point 1 Boulogne"),
            aes(x = Longitude.mean + 1.7, y = Latitude.mean,label = "Point 1 Boulogne" ), stat = 'unique',
            size = 3,color="black",fontface = "bold") +
  geom_text(data = filter(Table_sites,Code_point_Libelle == "At so"),
            aes(x = Longitude.mean + .7, y = Latitude.mean,label = "At so" ), stat = 'unique',
            size = 3,color="black",fontface = "bold") +
  # Baie de Seine
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Antifer ponton pétrolier"),
            aes(x = Longitude.mean + 2, y = Latitude.mean,label = "Antifer ponton pétrolier" ), stat = 'unique',
            size = 3,color="black",fontface = "bold") +
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Cabourg"),
            aes(x = Longitude.mean + 1, y = Latitude.mean,label = "Cabourg" ), stat = 'unique', 
            size = 3,color="black",fontface = "bold") +
  # Bretagne Nord
  geom_text(data = filter(Table_sites,Code_point_Libelle == "les Hébihens"),
            aes(x = Longitude.mean + 1.5, y = Latitude.mean,label = "les Hébihens" ), stat = 'unique',
            size = 3,color="black",fontface = "bold") +
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Loguivy"),
            aes(x = Longitude.mean - .9, y = Latitude.mean +0.2,label = "Loguivy" ),
            stat = 'unique', size = 3,color="black",fontface = "bold") +
  # Bretagne Sud
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Men er Roue"),
            aes(x = Longitude.mean - 1.3, y = Latitude.mean,label = "Men er Roue" ), stat = 'unique',
            size = 3,color="black",fontface = "bold") +
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Ouest Loscolo"),
            aes(x = Longitude.mean + 1.3, y = Latitude.mean,label = "Ouest Loscolo" ), stat = 'unique',
            size = 3,color="black",fontface = "bold") +
  # Pertuis
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Le Cornard"),
            aes(x = Longitude.mean + 1.1, y = Latitude.mean,label = "Le Cornard" ), stat = 'unique',
            size = 3,color="black",fontface = "bold") +
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Auger"),
            aes(x = Longitude.mean + .8, y = Latitude.mean,label = "Auger" ), stat = 'unique',
            size = 3,color="black",fontface = "bold") +
  # Arcachon
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Teychan bis"),
            aes(x = Longitude.mean + 1.17, y = Latitude.mean,label = "Teychan bis" ), stat = 'unique',
            size = 3,color="black",fontface = "bold") +
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Arcachon - Bouée 7"),
            aes(x = Longitude.mean + 2, y = Latitude.mean,label = "Arcachon - Bouée 7" ),
            stat = 'unique', size = 3,color="black",fontface = "bold") +
  
  geom_text(aes(x = 2.75, y = 48,label = "FRANCE" ), stat = 'unique', 
            size = 7,color="burlywood4",fontface = "bold") +
  
  geom_text(aes(x = -3.5, y = 46, label = 'Atlantic Ocean'), stat = 'unique', 
            size = 4,color = "lightblue4",fontface="italic",angle=-45) +
  geom_text(aes(x = -3.2, y = 50, label = 'English Channel'), stat = 'unique', 
            size = 4,color = "lightblue4",fontface="italic",angle=28) +
  geom_text(aes(x = 5, y = 42, label = 'Mediterranean Sea'), stat = 'unique', 
            size = 4,color = "lightblue4",fontface="italic",angle=8) +
  
  scale_fill_discrete(type = pheno_palette16, guide = 'none') +
  theme(panel.grid.major = element_line(color = 'gray10', size = .25),
        panel.background = element_rect(fill = '#BBD4F2'),panel.border = element_blank(),
        legend.background = element_rect(fill = "white"),legend.key = element_rect(fill = "white"),
        legend.text = element_text(size=8,color = "black"),
        legend.title = element_text(size = 8, color = "black",face="bold"),
        plot.background = element_rect(color = "white", linewidth = 0.001),
        legend.position = "bottom",plot.title = element_text(size=10,color ="black"),
        axis.text = element_text(size=8),axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))

# Save the map
# ggsave('MAP_Atl_DinoPhenology.tiff', dpi = 600, width = 200, height = 200,
# units = 'mm', compression = 'lzw')

# plotting the map (Only Seine Bay and Southern Brittany)
ggplot() + geom_polygon(data = Worldmap, aes(x = long, y = lat, group = group), 
                        fill = "white", color = 'gray10', size = .25)+
  coord_fixed(xlim=c(-5.5,9.5), ylim=c(41,51.5), ratio=1.4)+
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)',
       title = "Maps of the REPHY sampling stations")+
  theme_gdocs()+
  geom_point(data = Table_sites_Bio, aes(x = Longitude.mean, y = Latitude.mean, 
                                         fill = Code_point_Libelle), size =8,
             shape = 21, alpha = .75, stroke = .1, color = 'gray10')+
  # Baie de Seine
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Antifer ponton pétrolier"),
            aes(x = Longitude.mean + 2, y = Latitude.mean,label = "Antifer ponton pétrolier" ), stat = 'unique',
            size = 3,color="black",fontface = "bold") +
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Cabourg"),
            aes(x = Longitude.mean + 1, y = Latitude.mean,label = "Cabourg" ), stat = 'unique', 
            size = 3,color="black",fontface = "bold") +
  # Bretagne Sud
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Men er Roue"),
            aes(x = Longitude.mean - 1.3, y = Latitude.mean,label = "Men er Roue" ), stat = 'unique',
            size = 3,color="black",fontface = "bold") +
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Ouest Loscolo"),
            aes(x = Longitude.mean + 1.3, y = Latitude.mean,label = "Ouest Loscolo" ), stat = 'unique',
            size = 3,color="black",fontface = "bold") +
  
  geom_text(aes(x = 2.75, y = 48,label = "FRANCE" ), stat = 'unique', 
            size = 7,color="burlywood4",fontface = "bold") +
  
  geom_text(aes(x = -3.5, y = 46, label = 'Atlantic Ocean'), stat = 'unique', 
            size = 4,color = "lightblue4",fontface="italic",angle=-45) +
  geom_text(aes(x = -3.2, y = 50, label = 'English Channel'), stat = 'unique', 
            size = 4,color = "lightblue4",fontface="italic",angle=28) +
  geom_text(aes(x = 5, y = 42, label = 'Mediterranean Sea'), stat = 'unique', 
            size = 4,color = "lightblue4",fontface="italic",angle=8) +
  
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  theme(panel.grid.major = element_line(color = 'gray10', size = .25),
        panel.background = element_rect(fill = '#BBD4F2'),panel.border = element_blank(),
        legend.background = element_rect(fill = "white"),legend.key = element_rect(fill = "white"),
        legend.text = element_text(size=8,color = "black"),
        legend.title = element_text(size = 8, color = "black",face="bold"),
        plot.background = element_rect(color = "white", linewidth = 0.001),
        legend.position = "bottom",plot.title = element_text(size=10,color ="black"),
        axis.text = element_text(size=8),axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))

# Save the map
# ggsave('MAP_Bio_DinoPhenology.tiff', dpi = 600, width = 200, height = 200,
# units = 'mm', compression = 'lzw')

# Map of Southern Brittany area
ggplot() + geom_polygon(data = Worldmap, aes(x = long, y = lat, group = group), 
                        fill = "#869B7B", color = 'gray10', size = .25)+
  coord_fixed(xlim=c(-3.5,-1.5), ylim=c(47,48), ratio=1.4)+
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)')+
  # Nantes
  geom_point(aes(x = -1.55, y = 47.22), 
             fill = '#0A1635', size =8,
             shape = 21, alpha = .75, stroke = .1, color = 'gray10')+
  # Lorient
  geom_point(aes(x = -3.36, y = 47.75), 
             fill = '#0A1635', size =8,
             shape = 21, alpha = .75, stroke = .1, color = 'gray10')+
  # Vannes
  geom_point(aes(x = -2.76, y = 47.66),
             fill = '#0A1635', size =8,
             shape = 21, alpha = .75, stroke = .1, color = 'gray10')+
  # Lorient
  geom_text(aes(x = -3.36, y = 47.83, label = 'Lorient'), stat = 'unique', 
            size = 9,color = "gray10") +
  # Vannes
  geom_text(aes(x = -2.76, y = 47.73, label = 'Vannes'), stat = 'unique', 
            size = 9,color = "gray10") +
  # Nantes
  geom_text(aes(x = -1.575, y = 47.3, label = 'Nantes'), stat = 'unique', 
            size = 9,color = "gray10") +
  
  geom_text(aes(x = -3.025, y = 47.125, label = 'Atlantic Ocean'), stat = 'unique', 
            size = 10,color = "lightblue4", fontface="italic") +
  
  
  theme(panel.grid.major = element_line(color = 'gray10', size = .25),
        panel.background = element_rect(fill = '#BBD4F2'),panel.border = element_blank(),
        legend.background = element_rect(fill = "white"),legend.key = element_rect(fill = "white"),
        legend.text = element_text(size=8,color = "black"),
        legend.title = element_text(size = 8, color = "black",face="bold"),
        plot.background = element_rect(color = "white", linewidth = 0.001),
        legend.position = "bottom",plot.title = element_text(size=10,color ="black"),
        axis.text = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))

# Save the map
# ggsave('MAP_SBrittany_SFE2.tiff', dpi = 600, width = 230, height = 200,
# units = 'mm', compression = 'lzw')

### Only Cabourg and Ouest Loscolo (for COAST video)
pheno_palette2 <- c('orangered', '#5995E3')

# plotting the map (Only Seine Bay and Southern Brittany)
ggplot() + geom_polygon(data = Worldmap, aes(x = long, y = lat, group = group), 
                        fill = "white", color = 'gray10', size = .25)+
  coord_fixed(xlim=c(-5.5,9.5), ylim=c(41,51.5), ratio=1.4)+
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)',
       title = "Maps of the REPHY sampling stations")+
  theme_gdocs()+
  geom_point(data = subset(Table_sites_Bio, Code_point_Libelle %in% c('Cabourg', 'Ouest Loscolo')), 
             aes(x = Longitude.mean, y = Latitude.mean, 
                                         fill = Code_point_Libelle), size =8,
             shape = 21, alpha = .75, stroke = .1, color = 'gray10')+
  # Baie de Seine
  # geom_text(data = filter(Table_sites,Code_point_Libelle == "Antifer ponton pétrolier"),
  #           aes(x = Longitude.mean + 2, y = Latitude.mean,label = "Antifer ponton pétrolier" ), stat = 'unique',
  #           size = 3,color="black",fontface = "bold") +
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Cabourg"),
            aes(x = Longitude.mean + 1, y = Latitude.mean,label = "Cabourg" ), stat = 'unique', 
            size = 3,color="black",fontface = "bold") +
  # Bretagne Sud
  # geom_text(data = filter(Table_sites,Code_point_Libelle == "Men er Roue"),
  #           aes(x = Longitude.mean - 1.3, y = Latitude.mean,label = "Men er Roue" ), stat = 'unique',
  #           size = 3,color="black",fontface = "bold") +
  geom_text(data = filter(Table_sites,Code_point_Libelle == "Ouest Loscolo"),
            aes(x = Longitude.mean + 1.3, y = Latitude.mean,label = "Ouest Loscolo" ), stat = 'unique',
            size = 3,color="black",fontface = "bold") +
  
  geom_text(aes(x = 2.75, y = 48,label = "FRANCE" ), stat = 'unique', 
            size = 7,color="burlywood4",fontface = "bold") +
  
  geom_text(aes(x = -3.5, y = 46, label = 'Océan Atlantique'), stat = 'unique', 
            size = 4,color = "lightblue4",fontface="italic",angle=-45) +
  geom_text(aes(x = -3.2, y = 50, label = 'Manche'), stat = 'unique', 
            size = 4,color = "lightblue4",fontface="italic",angle=28) +
  
  scale_fill_discrete(type = pheno_palette2, guide = 'none') +
  theme(panel.grid.major = element_line(color = 'gray10', size = .25),
        panel.background = element_rect(fill = '#BBD4F2'),panel.border = element_blank(),
        legend.background = element_rect(fill = "white"),legend.key = element_rect(fill = "white"),
        legend.text = element_text(size=8,color = "black"),
        legend.title = element_text(size = 8, color = "black",face="bold"),
        plot.background = element_rect(color = "white", linewidth = 0.001),
        legend.position = "bottom",plot.title = element_text(size=10,color ="black"),
        axis.text = element_text(size=8),axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))

# Save the map
# ggsave('Plots/REPHY/MAP_DinoPhenology_video.tiff', dpi = 600, width = 200, height = 200,
# units = 'mm', compression = 'lzw')