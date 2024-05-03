#### Maps for the REPHY sites studied in Dinophysis phenology
## (Code adapted with permission from Jean-Yves DIAS)

# V. POCHIC
# 2024-05-03

#### Packages and functions ####
library(ggmap)
library(ggplot2)
library(ggthemes)
library(raster)
library(rasterVis)
library(tidyverse)
library(maps)
library(mapdata)
library(readr)
library(vctrs)


#### Import data ####

# Table with phytoplankton counts
Table_phyto_taxon <- read.csv2('Table1_phyto_taxon.csv', fileEncoding = "ISO-8859-1")

# Select only sampling sites used in our analysis
Table_phyto_geo <- filter(Table_phyto_taxon, Code_point_Libelle %in% 
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
  # We want only 1 row for each sampling site
  distinct(Code_point_Libelle, .keep_all = TRUE) %>%
  # Only keep the geographic variables (relevant)
  select(c('Code_point_Libelle', 'lon', 'lat')) %>%
  # Make the sampling site a factor
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  # Make sure latitude and longitude are coded as numeric
  mutate(across(c('lon', 'lat'), ~ as.numeric(.)))

# Make sure that the sampling site factor appears in the desired order
Table_phyto_geo <- Table_phyto_geo %>%
mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                        'Antifer ponton pétrolier', 'Cabourg',
                                        'Men er Roue', 'Ouest Loscolo',
                                        'Le Cornard', 'Auger',
                                        'Arcachon - Bouée 7', 'Teychan bis',
                                        'Parc Leucate 2', 'Bouzigues (a)',
                                        'Sète mer', 'Diana centre'))

# Clear useless data tables
rm(Table_phyto_taxon)

#### The map ####

# Coherent color palette
pheno_palette12 <- c('red3', 'orangered1', 'dodgerblue4', 'dodgerblue1', 
                     'chartreuse4', 'chartreuse2','goldenrod3', 'darkgoldenrod1',
                     'darkorchid4', 'darkorchid1', 'firebrick1', 'deeppink2'
)

# Making the map
Worldmap <- map_data('worldHires')


ggplot(Table_phyto_geo) + 
  geom_polygon(data = Worldmap, aes(x = long, y = lat, group = group), 
                                       fill = "gray", color = 'gray10', size = .25) +
  coord_fixed(xlim = c(-5.5, 9.5), ylim = c(41, 51.5), ratio = 1.4) +
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)') +
  theme_gdocs()+
  # All the points
  geom_point(aes(x = lon, y = lat, color = Code_point_Libelle), size = 8)+
  # And the text accompanying each point
  # Bassin d'Arcachon
  geom_text(data = filter(Table_phyto_geo, Code_point_Libelle == "Arcachon - Bouée 7"),
            aes(x = lon + 1.17, y = lat - .2, label = "Arcachon - Bouée 7" ), 
            stat = 'unique', 
            size = 3, color = "black", fontface = "bold") +
  geom_text(data = filter(Table_phyto_geo, Code_point_Libelle == "Teychan bis"),
            aes(x = lon + 1.17, y = lat,label = "Teychan bis" ), 
            stat = 'unique', 
            size = 3, color = "black", fontface = "bold") +
  # Bretagne Sud
  geom_text(data = filter(Table_phyto_geo, Code_point_Libelle == "Men er Roue"),
            aes(x = lon - 1.3, y = lat,label = "Men er Roue" ), 
            stat = 'unique', size = 3, color = "black", fontface = "bold") +
  geom_text(data = filter(Table_phyto_geo, Code_point_Libelle == "Ouest Loscolo"),
            aes(x = lon + 1.3, y = lat,label = "Ouest Loscolo" ), 
            stat = 'unique', size = 3,color = "black", fontface = "bold") +
  # Pertuis charentais
  geom_text(data = filter(Table_phyto_geo, Code_point_Libelle == "Le Cornard"),
            aes(x = lon + 1.1, y = lat,label = "Le Cornard" ), 
            stat = 'unique', size = 3, color = "black", fontface = "bold") +
  geom_text(data = filter(Table_phyto_geo, Code_point_Libelle == "Auger"),
            aes(x = lon + .8, y = lat,label = "Auger" ), 
            stat = 'unique', size = 3, color = "black", fontface = "bold") +
  # Baie de Seine
  geom_text(data = filter(Table_phyto_geo, Code_point_Libelle == "Antifer ponton pétrolier"),
            aes(x = lon + 2, y = lat,label = "Antifer ponton pétrolier" ), 
            stat = 'unique', size = 3, color = "black", fontface = "bold") +
  geom_text(data = filter(Table_phyto_geo,Code_point_Libelle == "Cabourg"),
            aes(x = lon + 1, y = lat,label = "Cabourg" ), 
            stat = 'unique', size = 3, color = "black", fontface = "bold") +
  # Mediterranee
  geom_text(data = filter(Table_phyto_geo, Code_point_Libelle == "Parc Leucate 2")[1,],
            aes(x = lon - 1.5, y = lat-0.11,label = "Parc Leucate 2" ), 
            stat = 'unique', 
            size = 3, color="black", fontface = "bold") +
  geom_text(data = filter(Table_phyto_geo, Code_point_Libelle == "Bouzigues (a)")[1,],
            aes(x = lon - 1.5, y = lat+0.15,label = "Bouzigues (a)" ), 
            stat = 'unique', size = 3, color = "black", fontface = "bold") +
  geom_text(data = filter(Table_phyto_geo, Code_point_Libelle == "Sète mer"),
            aes(x = lon - 1.1, y = lat-0.05,label = "Sète mer" ), 
            stat = 'unique', size = 3, color = "black", fontface = "bold") +
  geom_text(data = filter(Table_phyto_geo, Code_point_Libelle == "Diana centre")[1,],
            aes(x = lon - 1.4, y = lat,label = "Diana centre" ), 
            stat = 'unique', size = 3, color = "black", fontface = "bold") +
  # Text for the marine areas of the study
  geom_text(aes(x = -3.5, y = 46, label = 'Atlantic Ocean'), 
            size = 4, color = "lightblue4", fontface="italic", angle=-45) +
  geom_text(aes(x = -3.2, y = 50, label = 'English Channel'), 
            size = 4, color = "lightblue4", fontface="italic", angle=28) +
  geom_text(aes(x = 5, y = 42, label = 'Mediterranean Sea'), 
            size = 4, color = "lightblue4", fontface="italic", angle=8) +
  
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  theme(panel.grid.major = element_line(color = 'gray10', linewidth = .25), 
        panel.grid.minor = NULL, panel.ontop = FALSE,
        panel.background = element_rect(fill = 'lightblue1'))

# Save the map
# ggsave('DinoPhenology_map.tiff', height = 200, width = 200, units = 'mm',
#        dpi = 600, compression = 'lzw')
