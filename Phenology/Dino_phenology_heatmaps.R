#### Temperature and salinity heatmaps for Dinophysis phenology ##
## V. POCHIC
# 2024/04/30

### Packages ####

library(tidyverse)
library(viridis)
library(paletteer)
library(RColorBrewer)

### Import data ####
Table_hydro <- read.csv2('Table1_hydro_select.csv', header = TRUE, fileEncoding = "ISO-8859-1")

# Select only stations of interest from 2007 --> 2022
Table_hydro_select <- filter(Table_hydro, Code_point_Libelle %in% 
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
  # Only after 2007
  filter(Year >= 2007) %>%
  # Only events for which Temperature and Salinity were properly recorded
  filter(is.na(TEMP) == FALSE) %>%
  filter(is.na(SALI) == FALSE) %>%
  # Change to date format
  mutate(Date = ymd(Date)) %>%
  mutate(Month = month(Date)) %>%
  # Create a 'fortnight' variable to match the sampling frequency
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2))

##### Plotting temperatures #####

# First, reorder the factor so that sampling sites appear in the desired order
Table_hydro_select <- Table_hydro_select %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))
  
# Plot aesthetics
pheno_palette12 <- c('red3', 'orangered1', 'dodgerblue4', 'dodgerblue1', 
                     'chartreuse4', 'chartreuse2','goldenrod3', 'darkgoldenrod1',
                     'darkorchid4', 'darkorchid1', 'firebrick1', 'deeppink2'
)

# Little plot
ggplot(Table_hydro_select, aes(x = Date, y = TEMP, color = Code_point_Libelle)) +
  geom_point()  +
  scale_color_discrete(type = pheno_palette12) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y')

# Everything's fine
# Now let's plot that as a heatmap

ggplot(Table_hydro_select, aes(x = Date, y = 1, fill = TEMP)) +
  geom_tile()  +
  facet_wrap(facets = c('Code_point_Libelle'))
# Strange

# As dot plot by fortnight
ggplot(Table_hydro_select, aes(x = Fortnight, y = TEMP, color = Code_point_Libelle)) +
  geom_point()  +
  scale_color_discrete(type = pheno_palette12) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y')

# As box plot by fortnight
ggplot(Table_hydro_select) +
  geom_boxplot(aes(x = factor(Fortnight), y = TEMP, color = Code_point_Libelle
               ))  +
  geom_point(aes(x = Fortnight, y = TEMP, color = Code_point_Libelle),
             alpha = .25)  +
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  # Labels
  labs(title = 'Sea surface temperature',
       y = "Temperature (°C)",
       x = 'Fortnight',
       color = 'Sampling site'
  ) +
  # Divide this in little plots for each sampling site
  facet_wrap(facets = c('Code_point_Libelle')) +
  theme(
    panel.background = element_rect(fill = 'transparent', 
                                    linewidth = .3, 
                                    color = 'grey10'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = .2, color = 'grey60'),
    panel.grid.minor.y = element_line(linewidth = .1, color = 'grey60'),
    axis.text.x = element_text(size = 5)
  )

# Save this really nice plot
# ggsave('SST_12sites_Dinophenology.tiff', dpi = 300, width = 300, height = 200,
#        units = 'mm', compression = 'lzw')

##### Plotting salinity #####

# First, reorder the factor so that sampling sites appear in the desired order
Table_hydro_select <- Table_hydro_select %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

# Plot aesthetics
pheno_palette12 <- c('red3', 'orangered1', 'dodgerblue4', 'dodgerblue1', 
                     'chartreuse4', 'chartreuse2','goldenrod3', 'darkgoldenrod1',
                     'darkorchid4', 'darkorchid1', 'firebrick1', 'deeppink2'
)

# Little plot
ggplot(Table_hydro_select, aes(x = Date, y = SALI, color = Code_point_Libelle)) +
  geom_point()  +
  scale_color_discrete(type = pheno_palette12) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y')

# Everything's fine
# Now let's plot that as a heatmap

ggplot(Table_hydro_select, aes(x = Date, y = 1, fill = TEMP)) +
  geom_tile()  +
  facet_wrap(facets = c('Code_point_Libelle'))
# Strange

# As dot plot by fortnight
ggplot(Table_hydro_select, aes(x = Fortnight, y = SALI, color = Code_point_Libelle)) +
  geom_point()  +
  scale_color_discrete(type = pheno_palette12) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y')

# As box plot by fortnight
ggplot(Table_hydro_select) +
  geom_boxplot(aes(x = factor(Fortnight), y = SALI, color = Code_point_Libelle
  ))  +
  geom_point(aes(x = Fortnight, y = SALI, color = Code_point_Libelle),
             alpha = .25)  +
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  # Labels
  labs(title = 'Sea surface salinity',
       y = "Salinity",
       x = 'Fortnight',
       color = 'Sampling site'
  ) +
  # Divide this in little plots for each sampling site
  facet_wrap(facets = c('Code_point_Libelle')) +
  theme(
    panel.background = element_rect(fill = 'transparent', 
                                    linewidth = .3, 
                                    color = 'grey10'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = .2, color = 'grey60'),
    panel.grid.minor.y = element_line(linewidth = .1, color = 'grey60'),
    axis.text.x = element_text(size = 5)
  )

# Save this really nice plot
# ggsave('SSS_12sites_Dinophenology.tiff', dpi = 300, width = 300, height = 200,
#        units = 'mm', compression = 'lzw')

### Summarising by month ####

Table_hydro_monthly <- Table_hydro_select %>%
  group_by(Code_point_Libelle, Month) %>%
  summarise(# Temperature
    TEMP.med = median(TEMP), TEMP.mean = mean(TEMP),
    # Salinity
    SALI.med = median(SALI), SALI.mean = mean(SALI),
            .groups = 'keep')

# Now let's plot that as a heatmap

# Temperature
ggplot(Table_hydro_monthly, aes(x = Month, y = 1, fill = TEMP.med)) +
  geom_tile()  +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free') 

# Salinity
ggplot(Table_hydro_monthly, aes(x = Month, y = 1, fill = SALI.med)) +
  geom_tile()  +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free') 

### Summarising by fortnight ####

Table_hydro_fortnightly <- Table_hydro_select %>%
  group_by(Code_point_Libelle, Fortnight) %>%
  summarise(# Temperature
    TEMP.med = median(TEMP), TEMP.mean = mean(TEMP),
    # Salinity
    SALI.med = median(SALI), SALI.mean = mean(SALI),
    .groups = 'keep') %>%
  # filter out Fortnight 27 as there isn't enough measurements to calculate
  # a reliable median
  filter(Fortnight < 27)

# Write that down!
# write.csv2(Table_hydro_fortnightly, 'Table_hydro_fortnightly.csv',
#            row.names = FALSE, fileEncoding = "ISO-8859-1")

# Now let's plot that as a heatmap

# Temperature
ggplot(Table_hydro_fortnightly, aes(x = Fortnight, y = 1, fill = TEMP.med)) +
  geom_tile()  +
  scale_fill_distiller(palette = 'RdBu', direction = -1) +
  facet_wrap(facets = c('Code_point_Libelle'))

ggplot(Table_hydro_fortnightly, aes(x = Fortnight, y = 1, fill = TEMP.mean)) +
  geom_tile()  +
  scale_fill_distiller(palette = 'RdBu', direction = -1) +
  facet_wrap(facets = c('Code_point_Libelle'))

# Salinity
ggplot(Table_hydro_fortnightly, aes(x = Fortnight, y = 1, fill = SALI.med)) +
  geom_tile()  +
  scale_fill_cmocean(name = 'haline') +
  facet_wrap(facets = c('Code_point_Libelle'))

ggplot(Table_hydro_fortnightly, aes(x = Fortnight, y = 1, fill = SALI.mean)) +
  geom_tile()  +
  scale_fill_cmocean(name = 'haline') +
  facet_wrap(facets = c('Code_point_Libelle'))
