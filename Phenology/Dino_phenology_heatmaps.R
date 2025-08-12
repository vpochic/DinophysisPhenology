#### Temperature and salinity heatmaps for Dinophysis phenology ##
## V. POCHIC
# 2025/08/12
# --> Now with only the sites selected for phenology analysis
# (The Channel-Atlantic coast) and with Chl a (curated values - no hplc)

### Packages ####

library(tidyverse)
library(viridis)
library(paletteer)
library(RColorBrewer)
library(cmocean)

### Import data ####
Table_hydro <- read.csv2('Table1_hydro_select_20250812.csv', header = TRUE, fileEncoding = "ISO-8859-1")

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
  # Create a 'fortnight' variable to match the sampling frequency
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2))

### Importing stratification data ####

# This is just a bit trickier
# We need to import the data that sits in separate files then stick it together

### Stratification index averaged on the 14 days before the date
# Read all the files in a folder at once

# Create a list of file names (without the path)
list_SI_14d <- list.files(path = 'GAMAR_csv/SI_14d_before/', full.names = FALSE)
# Read all files in a dataframe with map_dfr
Stratif_14d <- map_dfr(list_SI_14d, ~ read_csv(
                                      file.path('GAMAR_csv/SI_14d_before/', .),
                                      id = 'name')) %>%
  
  # Fix a problem in column names : replace all spaces in column names with '_'
  rename_with(~ gsub(" ","_", .x), contains(" ")) %>%
  # Create 'Code_point_Libelle', the site id that is missing in the data
  mutate(Code_point_Libelle = ifelse(
    # Boulogne
    grepl('BOULOGNE', name, fixed = TRUE) == TRUE, 'Point 1 Boulogne',
    # At so
    ifelse(grepl('AT_SO', name, fixed = TRUE) == TRUE, 'At so',
    # Antifer
    ifelse(grepl('ANTIFER', name, fixed = TRUE) == TRUE, 'Antifer ponton pétrolier',
    # Cabourg
    ifelse(grepl('CABOURG', name, fixed = TRUE) == TRUE, 'Cabourg',
    # Les Hébihens
    ifelse(grepl('LES_HÉBIHENS', name, fixed = TRUE) == TRUE, 'les Hébihens',
    # Loguivy
    ifelse(grepl('LOGUIVY', name, fixed = TRUE) == TRUE, 'Loguivy',
    # Men er Roue
    ifelse(grepl('MEN_ER_ROUE', name, fixed = TRUE) == TRUE, 'Men er Roue',
    # Ouest Loscolo
    ifelse(grepl('LOSCOLO', name, fixed = TRUE) == TRUE, 'Ouest Loscolo',
    # Le Cornard
    ifelse(grepl('CORNARD', name, fixed = TRUE) == TRUE, 'Le Cornard',
    # Auger
    ifelse(grepl('AUGER', name, fixed = TRUE) == TRUE, 'Auger',
    # Arcachon
    ifelse(grepl('ARCACHON', name, fixed = TRUE) == TRUE, 'Arcachon - Bouée 7',
    # Teychan bis
    'Teychan bis')))))))))))) %>%
  # Get rid of the awful 'name' variable
  select(-c('name')) %>%
  # Remove duplicates
  distinct(Target_Date, Closest_Date, Code_point_Libelle,
           .keep_all = TRUE)

### Stratification index on the date of the REPHY sampling
# Exactly the same

# Create a list of file names (without the path)
list_SI_dates <- list.files(path = 'GAMAR_csv/SI_REPHY_dates/', full.names = FALSE)
# Read all files in a dataframe with map_dfr
Stratif_dates <- map_dfr(list_SI_dates, ~ read_csv(
  file.path('GAMAR_csv/SI_REPHY_dates/', .),
  id = 'name')) %>%
  # Fix a problem in column names : replace all spaces in column names with '_'
  rename_with(~ gsub(" ","_", .x), contains(" ")) %>%
  # Create 'Code_point_Libelle', the site id that is missing in the data
  mutate(Code_point_Libelle = ifelse(
    # Boulogne
    grepl('BOULOGNE', name, fixed = TRUE) == TRUE, 'Point 1 Boulogne',
    # At so
    ifelse(grepl('AT_SO', name, fixed = TRUE) == TRUE, 'At so',
           # Antifer
           ifelse(grepl('ANTIFER', name, fixed = TRUE) == TRUE, 'Antifer ponton pétrolier',
                  # Cabourg
                  ifelse(grepl('CABOURG', name, fixed = TRUE) == TRUE, 'Cabourg',
                         # Les Hébihens
                         ifelse(grepl('LES_HÉBIHENS', name, fixed = TRUE) == TRUE, 'les Hébihens',
                                # Loguivy
                                ifelse(grepl('LOGUIVY', name, fixed = TRUE) == TRUE, 'Loguivy',
                                       # Men er Roue
                                       ifelse(grepl('MEN_ER_ROUE', name, fixed = TRUE) == TRUE, 'Men er Roue',
                                              # Ouest Loscolo
                                              ifelse(grepl('LOSCOLO', name, fixed = TRUE) == TRUE, 'Ouest Loscolo',
                                                     # Le Cornard
                                                     ifelse(grepl('CORNARD', name, fixed = TRUE) == TRUE, 'Le Cornard',
                                                            # Auger
                                                            ifelse(grepl('AUGER', name, fixed = TRUE) == TRUE, 'Auger',
                                                                   # Arcachon
                                                                   ifelse(grepl('ARCACHON', name, fixed = TRUE) == TRUE, 'Arcachon - Bouée 7',
                                                                          # Teychan bis
                                                                          'Teychan bis')))))))))))) %>%
  # Get rid of the awful 'name' variable
  select(-c('name')) %>%
  # Remove duplicates
  distinct(Target_Date, Closest_Date, Code_point_Libelle,
           .keep_all = TRUE)

# Bind both datasets (making it easier for processing)
Table_stratif <- left_join(Stratif_14d, Stratif_dates, 
                           by = c('Target_Date', 
                                  'Closest_Date', 'Code_point_Libelle'),
                           suffix = c('',''))

# Save the stratification dataset
# write.csv2(Table_stratif, 'Stratif_index_GAMAR_12sites.csv', row.names = FALSE,
# fileEncoding = 'ISO-8859-1')

##### Plotting temperatures #####

# First, reorder the factor so that sampling sites appear in the desired order
Table_hydro_select <- Table_hydro_select %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis'#,
                                          # 'Parc Leucate 2', 'Bouzigues (a)',
                                          # 'Sète mer', 'Diana centre'
                                          ))
  
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

##### Plotting Chl a ####

### We make another table for chl a
Table_chla_select <- filter(Table_hydro, Code_point_Libelle %in% 
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
  # Only events for which Chl a was properly recorded
  filter(is.na(CHLOROA) == FALSE) %>%
  filter(CHLOROA > 0) %>%
  # Change to date format
  mutate(Date = ymd(Date)) %>%
  mutate(Month = month(Date)) %>%
  # Create a 'fortnight' variable to match the sampling frequency
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2))

### Now let's plot !

# First, reorder the factor so that sampling sites appear in the desired order
Table_chla_select <- Table_chla_select %>%
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
ggplot(Table_chla_select, aes(x = Date, y = CHLOROA, color = Code_point_Libelle)) +
  geom_point()  +
  scale_color_discrete(type = pheno_palette12) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y')

# Everything's fine
# Now let's plot that as a heatmap

ggplot(Table_chla_select, aes(x = Date, y = 1, fill = CHLOROA)) +
  geom_tile()  +
  facet_wrap(facets = c('Code_point_Libelle'))
# Strange

# As dot plot by fortnight
ggplot(Table_chla_select, aes(x = Fortnight,
                              # Convert chl a to log scale 
                              # (dealing with extreme values)
                              y = log10(CHLOROA), 
                              color = Code_point_Libelle)) +
  geom_point()  +
  scale_color_discrete(type = pheno_palette12) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y')

# As box plot by fortnight
ggplot(Table_chla_select) +
  geom_boxplot(aes(x = factor(Fortnight),
                   # Convert chl a to log scale 
                   # (dealing with extreme values)
                   y = log10(CHLOROA), 
                   color = Code_point_Libelle))  +
  geom_point(aes(x = Fortnight, y = log10(CHLOROA), color = Code_point_Libelle),
             alpha = .25)  +
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  # Labels
  labs(title = c(expression(paste('Chl ', italic('a')))),
       y = c(expression(paste('log[Chl ', italic('a'),' concentration]'))),
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

# As box plot by fortnight but with real values (not log-transformed)
ggplot(Table_chla_select) +
  geom_boxplot(aes(x = factor(Fortnight),
                   y = CHLOROA, 
                   color = Code_point_Libelle),
               outliers = FALSE)  +
  # geom_point(aes(x = Fortnight, y = CHLOROA, color = Code_point_Libelle),
  #            alpha = .25)  +
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  # Labels
  labs(title = c(expression(paste('Chl ', italic('a')))),
       y = c(expression(paste('Chl ', italic('a'),' concentration (mg.L'^'-1',')'))),
       x = 'Fortnight',
       color = 'Sampling site'
  ) +
  # Scale y axis
  # scale_y_continuous(limits = c(0,30)) +
  # Divide this in little plots for each sampling site
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
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
# ggsave('Chla_12sites_Dinophenology.tiff', dpi = 300, width = 300, height = 200,
#        units = 'mm', compression = 'lzw')

##### Plotting Stratification Index ####

### We make another table for SI
Table_stratif_select <- Table_stratif %>%
  # Change to date format
  mutate(Date = ymd_hms(Closest_Date)) %>%
  mutate(Month = month(Date)) %>%
  # Create a 'fortnight' variable to match the sampling frequency
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2))

### Now let's plot !

# First, reorder the factor so that sampling sites appear in the desired order
Table_stratif_select <- Table_stratif_select %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis'))

# Plot aesthetics
pheno_palette12 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646')

# Little plot on stratif (dates)
ggplot(Table_stratif_select, aes(x = Date, y = Stratification_Index, 
                                 color = Code_point_Libelle)) +
  geom_point()  +
  scale_color_discrete(type = pheno_palette12) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y')

# Everything's fine
# Now let's plot that as a heatmap

ggplot(Table_stratif_select, aes(x = Date, y = 1, fill = Stratification_Index)) +
  geom_tile()  +
  facet_wrap(facets = c('Code_point_Libelle'))
# Strange

# As dot plot by fortnight
ggplot(Table_stratif_select, aes(x = Fortnight,
                              # Convert chl a to log scale 
                              # (dealing with extreme values)
                              y = Stratification_Index, 
                              color = Code_point_Libelle)) +
  geom_point()  +
  scale_color_discrete(type = pheno_palette12) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y')

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
  # Chlorophyll a, ignore NAs
  CHLOROA.med = median(CHLOROA, na.rm = TRUE), 
  CHLOROA.mean = mean(CHLOROA, na.rm = TRUE),
  .groups = 'keep') %>%
  # filter out Fortnight 27 as there isn't enough measurements to calculate
  # a reliable median
  filter(Fortnight < 27)

# Write that down!
# write.csv2(Table_hydro_fortnightly, 'Table_hydro_fortnightly_20250812.csv',
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

# Chlorophyll a
ggplot(Table_hydro_fortnightly, aes(x = Fortnight, y = 1, fill = CHLOROA.med)) +
  geom_tile()  +
  scale_fill_cmocean(name = 'algae') +
  facet_wrap(facets = c('Code_point_Libelle'))

ggplot(Table_hydro_fortnightly, aes(x = Fortnight, y = 1, fill = CHLOROA.mean)) +
  geom_tile()  +
  scale_fill_cmocean(name = 'algae') +
  facet_wrap(facets = c('Code_point_Libelle'))


### Now for the stratification index

Table_stratif_fortnightly <- Table_stratif_select %>%
  group_by(Code_point_Libelle, Fortnight) %>%
  summarise(# Stratification_Index (on date)
    SI_date.med = median(Stratification_Index), 
    SI_date.mean = mean(Stratification_Index),
    # 14-Day_Average_SI (on the previous 14 days)
    SI_14d.med = median(`14-Day_Average_SI`), 
    SI_14d.mean = mean(`14-Day_Average_SI`),
    .groups = 'keep') %>%
  # filter out Fortnight 27 as there isn't enough measurements to calculate
  # a reliable median
  filter(Fortnight < 27)

# Write that down!
# write.csv2(Table_stratif_fortnightly, 'Table_stratif_fortnightly_20240923.csv',
#            row.names = FALSE, fileEncoding = "ISO-8859-1")

# Stratification index heatmaps
# On date
ggplot(Table_stratif_fortnightly, aes(x = Fortnight, y = 1, fill = SI_date.med)) +
  geom_tile()  +
  scale_fill_cmocean(name = 'tempo') +
  facet_wrap(facets = c('Code_point_Libelle'))
# 14 previous days
ggplot(Table_stratif_fortnightly, aes(x = Fortnight, y = 1, fill = SI_14d.med)) +
  geom_tile()  +
  scale_fill_cmocean(name = 'tempo') +
  facet_wrap(facets = c('Code_point_Libelle'))
