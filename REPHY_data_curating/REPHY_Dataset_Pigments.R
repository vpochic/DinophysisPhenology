#### Script REPHY : handling pigment coincentration data
# This data was extracted from the Quadrige DB separately on 2024/08/27
# It contains information about pigment concentrations measured with HPLC
# as part of the REPHY network, between 2007 and 2022-23

# V. POCHIC 2024/08/29

### Required packages ####

library(tidyverse)
library(paletteer)
library(cmocean)


#### Import data ####

# Open the file with validated entries
DataREPHY_pigments <- read.csv2('Data_REPHY_pigments_2007-2022_valid.csv')
# Careful, this one is not encoded in ISO-8859-1 like the rest

# Everything seems fine.

#### Formatting the dataframe with better column names ####

colnames(DataREPHY_pigments)[which(names(DataREPHY_pigments) == "Lieu...Libellé")] <- "Code_point_Libelle"
colnames(DataREPHY_pigments)[which(names(DataREPHY_pigments) == "Passage...Date")] <- "Date"
colnames(DataREPHY_pigments)[which(names(DataREPHY_pigments) == "Résultat...Unité...Libellé")] <- "Mesure_Unite"
colnames(DataREPHY_pigments)[which(names(DataREPHY_pigments) == "Résultat...Unité...Symbole")] <- "Mesure_Symbole"
colnames(DataREPHY_pigments)[which(names(DataREPHY_pigments) == "Résultat...Valeur.quantitative")] <- "Valeur_mesure"
colnames(DataREPHY_pigments)[which(names(DataREPHY_pigments) == "Prélèvement...Immersion...Valeur")] <- "Profondeur.metre"
colnames(DataREPHY_pigments)[which(names(DataREPHY_pigments) == "Prélèvement...Niveau.de.prélèvement...Libellé")] <- "Prelevement.niveau"
colnames(DataREPHY_pigments)[which(names(DataREPHY_pigments) == "Résultat...Paramètre...Code")] <- "Code.parametre"
colnames(DataREPHY_pigments)[which(names(DataREPHY_pigments) == "Résultat...Paramètre...Libellé")] <- "Parametre"
colnames(DataREPHY_pigments)[which(names(DataREPHY_pigments) == "Résultat...Niveau.de.qualité...Libellé")] <- "Qualite.resultat"
colnames(DataREPHY_pigments)[which(names(DataREPHY_pigments) == "Résultat...Date.de.validation")] <- "Date.validation"
colnames(DataREPHY_pigments)[which(names(DataREPHY_pigments) == "Prélèvement...Niveau.de.qualité...Libellé")] <- "Qualite.prelevement"
colnames(DataREPHY_pigments)[which(names(DataREPHY_pigments) == "Prélèvement...Identifiant")] <- "ID.interne.prelevement"
colnames(DataREPHY_pigments)[which(names(DataREPHY_pigments) == "Passage...Identifiant")] <- "ID.interne.passage"

#### Curate table to keep only desired variables ####
DataREPHY_pigments_select <- DataREPHY_pigments %>% # DataREPHY_pigments_select
  # Getting rid of some duplicates in the data identified by LERN
  # the duplicates have no validation date
  filter(Date.validation != '') %>%
  dplyr::select(c('Code_point_Libelle', 'Date',
                  'Prelevement.niveau', 'Profondeur.metre', 
                  'Code.parametre',
                  'Mesure_Unite', 'Mesure_Symbole', 'Valeur_mesure', 
                  'Qualite.prelevement', 'Qualite.resultat', 
                  'ID.interne.prelevement', 'ID.interne.passage',
                  'Date.validation')) %>%
# A better date format
  mutate(Date = dmy(Date)) %>%
  # Create Year and Month
  mutate(Month = month(Date, label = F)) %>%
  mutate(Year = year(Date)) %>%
  # create a calendar day variable
  mutate(Day = as.numeric(yday(Date))) %>%
  # Create a 'fortnight' variable to match the sampling frequency
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2)) %>%
  # If no Year is indicated, discard
  filter(!is.na(Year))

#### Tidying table structure ####

# We want our data in wide format with one row for each sample and 1 column
# for each pigment
Table_pigments_select <- DataREPHY_pigments_select %>%
  # group
  group_by(Code_point_Libelle, Date, Year, Month,
           Prelevement.niveau, Profondeur.metre,
           ID.interne.prelevement, ID.interne.passage) %>%
  # filter out one row because it has a duplicate measurement for Chl a
  # (CHLOROA for Géfosse on 2020/01/27)
  # (We keep the one with the latest validation date)
  filter(!((ID.interne.prelevement == 61845366) 
           & Date.validation == '14/04/2021')) %>%
  pivot_wider(names_from = Code.parametre, values_from = Valeur_mesure)

# Check spatial and temporal extent of the dataset
unique(Table_pigments_select$Code_point_Libelle)
unique(Table_pigments_select$Year)

# Save the curated table
# write.csv2(Table_pigments_select, 'Table_pigments_2007-2022.csv',
#            row.names = FALSE)

#### Plotting pigment time series ####

# Define a color palette
pheno_palette13 <- c(# Pas de Calais
                     'black', 'sienna4', 'tan3',
                     # Baie de Seine
                     'red3', 'red', 'orangered',
                     # Other Normandy
                     '#642C3A', '#DEB1CC', '#FC4D6B',
                     # les Hébihens, Concarneau large
                     '#0A1635', '#2B4561', 
                     # Men er Roue, Ouest Loscolo
                     '#2156A1', '#5995E3')

### Alloxanthin ####
## All data points ####

# Computing the ratio of alloxanthin by chlorophyll a
Table_pigments_plot <- Table_pigments_select %>%
  mutate(Allo.ratio = Allo/CHLOROA) %>%
  # select only surface measurements (there are only surface measurements
  # in the table)
  filter(Prelevement.niveau == 'Surface (0-1m)') %>%
  # Reorder Code_point_Libelle so sites appear in the desired order
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Dunkerque',
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Ouistreham 1 mille',
                                          'Cabourg', 'Géfosse', 'Digue de Querqueville',
                                          'Champeaux', 'les Hébihens',
                                          'Concarneau large',
                                          'Men er Roue', 'Ouest Loscolo'))

### Alloxanthin concentration
# Plot the alloxanthin concentration as function of date
ggplot(Table_pigments_plot)+
  geom_point(aes(x=Date, y=Allo, color  = Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  scale_color_discrete(type = pheno_palette13, guide = 'none') +
  theme_classic() +
  labs(y="Alloxanthin (µg/L)", x="Date")

# Note : there are some NAs in the dataset

# Save the plot (for sending it to Maud :)
# ggsave('REPHY_timeseries_Alloxanthin.tiff', height = 150,
#        width = 200, unit = 'mm', compression = 'lzw')

### Alloxanthin/chlorophyll a ratio
# Plot the alloxanthin/chlorophyll a ratio as function of date
ggplot(Table_pigments_plot)+
  geom_point(aes(x=Date, y=Allo.ratio, color  = Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  scale_color_discrete(type = pheno_palette13, guide = 'none') +
  theme_classic() +
  labs(y="Alloxanthin/Chlorophyll a", x="Date")

# Note : there are some NAs in the dataset

# Save the plot (for sending it to Maud :)
# ggsave('REPHY_timeseries_Alloxanthin_ratio.tiff', height = 150,
#        width = 200, unit = 'mm', compression = 'lzw')

# The seasonal signal for alloxanthin is probably driven in part by 
# chlorophyll a (i.e., by the bloom)

### Alloxanthin concentration in one composite year
# Plot the alloxanthin concentration as function of calendar day
ggplot(Table_pigments_plot)+
  geom_point(aes(x=Day, y=Allo, color  = Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  scale_color_discrete(type = pheno_palette13, guide = 'none') +
  theme_classic() +
  labs(y="Alloxanthin concentration (µg/L)", x="Calendar day")

# Save the plot (for sending it to Maud :)
# ggsave('REPHY_seasonality_Alloxanthin.tiff', height = 150,
#        width = 200, unit = 'mm', compression = 'lzw')

### Alloxanthin/chlorophyll a ratio in one composite year
# Plot the alloxanthin/chlorophyll a ratio as function of calendar day
ggplot(Table_pigments_plot)+
  geom_point(aes(x=Day, y=Allo.ratio, color  = Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  scale_color_discrete(type = pheno_palette13, guide = 'none') +
  theme_classic() +
  labs(y="Alloxanthin/Chlorophyll a", x="Calendar day")

# Save the plot (for sending it to Maud :)
# ggsave('REPHY_seasonality_Alloxanthin_ratio.tiff', height = 150,
#        width = 200, unit = 'mm', compression = 'lzw')

## Cropping some data points ####

# Based on the plots in the previous section, we crop out the highest 
# data points when they alter the visualisation

Table_pigments_plot_crop <- Table_pigments_plot %>%
  filter(ifelse(
    # Dunkerque
    Code_point_Libelle == 'Point 1 Dunkerque' & Allo < 0.5, TRUE,
    # Boulogne
    ifelse(Code_point_Libelle == 'Point 1 Boulogne' & Allo.ratio < 0.1, TRUE,
    # Antifer
    ifelse(Code_point_Libelle == 'Antifer ponton pétrolier' & Allo.ratio < 0.3, TRUE,
    # Géfosse
    ifelse(Code_point_Libelle == 'Géfosse' & Allo.ratio < 0.3, TRUE,
    # Digue de Querqueville
    ifelse(Code_point_Libelle == 'Digue de Querqueville' & Allo.ratio < 0.2, TRUE,
    # les Hébihens
    ifelse(Code_point_Libelle == 'les Hébihens' & Allo.ratio < 2, TRUE,
    # all the other points have no outliers
    ifelse(Code_point_Libelle %in% c('At so', 'Ouistreham 1 mille', 'Cabourg',
                                     'Champeaux', 'Concarneau large', 
                                     'Men er Roue', 'Ouest Loscolo'), TRUE, 
           FALSE))))))))

### Now the same plots but with the cropped tops

### Alloxanthin concentration cropped
# Plot the alloxanthin concentration as function of date
ggplot(Table_pigments_plot_crop)+
  geom_point(aes(x=Date, y=Allo, color  = Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  scale_color_discrete(type = pheno_palette13, guide = 'none') +
  theme_classic() +
  labs(y="Alloxanthin (µg/L)", x="Date")

# Note : there are some NAs in the dataset

# Save the plot (for sending it to Maud :)
# ggsave('REPHY_timeseries_Alloxanthin_crop.tiff', height = 150,
#        width = 200, unit = 'mm', compression = 'lzw')

### Alloxanthin/chlorophyll a ratio
# Plot the alloxanthin/chlorophyll a ratio as function of date
ggplot(Table_pigments_plot_crop)+
  geom_point(aes(x=Date, y=Allo.ratio, color  = Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  scale_color_discrete(type = pheno_palette13, guide = 'none') +
  theme_classic() +
  labs(y="Alloxanthin/Chlorophyll a", x="Date")

# Note : there are some NAs in the dataset

# Save the plot (for sending it to Maud :)
# ggsave('REPHY_timeseries_Alloxanthin_ratio_crop.tiff', height = 150,
#        width = 200, unit = 'mm', compression = 'lzw')

# The seasonal signal for alloxanthin is probably driven in part by 
# chlorophyll a (i.e., by the bloom)

### Alloxanthin concentration in one composite year
# Plot the alloxanthin concentration as function of calendar day
ggplot(Table_pigments_plot_crop)+
  geom_point(aes(x=Day, y=Allo, color  = Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  scale_color_discrete(type = pheno_palette13, guide = 'none') +
  theme_classic() +
  labs(y="Alloxanthin concentration (µg/L)", x="Calendar day")

# Save the plot (for sending it to Maud :)
# ggsave('REPHY_seasonality_Alloxanthin.tiff', height = 150,
#        width = 200, unit = 'mm', compression = 'lzw')

### Alloxanthin/chlorophyll a ratio in one composite year
# Plot the alloxanthin/chlorophyll a ratio as function of calendar day
ggplot(Table_pigments_plot_crop)+
  geom_point(aes(x=Day, y=Allo.ratio, color  = Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  scale_color_discrete(type = pheno_palette13, guide = 'none') +
  theme_classic() +
  labs(y="Alloxanthin/Chlorophyll a", x="Calendar day")

# Save the plot (for sending it to Maud :)
# ggsave('REPHY_seasonality_Alloxanthin_ratio_crop.tiff', height = 150,
#        width = 200, unit = 'mm', compression = 'lzw')

### Chlorophyll a ####
## All data points ####

### chlorophyll a concentration
# Plot the chlorophyll a concentration as function of date
ggplot(Table_pigments_plot)+
  geom_point(aes(x=Date, y=CHLOROA, color  = Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  scale_color_discrete(type = pheno_palette13, guide = 'none') +
  theme_classic() +
  labs(y="Chlorophyll a (µg/L)", x="Date")

# Note : there are some NAs in the dataset

# Save the plot (for sending it to Maud :)
# ggsave('REPHY_timeseries_CHLOROA.tiff', height = 150,
#        width = 200, unit = 'mm', compression = 'lzw')

### Chlorophyll a concentration in one composite year
# Plot the chlorophyll a concentration as function of calendar day
ggplot(Table_pigments_plot)+
  geom_point(aes(x=Day, y=CHLOROA, color  = Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  scale_color_discrete(type = pheno_palette13, guide = 'none') +
  theme_classic() +
  labs(y="Chlorophyll a concentration (µg/L)", x="Calendar day")

# Save the plot (for sending it to Maud :)
# ggsave('REPHY_seasonality_CHLOROA.tiff', height = 150,
#        width = 200, unit = 'mm', compression = 'lzw')



#### Summarising by fortnight ####

Table_pigments_fortnightly <- Table_pigments_plot %>%
  group_by(Code_point_Libelle, Fortnight) %>%
  summarise(# Allox
    Allo.med = median(Allo, na.rm = TRUE), 
    Allo.mean = mean(Allo, na.rm = TRUE),
    # Allox/CHLOROA
    Allo.ratio.med = median(Allo.ratio, na.rm = TRUE), 
    Allo.ratio.mean = mean(Allo.ratio, na.rm = TRUE),
    # Chlorophyll a, ignore NAs
    CHLOROA.med = median(CHLOROA, na.rm = TRUE), 
    CHLOROA.mean = mean(CHLOROA, na.rm = TRUE),
    .groups = 'keep') %>%
  # filter out Fortnight 27 as there isn't enough measurements to calculate
  # a reliable median
  filter(Fortnight < 27)

# Write that down!
# write.csv2(Table_pigments_fortnightly, 'Table_hydro_fortnightly_20240905.csv',
#            row.names = FALSE, fileEncoding = "ISO-8859-1")

# As boxplots
# Alloxanthin
ggplot(Table_pigments_plot)+
  geom_boxplot(aes(x=as.factor(Fortnight), y=Allo, color  = Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  scale_color_discrete(type = pheno_palette13, guide = 'none') +
  theme_classic() +
  labs(y="Alloxanthin by fortnight (µg/L)", x="Fortnight")

# Chlorophyll a
ggplot(Table_pigments_plot)+
  geom_boxplot(aes(x=as.factor(Fortnight), y=CHLOROA, color  = Code_point_Libelle), 
               alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  scale_color_discrete(type = pheno_palette13, guide = 'none') +
  theme_classic() +
  labs(y="Chlorophyll a by fortnight (µg/L)", x="Fortnight")

# As heatmaps
# Now let's plot that as a heatmap

# Alloxanthin
ggplot(Table_pigments_fortnightly, aes(x = Fortnight, y = 1, fill = Allo.med)) +
  geom_tile()  +
  scale_fill_distiller(palette = 'RdBu', direction = -1) +
  facet_wrap(facets = c('Code_point_Libelle'))

ggplot(Table_pigments_fortnightly, aes(x = Fortnight, y = 1, fill = Allo.mean)) +
  geom_tile()  +
  scale_fill_distiller(palette = 'RdBu', direction = -1) +
  facet_wrap(facets = c('Code_point_Libelle'))

# Chlorophyll a
ggplot(Table_pigments_fortnightly, aes(x = Fortnight, y = 1, fill = CHLOROA.med)) +
  geom_tile()  +
  scale_fill_cmocean(name = 'algae') +
  facet_wrap(facets = c('Code_point_Libelle'))

ggplot(Table_pigments_fortnightly, aes(x = Fortnight, y = 1, fill = CHLOROA.mean)) +
  geom_tile()  +
  scale_fill_cmocean(name = 'algae') +
  facet_wrap(facets = c('Code_point_Libelle'))
