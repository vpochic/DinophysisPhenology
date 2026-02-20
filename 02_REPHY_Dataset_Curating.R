#### Script REPHY: curating data ###
# Author: V. POCHIC 
# Last modif: 2026/02/20

## General information ####

## Description: 
# This script takes REPHY datasets downloaded directly from the REPHY dataset
# on SEANOE (https://doi.org/10.17882/47248) and processes them to versions used 
# in our analysis.

## Files required:
# Data/REPHY/REPHY_Manche_Atlantique_1987-2022.csv
# Data/REPHY/REPHY_Med_1987-2022.csv
# (Downloaded directly from the REPHY dataset on SEANOE)
# Data/REPHY/Data_REPHY_pigments_2007-2022_valid.csv
# (Obtained from the Q2 database that contains all REPHY data)

## Outputs:
## Data  files:
# Data/REPHY_outputs/...
# Table_phyto_taxon.csv ; Table_hydro.csv
# Season_Dino.csv ; Season_Meso.csv ; Season_Allo.csv
## Figure:
# Figure S2 - Seasonality of Dinophysis taxa

### Required packages ####

library(tidyverse)

####-----------------------------------------------------------------------####
#### Import data ####

# Open with the correct file encoding allows to preserve accents
DataREPHY_MA <- read.csv2('Data/REPHY/REPHY_Manche_Atlantique_1987-2022.csv', fileEncoding = "ISO-8859-1")
DataREPHY_Med <- read.csv2('Data/REPHY/REPHY_Med_1987-2022.csv', fileEncoding = "ISO-8859-1")

# Bind MA and Med datasets
DataREPHY <- bind_rows(DataREPHY_MA, DataREPHY_Med) %>%
  # This line allows to remove the problematic row that separates the hydro and phyto datasets
  filter(Passage...Mois != 'Passage : Mois')

# Clear some space
rm(DataREPHY_MA, DataREPHY_Med)

#### Formatting the dataframe with better column names ####

# Change column names
colnames(DataREPHY)[which(names(DataREPHY) == "Lieu.de.surveillance...Mnémonique")] <- "Code_point_Mnemonique"
colnames(DataREPHY)[which(names(DataREPHY) == "Lieu.de.surveillance...Libellé")] <- "Code_point_Libelle"
colnames(DataREPHY)[which(names(DataREPHY) == "Passage...Date")] <- "Date"
colnames(DataREPHY)[which(names(DataREPHY) == "Coordonnées.passage...Coordonnées.minx")] <- "lon"
colnames(DataREPHY)[which(names(DataREPHY) == "Coordonnées.passage...Coordonnées.miny")] <- "lat"
colnames(DataREPHY)[which(names(DataREPHY) == "Résultat...Libellé.unité.de.mesure.associé.au.quintuplet")] <- "Mesure_Unite"
colnames(DataREPHY)[which(names(DataREPHY) == "Résultat...Symbole.unité.de.mesure.associé.au.quintuplet")] <- "Mesure_Symbole"
colnames(DataREPHY)[which(names(DataREPHY) == "Résultat...Nom.du.taxon.référent")] <- "Taxon"
colnames(DataREPHY)[which(names(DataREPHY) == "Résultat...Libellé.du.groupe.de.taxon")] <- "Groupe_Taxon"
colnames(DataREPHY)[which(names(DataREPHY) == "Résultat...Valeur.de.la.mesure")] <- "Valeur_mesure"
colnames(DataREPHY)[which(names(DataREPHY) == "Prélèvement...Immersion")] <- "Profondeur.metre"
colnames(DataREPHY)[which(names(DataREPHY) == "Prélèvement...Niveau")] <- "Prelevement.niveau"
colnames(DataREPHY)[which(names(DataREPHY) == "Résultat...Code.paramètre")] <- "Code.parametre"
colnames(DataREPHY)[which(names(DataREPHY) == "Résultat...Libellé.paramètre")] <- "Parametre"
colnames(DataREPHY)[which(names(DataREPHY) == "Résultat...Libellé.méthode")] <- "Methode.analyse"
colnames(DataREPHY)[which(names(DataREPHY) == "Résultat...Niveau.de.qualité")] <- "Qualite.resultat"
colnames(DataREPHY)[which(names(DataREPHY) == "Prélèvement...Niveau.de.qualité")] <- "Qualite.prelevement"
colnames(DataREPHY)[which(names(DataREPHY) == "Passage...Service.saisisseur...Libellé")] <- "Service.saisie"
colnames(DataREPHY)[which(names(DataREPHY) == "Passage...Heure")] <- "Heure"
colnames(DataREPHY)[which(names(DataREPHY) == "Résultat...Service.analyste...Code")] <- "Service.analyse"
colnames(DataREPHY)[which(names(DataREPHY) == "Prélèvement...Service.préleveur...Code")] <- "Service.prelevement"
colnames(DataREPHY)[which(names(DataREPHY) == "Prélèvement...Identifiant.interne")] <- "ID.interne.prelevement"
colnames(DataREPHY)[which(names(DataREPHY) == "Passage...Identifiant.interne")] <- "ID.interne.passage"

#### Curate table to keep only desired variables ####
DataREPHY <- DataREPHY %>%
  dplyr::select(c('Code_point_Mnemonique', 'Code_point_Libelle', 'Date', 
                  'Heure', 'lon', 'lat', 'Mesure_Unite', 'Mesure_Symbole', 'Taxon', 'Valeur_mesure', 'Methode.analyse',
                  'Prelevement.niveau', 'Profondeur.metre', 'Code.parametre', 'Parametre', 
                  'Qualite.prelevement', 'Qualite.resultat', 'ID.interne.prelevement', 'ID.interne.passage'))

# Modifying date format so that it gives the year and month, and getting rid of rows with no Year value
DataREPHY <- DataREPHY %>%
  mutate(Date = dmy(Date)) %>%
  # modifies the date format
  mutate(Day = day(Date)) %>%
  mutate(Month = month(Date, label = F)) %>%
  mutate(Year = year(Date)) %>%
  filter(!is.na(Year))

# Transform the measured values into numbers
# Caution! The decimal separator is a comma, 
# we need first to transform it to a full stop to avoid fuckery
DataREPHY$Valeur_mesure <-  str_replace_all(DataREPHY$Valeur_mesure, ',', '.')
DataREPHY <- DataREPHY %>%
  mutate(Valeur_mesure = as.numeric(Valeur_mesure))

#### Tidying data table structure ####

# As of now, the hydrology measurements and phytoplankton counts are mixed up
# in the table
# Separate the table into 2 : 1 for hydrology and the other for phytoplankton
# For hydro measurements, the 'Taxon' field is empty
Table_hydro <- DataREPHY %>%
  filter(Taxon == "")

Table_phyto <- DataREPHY %>%
  filter(Taxon != "")

### Manipulating the phyto table ###

# Applying quality control and selecting only FLORTOT
Table_phyto <- Table_phyto %>%
  filter(Code.parametre == 'FLORTOT') %>%
  # Keeping only surface (or near surface) sampling
  filter(Prelevement.niveau == 'Surface (0-1m)' | Prelevement.niveau == '2 metres' |
           Prelevement.niveau == 'de 3 à 5 metres') %>%
  # Grouping by taxon
  group_by(Code_point_Libelle, lon, lat, Year, Month, Date, Heure,
           Code.parametre, ID.interne.passage, Taxon) %>%
  summarise(Comptage = sum(Valeur_mesure), .groups = 'keep')

# Save the phytoplankton taxon table
# write.csv2(Table_phyto, 'Data/REPHY_outputs/Table_phyto_taxon.csv',
# row.names = FALSE, fileEncoding = "ISO-8859-1")

### Manipulating the hydro table ###
# Spread the hydrology measurements

Table_hydro <- Table_hydro %>%
  select(c('ID.interne.passage', 'Qualite.resultat', 'Code.parametre', 'Valeur_mesure', 'Methode.analyse', 
           'Prelevement.niveau','Profondeur.metre', 'Date', 'Day',
           'Month', 'Year', 'Heure', 'Code_point_Libelle', 'Code_point_Mnemonique',
           'lon', 'lat')) %>%
  # Only surface or near surface samples
  filter(Prelevement.niveau == 'Surface (0-1m)' | Prelevement.niveau == '2 metres' |
           Prelevement.niveau == 'de 3 à 5 metres') %>%
  # There is  something shifty here, regarding multiple CHLOROA measurements at certain stations,
  # made with different methods. This is due to additional programs of measuring pigments with
  # HPLC is some sampling sites.
  #To fix that, we're going to create another category, named CHLOROA_HPLC, for measures of 
  # chl a in HPLC, that come on top of the spectrometric method in many cases
  mutate(Code.parametre = ifelse(
    Methode.analyse == 'Chromatographie liquide - pigments phytoplanctoniques (Van Heukelem et Thomas 2001)' &
      Code.parametre == 'CHLOROA',
    'CHLOROA_HPLC', Code.parametre)) %>%
  # Now grouping by sampling ID and parameter
  group_by(Code_point_Libelle, lon, lat, Year, Month, Date, ID.interne.passage, Code.parametre) %>%
  # When multiple values available (a few rare cases with close values, and some
  # even rarer full duplicates), we take the average of the measurements
  summarise(Valeur_mesure = mean(Valeur_mesure), .groups = 'keep') %>%
  # And then pivot wider
  pivot_wider(names_from = Code.parametre, values_from = Valeur_mesure)

# Save the hydrology table
# write.csv2(Table_hydro, 'Data/REPHY_outputs/Table_hydro.csv',
# row.names = FALSE, fileEncoding = "ISO-8859-1")



####-----------------------------------------------------------------------####
#### Compute smaller specific datasets for organisms/pigments of interest ####

# In this section, we'll create 3 datasets (Season_Dino, Season_Meso and
# Season_Allo) for our phenology analysis of Dinophysis, Mesodinium and
# Alloxanthin.

#### Import data ####

# Import table with phytoplankton counts
Table_phyto_taxon <- read.csv2('Data/REPHY_outputs/Table_phyto_taxon.csv', 
                               fileEncoding = "ISO-8859-1", header = TRUE)

# Select only sampling sites of interest from 2007 --> 2022
Table_phyto_select <- filter(Table_phyto_taxon, Code_point_Libelle %in% 
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
                             'Arcachon - Bouée 7', 'Teychan bis',
                             # Mediterranée
                             'Parc Leucate 2', 'Bouzigues (a)', 'Sète mer',
                             'Diana centre')
) %>%
  # Only years in the period [2007-2022]
  filter(Year >= 2007 & Year <= 2022)

# Import table with pigment data
# This is a small additional dataset necessary because some pigment data
# are lacking in the official release
DataREPHY_pigments <- read.csv2('Data/REPHY/Data_REPHY_pigments_2007-2022_valid.csv')
# Careful, this one is not encoded in ISO-8859-1 like the rest

#### Dinophysis and Mesodinium ####

## More zeros

### The goal here is to create a table where all absences of detection are 
# listed as zeros
Table_phyto_zeros <- Table_phyto_select %>%
  # 'Comptage' is the cell count
  pivot_wider(names_from = Taxon, values_from = Comptage, values_fill = 0) %>%
  # Create Latitude and Longitude variables for selection (because lon and lat
  # introduce a bunch of unwanted taxa )
  mutate(Longitude = as.numeric(lon)) %>%
  mutate(Latitude = as.numeric(lat)) %>%
  # We keep only Dinophysis, Mesodinium and the metadata columns
  select(starts_with('Dinophysis') | starts_with('Mesodinium') |
           contains(c('Code_point_Libelle', 'Year', 'Month',
                      'Date', 'Heure', 'Code.parametre', 
                      'Longitude', 'Latitude', 
                      'ID.interne.passage'))) %>%
  # create a day of the year variable and get the Date in a 'Date' format
  mutate(Date = ymd(Date)) %>%
  mutate(Day = as.numeric(yday(Date)))

  #---# Dinophysis #---#

Season_Dino <- Table_phyto_zeros %>%
  #--- Specificity in data from Arcachon sampling sites ---#
  # Note that for the 2 sites in Arcachon, counts are sometimes done in 100 mL, 
  # so we need to filter out those counts.
  # For this, we go and look for any count values IN ANY DINOPHYSIS TAXON that
  # cannot correspond to 10 mL counts (i.e., not multiples of 100)
  # We replace these values with 0. Note that for a given date, we can have
  # one Dinophysis species counted in 10 mL and another in 100 mL
  # That's why we have to proceed taxon by taxon
  mutate(across(.cols = c('Dinophysis', 'Dinophysis + phalacroma', 'Dinophysis acuta',
                          'Dinophysis acuminata', 'Dinophysis caudata', 'Dinophysis tripos',
                          'Dinophysis sacculus', 'Dinophysis fortii'), 
                .fns = ~ ifelse(. %% 100 != 0, 0, .))) %>%
  # We then create a count variable for Dinophysis as a genus
  mutate(Dinophysis_genus = rowSums(across(contains('Dinophysis')))) %>%
  #--- Create a 'true count' variable ---#
  # this variable corresponds to the number of cells that were actually
  # counted by the operator (in 10 mL). This variable will follow a Poisson
  # distribution, contrary to 'Dinophysis_genus' because the conversion from
  # 10 mL to 1L (*100) in the original data prevents some intermediate values 
  # (e.g., 150 cells.L-1)
  mutate(true_count = Dinophysis_genus/100) %>%
  # Filter out exceptionnaly high counts (>500 cells observed in 10 mL)
  # this represents 3 events (2 in Antifer and 1 in Cabourg)
  filter(true_count < 500)

## Let's save that precious precious file
# write.csv2(Season_Dino, 'Data/REPHY_outputs/Season_Dino.csv', row.names = FALSE,
#            fileEncoding = "ISO-8859-1")

#---# Mesodinium #---#

Season_Meso <- Table_phyto_zeros %>%
  #--- Restricted spatio-temporal range ---#
  # In the REPHY protocol, counting Mesodinium is not mandatory. Therefore, the
  # data on this genus is heterogeneous. We will select only the 4 sampling sites
  # and the time period for which we are sure that the identification and counting
  # were done systematically and precisely.
  filter(Code_point_Libelle %in% c('Antifer ponton pétrolier', 'Cabourg',
                                   'Ouest Loscolo', 'Men er Roue')) %>%
  filter(Year >= 2016) %>%
  # We then create a count variable for Dinophysis as a genus
  mutate(Mesodinium_genus = rowSums(across(contains('Mesodinium')))) %>%
  #--- Create a 'true count' variable ---#
  # this variable corresponds to the number of cells that were actually
  # counted by the operator (in 10 mL). This variable will follow a Poisson
  # distribution, contrary to 'Dinophysis_genus' because the conversion from
  # 10 mL to 1L (*100) in the original data prevents some intermediate values 
  # (e.g., 150 cells.L-1)
  mutate(true_count = Mesodinium_genus/100) %>%
  # Filter out exceptionnaly high counts (>500 cells observed in 10 mL)
  # this represents 1 event (1 in Cabourg)
  filter(true_count < 500)

## Again, let's save that precious file
# write.csv2(Season_Meso, 'Data/REPHY_outputs/Season_Meso.csv', row.names = FALSE,
#            fileEncoding = "ISO-8859-1")

#### Alloxanthin ####

# Now we're going to take the alloxanthin concentration data from the pigments
# table. Here we go.

# First, we need to curate this table in the same fashion we did for the big
# REPHY tables.

#--- Better column names ---#

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

#--- Curate table to keep only desired variables ---#
DataREPHY_pigments <- DataREPHY_pigments %>%
  # Getting rid of some duplicates in the data identified by IFREMER/LERN
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
  # create a day of the year variable
  mutate(Day = as.numeric(yday(Date))) %>%
  # If no Year is indicated, discard
  filter(!is.na(Year))

#--- Tidying table structure ---#

# We want our data in wide format with one row for each sample and 1 column
# for each pigment
Table_pigments <- DataREPHY_pigments %>%
  # group
  group_by(Code_point_Libelle, Date, Year, Month,
           Prelevement.niveau, Profondeur.metre,
           ID.interne.prelevement, ID.interne.passage) %>%
  # filter out one row because it has a duplicate measurement for Chl a
  # (CHLOROA for Géfosse on 2020/01/27)
  # (We keep the one with the latest validation date)
  filter(!((ID.interne.prelevement == 61845366) 
           & Date.validation == '14/04/2021')) %>%
  # pivot wider
  pivot_wider(names_from = Code.parametre, values_from = Valeur_mesure)

#--- Specific dataset for alloxanthin ---#

Season_Allo <- Table_pigments %>%
  #--- Sam restricted spatio-temporal coverage as Mesodinium ---#
  filter(Code_point_Libelle %in% c('Antifer ponton pétrolier', 'Cabourg',
                                   'Ouest Loscolo', 'Men er Roue')) %>%
  filter(Year >= 2016) %>%
  # ungroup
  ungroup() %>%
  # keep only alloxanthin and chlorophyll a variables
  select(c('Code_point_Libelle', 'Year', 'Month', 'Day',
           'Date', 'ID.interne.passage',
           'Allo', 'CHLOROA'))

## Write the thing down
# write.csv2(Season_Allo, 'Data/REPHY_outputs/Season_Allo.csv', row.names = FALSE,
#            fileEncoding = "ISO-8859-1")
  

#### Proportion of Dinophysis taxa per fortnight -> Fig. S2 ####

# An interesting information in the dataset is when each Dinophysis species occur
# (when identified) and in which propotion compared to all Dinophysis counts.

# We'll go through the trouble of starting again from Season_Dino as we need to
# keep the correction made for the Arcachon sampling sites
Season_Dino_longer <- Season_Dino %>%
  # Pivot longer with values to 'Comptage' (cell count)
  pivot_longer(cols = contains(c('Dinophysis')), 
               names_to = 'Taxon', 
               values_to = 'Comptage') %>%
  # Remove Dinophysis_genus
  filter(Taxon != 'Dinophysis_genus') %>%
  # Remove these values that were calculated on Dinophysis_genus
  select(-c('true_count')) %>%
  # Create a 'Fortnight' variable as sampling occurs fortnightly
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2))

# A plot-ready version of the table
Season_Dino_fortnight <- Season_Dino_longer %>%
  # Group data by fortnight to get a 1-year seasonality
  group_by(Code_point_Libelle, Fortnight, Taxon) %>%
  # Summarise by fortnight
  summarise(Comptage = sum(Comptage), .groups = 'keep')
  
  
#--- Now plot ---#

## Color palette
taxo_palette8 <- c(# Dinophysis  
  '#11203E',
  # Dinophysis acuminata-complex (acuminata, sacculus, fortii)
  # Though I'm not really sure anymore than D. fortii is part of the complex...
  '#435E7B', '#7F96B6', '#BBD4F2',
  # Dinophysis acuta
  '#8D6456',
  # Dinophysis caudata and tripos
  '#FBA823', '#FBB646',
  # Dinophysis + phalacroma
  '#FF6448'
)

# Changing factor order so sampling sites appear in the order we want
Season_Dino_fortnight <- Season_Dino_fortnight %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Loguivy', 'les Hébihens',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

# And we do the same for the Dinophysis taxa
Season_Dino_fortnight <- Season_Dino_fortnight %>%
  mutate(Taxon = as_factor(Taxon)) %>%
  mutate(Taxon = fct_relevel(Taxon,
                             'Dinophysis', 'Dinophysis acuminata',
                             'Dinophysis sacculus', 'Dinophysis fortii',
                             'Dinophysis acuta', 'Dinophysis caudata',
                             'Dinophysis tripos',
                             'Dinophysis + phalacroma')) %>%
  # Add a 'true count' variable that identifies the number of cells actually 
  # observed and counted
  mutate(true_count = Comptage/100)

## Plot command
ggplot(Season_Dino_fortnight, aes(x = Fortnight, y = true_count, fill = Taxon)) +
  geom_col() +
  facet_wrap(facets = 'Code_point_Libelle', scales = 'free_y') +
  # axis
  scale_x_continuous(breaks = c(1, 8, 16, 20, 26),
                     labels = c('1', '8', '16', '20', '26')) +
  # Labels
  labs(title = "Observations of Dinophysis taxa (2007-2022)",
       y = "Sum of observed Dinophysis cells",
       x = 'Fortnight',
       fill = 'Taxon:'
  ) +
  # Legend labels
  scale_fill_discrete(labels = c('Dinophysis', 'D. acuminata',
                                 'D. sacculus', 'D. fortii',
                                 'D. acuta', 'D. caudata',
                                 'D. tripos', 'Dinophysis + Phalacroma'),
                      # Color palette
                      type = taxo_palette8,
                      # guide legend
                      guide = guide_legend(nrow = 2)) +
  # Theme
  theme_classic() +
  theme(# title and facets text
    title = element_text(size = 11),
    strip.text = element_text(size = 7.5),
    # axes text
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7.5),
    # legend
    legend.key.height = unit(.75, 'cm'),
    legend.key.width = unit(.6, 'cm'),
    legend.title = element_text(size=9),
    legend.text = element_text(size=8),
    legend.position = 'bottom',
    legend.background = element_rect(color = 'black', linewidth = .25)
  )

## Save this awesome plot as Fig S2
# ggsave('Plots/REPHY/FigS2_Taxa_16sites_Dinophenology_fortnight.tiff', dpi = 300, width = 164, height = 164,
#        units = 'mm', compression = 'lzw')

### The plot shows the prevalence of cells identified at the genus level
# ('Dinophysis'), that explains our choice to work at the genus level


####-------------------------- End of script ------------------------------####