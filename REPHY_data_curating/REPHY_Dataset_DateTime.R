#### Script REPHY : getting dates for hydrological model data extraction
# V. POCHIC 2024/05/24

# Packages

library(tidyverse)
library(stringr)
library(FactoMineR)
library(factoextra)
library(vegan)
library(paletteer)
library(mgcv)
library(gratia)
library(scatterpie)


#### Import data ####

# Open with the correct file encoding allows to preserve accents
DataREPHY_MA <- read.csv2('REPHY_Manche_Atlantique_1987-2022.csv', fileEncoding = "ISO-8859-1")
DataREPHY_Med <- read.csv2('REPHY_Med_1987-2022.csv', fileEncoding = "ISO-8859-1")

# Creating a vector to import all columns as characters
classes_char <- rep('character', 56)
# CAREFUL! ALL 2023 DATA LACK TEMPERATURE! Problem seen with Maud on 2024/01/09. To be continued.
DataREPHY_2023 <- read.csv2('Extraction SEANOE_REPHY_phyto-Hydro_Manche_Atl-Med validé 19122023.csv', 
                            fileEncoding = "ISO-8859-1", colClasses = classes_char)

# Merge the first 2 datasets
DataREPHY <- bind_rows(DataREPHY_MA, DataREPHY_Med) %>%
  # This line allows to remove the problematic row that separates the hydro and phyto datasets
  filter(Passage...Mois != 'Passage : Mois')

# Binding the 2 datasets
DataREPHY_8723 <- bind_rows(DataREPHY, DataREPHY_2023)


### Load tables used to enrich the dataset
# Load table "Zones_marines"
ZM <- read.csv('Zones_marines.csv', sep = ';', header = TRUE)

# Load table "Liste_phylum.classe"
PhyClasse <- read.csv('Liste_phylum.classe_REPHY.csv', sep =';', header = TRUE, fileEncoding = 'ISO-8859-1')

#### Formatting the dataframe with better column names ####

# Extracting the numeric code for the ZM
Table1 <- DataREPHY_8723 %>%
  mutate(ZM_Quadrige_Numero = as.numeric(str_extract(Lieu.de.surveillance...Entité.de.classement...Libellé, '[:alnum:]+')))

# Change column names (to match ZM)
colnames(Table1)[which(names(Table1) == "Lieu.de.surveillance...Mnémonique")] <- "Code_point_Mnemonique"
colnames(Table1)[which(names(Table1) == "Lieu.de.surveillance...Libellé")] <- "Code_point_Libelle"
colnames(Table1)[which(names(Table1) == "Passage...Date")] <- "Date"
colnames(Table1)[which(names(Table1) == "Coordonnées.passage...Coordonnées.minx")] <- "lon"
colnames(Table1)[which(names(Table1) == "Coordonnées.passage...Coordonnées.miny")] <- "lat"
colnames(Table1)[which(names(Table1) == "Résultat...Libellé.unité.de.mesure.associé.au.quintuplet")] <- "Mesure_Unite"
colnames(Table1)[which(names(Table1) == "Résultat...Symbole.unité.de.mesure.associé.au.quintuplet")] <- "Mesure_Symbole"
colnames(Table1)[which(names(Table1) == "Résultat...Nom.du.taxon.référent")] <- "Taxon"
colnames(Table1)[which(names(Table1) == "Résultat...Libellé.du.groupe.de.taxon")] <- "Groupe_Taxon"
colnames(Table1)[which(names(Table1) == "Résultat...Valeur.de.la.mesure")] <- "Valeur_mesure"
colnames(Table1)[which(names(Table1) == "Prélèvement...Immersion")] <- "Profondeur.metre"
colnames(Table1)[which(names(Table1) == "Prélèvement...Niveau")] <- "Prelevement.niveau"
colnames(Table1)[which(names(Table1) == "Résultat...Code.paramètre")] <- "Code.parametre"
colnames(Table1)[which(names(Table1) == "Résultat...Libellé.paramètre")] <- "Parametre"
colnames(Table1)[which(names(Table1) == "Résultat...Niveau.de.qualité")] <- "Qualite.resultat"
colnames(Table1)[which(names(Table1) == "Prélèvement...Niveau.de.qualité")] <- "Qualite.prelevement"
colnames(Table1)[which(names(Table1) == "Passage...Service.saisisseur...Libellé")] <- "Service.saisie"
colnames(Table1)[which(names(Table1) == "Passage...Heure")] <- "Heure"
colnames(Table1)[which(names(Table1) == "Résultat...Service.analyste...Code")] <- "Service.analyse"
colnames(Table1)[which(names(Table1) == "Prélèvement...Service.préleveur...Code")] <- "Service.prelevement"
colnames(Table1)[which(names(Table1) == "Prélèvement...Identifiant.interne")] <- "ID.interne.prelevement"
colnames(Table1)[which(names(Table1) == "Passage...Identifiant.interne")] <- "ID.interne.passage"

#### Curate table to keep only desired variables ####
Table1 <- Table1 %>%
  dplyr::select(c('ZM_Quadrige_Numero', 'Code_point_Mnemonique', 'Code_point_Libelle', 'Date', 
                  'Heure', 'lon', 'lat', 'Mesure_Unite', 'Mesure_Symbole', 'Taxon', 'Valeur_mesure', 
                  'Prelevement.niveau', 'Profondeur.metre', 'Code.parametre', 'Parametre', 
                  'Qualite.prelevement', 'Qualite.resultat', 'ID.interne.prelevement', 'ID.interne.passage'))

# Modifying date format so that it gives the year and month, and getting rid of rows with no Year value
Table1 <- Table1 %>%
  # Compute the date
  mutate(Date = dmy(Date)) %>%
  # Convert date back to character for date and time concatenation
  mutate(Date = as.character(Date)) %>%
  # Concatenate Date and Heure to have an "unambiguous" format
  mutate(DateTime = paste(Date, Heure, sep = ' ')) %>%
  mutate(DateTime = ymd_hms(DateTime)) %>%
  # Exclude those dates that do not have an hour
  filter(is.na(DateTime) == FALSE) %>%
  # Convert Date back to date format
  mutate(Date = ymd(Date)) %>%
  # Create Year, Month, Hour and Minute variables
  mutate(Day = day(Date)) %>%
  mutate(Month = month(Date, label = F)) %>%
  mutate(Year = year(Date)) %>%
  # Compute the hour and minutes
  mutate(Time_hour = hour(DateTime)) %>%
  mutate(Time_minute = minute(DateTime)) %>%
  filter(!is.na(Year))

### Now let's clean this  a bit

#### Cleaning up ####

# We only want samplings corresponding to the FLORTOT protocol
Table1_select <- Table1 %>%
  filter(Code.parametre == 'FLORTOT')

# Spread the measurements

Table1_select <- Table1_select %>%
  select(c(# Geographical data
    'Code_point_Libelle', 'Code_point_Mnemonique', 'lon', 'lat',
    # Time data
    'DateTime', 'Date', 'Day', 'Month', 'Year', 'Time_hour', 'Time_minute',
    # Sample ID and data
    'ID.interne.passage', 'Qualite.prelevement', 'Code.parametre', 
    'Valeur_mesure', 'Prelevement.niveau',
    )) %>%
  filter(Code_point_Libelle %in% c( # Mer du Nord
    'Point 1 Boulogne', 'At so',
    # Baie de Seine
    'Antifer ponton pétrolier', 'Cabourg',
    # Bretagne Nord
    'Loguivy', 'les Hébihens',
    # Bretagne Sud
    'Men er Roue', 'Ouest Loscolo',
    # Pertuis charentais
    'Auger', 'Le Cornard',
    # Arcachon
    'Arcachon - Bouée 7', 'Teychan bis',
    # Mediterranée
    'Parc Leucate 2', 'Bouzigues (a)', 'Sète mer',
    'Diana centre', '22B - Toulon gde rade')) %>%
  # Only keep (roughly) surface level measurements
  filter(Prelevement.niveau == 'Surface (0-1m)' | Prelevement.niveau == '2 metres' |
           Prelevement.niveau == 'de 3 à 5 metres') %>%
  group_by(Code_point_Libelle, lon, lat, Year, Month, Date, Day,
           DateTime, Time_hour, Time_minute,
           ID.interne.passage, Code.parametre) %>%
  # Mutate Valeur_mesure as numeric so we avoid a warning issued by the mean()
  # function
  mutate(Valeur_mesure = as.numeric(Valeur_mesure)) %>%
  # The pivot operation's only purpose here is to make samples unique,
  # it is not valid scientifically speaking (averaging several phytoplankton taxa)
  summarise(Valeur_mesure = mean(Valeur_mesure), .groups = 'keep') %>%
  pivot_wider(names_from = Code.parametre, values_from = Valeur_mesure)

# Select only relevant columns to have a clean table of measurements
Table_dates <- Table1_select %>%
  dplyr::select(c( #Space info
    'Code_point_Libelle', 'lon', 'lat',
    # Time info
    'Date', 'Year', 'Month', 'Day',
    'DateTime', 'Time_hour', 'Time_minute',
    # Sample info
    'ID.interne.passage'))

# Filter this to fit the time period of our study
Table_dates <- filter(Table_dates, Year >= 2007 & Year <= 2022)

# Write that down WRITE THAT DOWN!
write.csv2(Table_dates, 'Table_dates_REPHY_select_20240524.csv', 
           row.names = FALSE)

# Make a little sample
Table_dates_sample <- filter(Table_dates,
                             (Code_point_Libelle == 'Ouest Loscolo' |
                               Code_point_Libelle == 'Men er Roue') &
                               Year == 2018 &
                               Month >= 4 &
                               Month <= 6)

# Write the sample down in a table
write.csv2(Table_dates_sample, 'Table_dates_REPHY_sample_20240524.csv', 
           row.names = FALSE)
