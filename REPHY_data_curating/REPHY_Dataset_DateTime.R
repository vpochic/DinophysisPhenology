#### Script REPHY : getting dates for hydrological model data extraction
# V. POCHIC 2024/08/19

# Packages

library(tidyverse)
library(stringr)
library(FactoMineR)
library(factoextra)
library(vegan)


#### Import data ####

# Import the table corresponding to DinoZeros 
# (table of all Dinophysis counts including zeros, in the 16 sites)
Season_Dino <- read.csv2('Season_Dino.csv', header = TRUE, 
                         fileEncoding = 'ISO-8859-1')

#### Curate data (mainly to give an appropriate time format) ####

# Modifying date format so that it gives the year and month, and getting rid of rows with no Year value
Table_extractions <- Season_Dino %>%
  # Compute the date
  mutate(Date = ymd(Date)) %>%
  # Convert date back to character for date and time concatenation
  mutate(Date = as.character(Date)) %>%
  # Concatenate Date and Heure to have an "unambiguous" format
  mutate(DateTime = paste(Date, Heure, sep = ' ')) %>%
  mutate(DateTime = ymd_hms(DateTime)) %>%
  # A warning followed by a little inspection tells us there is 1 
  # event with no hour indicated. It is discarded at this point.
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

# Only select the relevant (space and time)

Table_dates_select <- Table_extractions %>%
  select(c(# Geographical data
    'Code_point_Libelle', 'Longitude', 'Latitude',
    # Time data
    'DateTime', 'Date', 'Day', 'Month', 'Year', 'Time_hour', 'Time_minute',
    # Identification code
    'ID.interne.passage'
    ))

# Here we do a little test: we concatenate 'Code_point_Libelle' (site code) and
# 'DateTime, and check that the number of UNIQUE values in this column matches
# the length of the whole dataset (= number of ID.interne.passage, that are 
# attributed to each sampling)
# If true, this means each observation in the dataset corresponds to one 
# observation in real life (no duplicates)
# We have to do this procedure because there are genuine duplicates in 
# DateTime : 2 different sites sampled simultaneously

test <- Table_dates_select %>%
  mutate(DateTime = as.character(DateTime)) %>%
  mutate(test_column = paste(Code_point_Libelle, DateTime, sep = ' '))

length(unique(test$test_column))
length(Table_dates_select$ID.interne.passage)
# It's alright, we're good
rm(test)
  

# Write that down WRITE THAT DOWN!
write.csv2(Table_dates_select, 'Table_dates_REPHY_select_20240819.csv',
           row.names = FALSE)

# Make a little sample
Table_dates_sample <- filter(Table_dates_select,
                             (Code_point_Libelle == 'Ouest Loscolo' |
                               Code_point_Libelle == 'Men er Roue') &
                               Year == 2018 &
                               Month >= 4 &
                               Month <= 6)

# Write the sample down in a table
# write.csv2(Table_dates_sample, 'Table_dates_REPHY_sample_20240819.csv',
#            row.names = FALSE)

#### A hydrology table for checking model results ####

# Import the hydrology table (with all depth levels measured)
Table_hydro <- read.csv2('Table1_hydro_models.csv', header = TRUE,
                         fileEncoding = 'ISO-8859-1') %>%
  # Filter to only keep the interesting sites
  filter(Code_point_Libelle %in% 
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
  # Rename lon and lat into Longitude and Latitude, and Heure into Hour
  rename(Longitude = lon, Latitude = lat, Hour = Heure) %>%
  ### Date and Time changes
  # Compute the date
  mutate(Date = ymd(Date)) %>%
  # Convert date back to character for date and time concatenation
  mutate(Date = as.character(Date)) %>%
  # Concatenate Date and Heure to have an "unambiguous" format
  mutate(DateTime = paste(Date, Hour, sep = ' ')) %>%
  mutate(DateTime = ymd_hms(DateTime)) %>%
  # A warning followed by a little inspection tells us there are 2
  # events with no hour indicated. They are discarded at this point.
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
  filter(!is.na(Year)) %>%
  # Only events for which Temperature and Salinity were properly recorded
  filter(is.na(TEMP) == FALSE) %>%
  filter(is.na(SALI) == FALSE) %>%
  # Select some columns
  select(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,47,
         52,53,54,55)

# check the number of unique IDs (ID.interne.prelevement)
length(unique(Table_hydro$ID.interne.prelevement))
# It's equal to the number of rows in the table: hurray!

# Write the table down and send it to the modellers!
# write.csv2(Table_hydro, 'Table_hydro_REPHY_multidepth_20240826.csv', row.names = FALSE)

# Note : the number of unique IDs in this table is smaller than in the first one
# because it doesn't include the Mediterranean sites.
# When they are excluded, the table with only surface measurements is indeed
# smaller (4209 unique IDs)
