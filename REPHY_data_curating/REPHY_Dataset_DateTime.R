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
                         fileEncoding = "ISO-8859-1")

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
