### ERA5 data: managing datasets ###

# V. POCHIC, 2025/10/09

# A script for managing ERA 5 datasets and getting them into the right format
# for further integration in our analyses.

# The data used in this script were extracted from the ERA5 atmospheric model.
# There are 4 variables from the model, associated with the time and lon/lat.

# These variables are :
# 1- ssr (Surface Solar Radiation), the solar energy flux at sea surface level 
# in (J/m2)
# 2- tcc, Total Cloud Cover, in proportion from 0 to 1
# 3- u10, U zonal component of wind (West to East), 10m above sea surface, 
# in m/s
# 4- v10, V zonal component of wind (South to North), 10m above sea surface, 
# in m/s


# Packages ####
library(tidyverse)

# Import data ####
# Define the path to get the data files

path<-c('Data/Models/ERA5/Extraction_rephy_era')

# List all files in the folders
file_names <- path %>%
  lapply(list.files,full.names=TRUE) %>%  # Apply list.files() to the folder
  unlist()            # Flatten the list into a vector


## Import and merge the different files together (this step can take some time)
# Read all the csv files into a list
dataset <- lapply(file_names,function(i){read.csv2(i,
                                                header = TRUE,
                                                fileEncoding = 'ISO-8859-1')
  }) %>%
  # Name each dataframe in the list 
  set_names(file_names) %>% 
  # Merge into a big global dataframe, while keeping the original column names
  plyr::ldply(.,function(df) as_tibble(df)) %>% 
  # Create a column to identify the sampling site (instead of folder name)
  mutate(.id = basename(sub('_era5.csv', '', .id))) %>%
  mutate(Code_point_Libelle = sub('rephy', '', .id)) %>%
  # Remove the .id column
  select(-.id)

# Works fine

## Working a bit around the dataset to make it easier to handle

# First, transform the date in a lubridate date format
dataset <- dataset %>%
  mutate(Date = ymd_hms(time)) %>%
  mutate(Hour = hour(Date)) %>%
  mutate(Day = yday(Date)) %>%
  mutate(Year = year(Date)) %>%
  mutate(Week = week(Date)) %>%
  # the fortnight variable is based on Week
  mutate(Fortnight = ceiling(Week/2))

# Then, make sure the variables are numeric
dataset <- dataset %>%
  mutate(across(c(latitude, longitude, ssr, tcc, u10, v10), 
         ~ as.numeric(.))) %>%
  # Get the names of sampling points in the right format
  mutate(Code_point_Libelle = ifelse(grepl('ANTIFER', Code_point_Libelle) == TRUE, 
                   'Antifer ponton pétrolier',
            ifelse(grepl('CABOURG', Code_point_Libelle) == TRUE, 
                   'Cabourg',
            ifelse(grepl('BOULOGNE', Code_point_Libelle) == TRUE, 
                   'Point 1 Boulogne',
            ifelse(grepl('AT_SO', Code_point_Libelle) == TRUE, 
                   'At so',
            ifelse(grepl('LOGUIVY', Code_point_Libelle) == TRUE, 
                   'Loguivy',
            ifelse(grepl('LES_H', Code_point_Libelle) == TRUE, 
                   'les Hébihens',
            ifelse(grepl('MEN_ER', Code_point_Libelle) == TRUE, 
                   'Men er Roue',
            ifelse(grepl('LOSCOLO', Code_point_Libelle) == TRUE, 
                   'Ouest Loscolo',
            ifelse(grepl('CORNARD', Code_point_Libelle) == TRUE, 
                   'Le Cornard',
            ifelse(grepl('AUGER', Code_point_Libelle) == TRUE, 
                   'Auger',
            ifelse(grepl('ARCACHON', Code_point_Libelle) == TRUE, 
                   'Arcachon - Bouée 7', 'Teychan bis'))))))))))))

# It's quite ugly, but it seems to work

# Modify the dataset for ssr to exclude night time (ssr=0)
dataset_ssr <- select(dataset, -c('tcc','u10','v10')) %>%
  filter(ssr > 0)

# And exclude ssr from the other dataset
dataset <- select(dataset, -c('ssr'))

## Aggregating data (daily, fortnightly...) ####

# u10, v10, tcc ####

# Currently the data is in hourly format. This is a bit too detailed for us.
# We will aggregate the data over longer periods (day, week, fortnight...)

## Day
# Mean
# Compute the mean of variables over daily periods
dataset_daily_mean <- dataset %>%
 # Group by day, year and sampling site
 group_by(Day, Year, Code_point_Libelle) %>%
 # Ditch other variables
 select(-c(time, Date, Week, Fortnight)) %>%
 # Take the mean of other variables
 summarise(across(c(latitude, longitude, tcc, u10, v10), 
                  ~ mean(.)))

# Median
dataset_daily_median <- dataset %>%
  # Group by day, year and sampling site
  group_by(Day, Year, Code_point_Libelle) %>%
  # Ditch other variables
  select(-c(time, Date, Week, Fortnight)) %>%
  # Take the mean of other variables
  summarise(across(c(latitude, longitude, tcc, u10, v10), 
                   ~ median(.)))


## Week
# Mean
# Compute the mean of variables over daily periods
dataset_weekly_mean <- dataset %>%
  # Group by day, year and sampling site
  group_by(Week, Year, Code_point_Libelle) %>%
  # Ditch other variables
  select(-c(time, Date, Day, Fortnight)) %>%
  # Take the mean of other variables
  summarise(across(c(latitude, longitude, tcc, u10, v10), 
                   ~ mean(.)))

# Median
dataset_weekly_median <- dataset %>%
  # Group by day, year and sampling site
  group_by(Week, Year, Code_point_Libelle) %>%
  # Ditch other variables
  select(-c(time, Date, Day, Fortnight)) %>%
  # Take the mean of other variables
  summarise(across(c(latitude, longitude, tcc, u10, v10), 
                   ~ median(.)))


## Fortnight
# Mean
dataset_fortnight_mean <- dataset %>%
  # Group by day, year and sampling site
  group_by(Fortnight, Year, Code_point_Libelle) %>%
  # Ditch other variables
  select(-c(time, Date, Week, Day)) %>%
  # Take the mean of other variables
  summarise(across(c(latitude, longitude, tcc, u10, v10), 
                   ~ mean(.)))

# Median
dataset_fortnight_median <- dataset %>%
  # Group by day, year and sampling site
  group_by(Fortnight, Year, Code_point_Libelle) %>%
  # Ditch other variables
  select(-c(time, Date, Week, Day)) %>%
  # Take the mean of other variables
  summarise(across(c(latitude, longitude, tcc, u10, v10), 
                   ~ median(.)))


# ssr ####

# Currently the data is in hourly format. This is a bit too detailed for us.
# We will aggregate the data over longer periods (day, week, fortnight...)

## Day
# Mean
# Compute the mean of variables over daily periods
dataset_ssr_daily_mean <- dataset_ssr %>%
  # Group by day, year and sampling site
  group_by(Day, Year, Code_point_Libelle) %>%
  # Ditch other variables
  select(-c(time, Date, Week, Fortnight)) %>%
  # Take the mean of other variables
  summarise(across(c(latitude, longitude, ssr), 
                   ~ mean(.)))

# Median
dataset_ssr_daily_median <- dataset_ssr %>%
  # Group by day, year and sampling site
  group_by(Day, Year, Code_point_Libelle) %>%
  # Ditch other variables
  select(-c(time, Date, Week, Fortnight)) %>%
  # Take the mean of other variables
  summarise(across(c(latitude, longitude, ssr), 
                   ~ median(.)))


## Week
# Mean
# Compute the mean of variables over daily periods
dataset_ssr_weekly_mean <- dataset_ssr %>%
  # Group by day, year and sampling site
  group_by(Week, Year, Code_point_Libelle) %>%
  # Ditch other variables
  select(-c(time, Date, Day, Fortnight)) %>%
  # Take the mean of other variables
  summarise(across(c(latitude, longitude, ssr), 
                   ~ mean(.)))

# Median
dataset_ssr_weekly_median <- dataset_ssr %>%
  # Group by day, year and sampling site
  group_by(Week, Year, Code_point_Libelle) %>%
  # Ditch other variables
  select(-c(time, Date, Day, Fortnight)) %>%
  # Take the mean of other variables
  summarise(across(c(latitude, longitude, ssr), 
                   ~ median(.)))


## Fortnight
# Mean
dataset_ssr_fortnight_mean <- dataset_ssr %>%
  # Group by day, year and sampling site
  group_by(Fortnight, Year, Code_point_Libelle) %>%
  # Ditch other variables
  select(-c(time, Date, Week, Day)) %>%
  # Take the mean of other variables
  summarise(across(c(latitude, longitude, ssr), 
                   ~ mean(.)))

# Median
dataset_ssr_fortnight_median <- dataset_ssr %>%
  # Group by day, year and sampling site
  group_by(Fortnight, Year, Code_point_Libelle) %>%
  # Ditch other variables
  select(-c(time, Date, Week, Day)) %>%
  # Take the mean of other variables
  summarise(across(c(latitude, longitude, ssr), 
                   ~ median(.)))


## Taking a look ####

# Let's take a look at that data  

# We will plot the daily mean of ssr over the whole study period,
# for each sampling point

dataset_daily_mean_plot <- dataset_daily_mean %>%
  # Adding back the date variable
  # To do this, we will exploit the fact that adding a number to an object of type
  # 'Date' just adds the appropriate number of days.
  # Create a month and day character variable, that is always '01-01'
  mutate(MonthDay = '01-01') %>%
  # Then unite it with Year to create something that will resemble a ymd date
  mutate(CharYear = as.character(Year)) %>%
  unite(Date, CharYear, MonthDay, sep = '-') %>%
  # And make it a 'Date' object
  mutate(Date = ymd(Date)) %>%
  # Now, add Day-1 to each
  mutate(Date = Date + (as.numeric(Day)-1)) %>%
  # Code_point_Libelle as factor
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis'))

# Same for ssr dataset
dataset_ssr_daily_mean_plot <- dataset_ssr_daily_mean %>%
  # Adding back the date variable
  # To do this, we will exploit the fact that adding a number to an object of type
  # 'Date' just adds the appropriate number of days.
  # Create a month and day character variable, that is always '01-01'
  mutate(MonthDay = '01-01') %>%
  # Then unite it with Year to create something that will resemble a ymd date
  mutate(CharYear = as.character(Year)) %>%
  unite(Date, CharYear, MonthDay, sep = '-') %>%
  # And make it a 'Date' object
  mutate(Date = ymd(Date)) %>%
  # Now, add Day-1 to each
  mutate(Date = Date + (as.numeric(Day)-1)) %>%
  # Code_point_Libelle as factor
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis'))

# Color palette
pheno_palette12 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646')

# Here we plot the tcc, but we can look at any of the 4 variables (ssr, tcc, 
# u10 or v10) with the same plot
ggplot(dataset_daily_mean_plot, 
       aes(x = Date, y = tcc, color = Code_point_Libelle)) +
  geom_point(size = .5, alpha = .8) +
  # aesthetics
  facet_wrap(facets = 'Code_point_Libelle') +
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  labs(y = 'Total Cloud Cover (tcc)') +
  theme_classic()

ggplot(dataset_ssr_daily_mean_plot, 
       aes(x = Day, y = ssr, color = Code_point_Libelle)) +
  geom_point(size = .5, alpha = .8) +
  # aesthetics
  facet_wrap(facets = 'Code_point_Libelle') +
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  labs(y = 'Surface Solar Radiation (ssr)') +
  theme_classic()

# Nice. We can be quite happy with that.

### Saving the new aggregated datasets ####

# Now we can save the daily/weekly/fortnightly aggregated datasets so we can
# work with them later.

# ## Daily datasets
# # Mean
# write.csv2(dataset_daily_mean, file = 'Data/Models/ERA5/Outputs/era5_dataset_daily_mean.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')
# # Median
# write.csv2(dataset_daily_median, file = 'Data/Models/ERA5/Outputs/era5_dataset_daily_median.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')
#
# ## Weekly datasets
# # Mean
# write.csv2(dataset_weekly_mean, file = 'Data/Models/ERA5/Outputs/era5_dataset_weekly_mean.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')
# # Median
# write.csv2(dataset_weekly_median, file = 'Data/Models/ERA5/Outputs/era5_dataset_weekly_median.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')
#
# ## Fortnightly datasets
# # Mean
# write.csv2(dataset_fortnight_mean, file = 'Data/Models/ERA5/Outputs/era5_dataset_fortnight_mean.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')
# # Median
# write.csv2(dataset_fortnight_median, file = 'Data/Models/ERA5/Outputs/era5_dataset_fortnight_median.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')

# Same with the ssr datasets
# ## Daily datasets
# # Mean
# write.csv2(dataset_ssr_daily_mean, file = 'Data/Models/ERA5/Outputs/era5_dataset_ssr_daily_mean.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')
# # Median
# write.csv2(dataset_ssr_daily_median, file = 'Data/Models/ERA5/Outputs/era5_dataset_ssr_daily_median.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')
#
# ## Weekly datasets
# # Mean
# write.csv2(dataset_ssr_weekly_mean, file = 'Data/Models/ERA5/Outputs/era5_dataset_ssr_weekly_mean.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')
# # Median
# write.csv2(dataset_ssr_weekly_median, file = 'Data/Models/ERA5/Outputs/era5_dataset_ssr_weekly_median.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')
#
# ## Fortnightly datasets
# # Mean
# write.csv2(dataset_ssr_fortnight_mean, file = 'Data/Models/ERA5/Outputs/era5_dataset_ssr_fortnight_mean.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')
# # Median
# write.csv2(dataset_ssr_fortnight_median, file = 'Data/Models/ERA5/Outputs/era5_dataset_ssr_fortnight_median.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')
