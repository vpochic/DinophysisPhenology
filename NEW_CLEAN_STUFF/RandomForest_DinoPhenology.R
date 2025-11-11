###### Random forest model for Dinophysis phenology ###
## V. POCHIC
# 2025-11-05

# The goal here is to apply a random forest model to our data on Dinophysis
# accumulation/loss rates, in order to identify the variables that are best
# predictors for the phenology

# Based on the tutorial on machine learning by Bede Davies
# https://bedeffinianrowedavies.com/statisticsmltutorials/intromachinelearning

# Careful! We need data files computed with the following scripts:
# GAM_Dino_all_sites-mgcv.R (for the GAM models)
# Deriv_GAMs_DinoPhenology_REPHY.R (for the derivatives = acc. rates)
# Dino_phenology_envdata.R (for the environmental variables)

### Packages ####
library(tidymodels)
library(tidyverse)
library(vip)
library(ranger)
library(DescTools)
library(viridis)
library(ggpointdensity)

### Import data ####
## Part 1: Dinophysis sampling data ####
# That's the file Season_Dino.csv with all sampling events in all sites for
# 2007-2022
Season_Dino <- read.csv2('Data/REPHY_outputs/Season_Dino_20250604.csv', header = TRUE, 
                         fileEncoding = 'ISO-8859-1')


Season_Dino_12sites <- filter(Season_Dino, Code_point_Libelle %in% 
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
                                  'Arcachon - Bouée 7', 'Teychan bis')) %>%
  # Get rid of the one case in which Heure is ""
  # (Antifer, 2020/05/27, 0 Dinophysis)
  filter(Heure != "") %>%
  # Get a date-time variable by pasting together date and time
  mutate(DateTime = paste0(Date, ' ', Heure)) %>%
  # Get the Target date in the right format for corresponding with the stratif
  # dataset
  mutate(Target_Date = ymd_hms(DateTime)) %>%
  # And get the Date variable in the right format
  mutate(Date = ymd(Date))

## Part 2: Mesodinium sampling data ####
# That's the file Season_Meso.csv with all sampling events in all sites for
# 2007-2022
Season_Meso <- read.csv2('Data/REPHY_outputs/Season_Meso_20250604.csv', 
                         header = TRUE, fileEncoding = 'ISO-8859-1')


Season_Meso_12sites <- filter(Season_Meso, Code_point_Libelle %in% 
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
                                  'Arcachon - Bouée 7', 'Teychan bis')) %>%
  # Get rid of the potential cases in which Heure is ""
  filter(Heure != "") %>%
  # Get a date-time variable by pasting together date and time
  mutate(DateTime = paste0(Date, ' ', Heure)) %>%
  # Get the Target date in the right format for corresponding with the stratif
  # dataset
  mutate(Target_Date = ymd_hms(DateTime)) %>%
  # And get the Date variable in the right format
  mutate(Date = ymd(Date))

# Same number of events as the Dinophysis tibble: seems ok

## Part 3: Pigment data ####
# We want to get the HPLC pigment data to look at the alloxanthin

Table_pigments <- read.csv2('Data/REPHY_outputs/Table_pigments_2007-2022.csv', 
                            header = TRUE, fileEncoding = 'ISO-8859-1')

# Only keep the 4 sites we will use for our second RF model
Table_pigments_select <- filter(Table_pigments, Code_point_Libelle %in% 
                                c(# Baie de Seine
                                  'Antifer ponton pétrolier', 'Cabourg',
                                  # Bretagne Sud
                                  'Men er Roue', 'Ouest Loscolo')) %>%
  # And get the Date variable in the right format
  mutate(Date = ymd(Date))

# This dataset is a bit different from the others.

## Part 4: gam derivative data ####

# We import the dataframe of the derivative of the GAM.
# From the script Deriv_GAMs_DinoPhenology...
gam_Dino.d <- read.csv2('Data/GAM_outputs/Derivatives_GAM_Dino_12sites_20241203.csv', 
                        header = TRUE, fileEncoding = 'ISO-8859-1')

# Or with the multiyear data
gam_Dino.d_multiyear <- read.csv2('Data/GAM_outputs/gam_Dino_multiyear_deriv_20241204.csv', 
                                  header = TRUE, fileEncoding = 'ISO-8859-1') %>%
  # Put Date in date format
  mutate(Date = ymd(Date))

## Part 5: environmental data ####
### Hydrological parameters form the REPHY dataset
Table_hydro <- read.csv2('Data/REPHY_outputs/Table1_hydro_select.csv', header = TRUE, 
                         fileEncoding = "ISO-8859-1")

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
  mutate(Day = as.numeric(yday(Date))) %>%
  # Create a 'fortnight' variable to match the sampling frequency
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2))

# Checking for NAs
# Temperature
isNA_TEMP <- filter(Table_hydro_select, is.na(TEMP) == TRUE)
# Clean
# Salinity
isNA_SALI <- filter(Table_hydro_select, is.na(SALI) == TRUE)
# Clean
# Chlorophyll a
isNA_CHLOROA <- filter(Table_hydro_select, is.na(CHLOROA) == TRUE)
# Not clean. 733 obs with CHLOROA = NA, from 2007 to 2022, in several sites

## Part 6: Model data ####

### Stratification index from the GAMAR model ##
Table_stratif <- read.csv2('Data/Models/GAMAR/Outputs/Stratif_index_GAMAR_12sites.csv',
                           header = TRUE, fileEncoding = 'ISO-8859-1') %>%
  # Getting it in date format
  mutate(Target_Date = ymd_hms(Target_Date))

# The stratif file is much heavier because we have ALL THE VALUES FOR ALL THE
# DATES ON ALL SITES
# By joining our files and ditching what's superfluous, we should manage

### Atmospheric data from the ERA5 model ##

# These data are in 3 period formats (daily, weekly, fortnightly) and in 2 stats
# format (median over the period or mean over the period)

## Daily wind and tcc data
# Mean
Table_era5_daily_mean <- read.csv2('Data/Models/ERA5/Outputs/era5_dataset_daily_mean.csv', 
                                   header = TRUE, fileEncoding = 'ISO-8859-1')

# Surface Solar Radiation
Table_ssr_daily_mean <- read.csv2('Data/Models/ERA5/Outputs/era5_dataset_ssr_daily_mean.csv', 
                                  header = TRUE, fileEncoding = 'ISO-8859-1')

#### Joining datasets ####

## First, let's join Season_Dino and stratif to see if this goes according to plan
# (1 sampling event = 1 row in the stratif dataset)
Table_data_RF <- left_join(Season_Dino_12sites, Table_stratif, 
                           by = c('Target_Date', 'Code_point_Libelle'),
                           suffix = c('',''))

# It works. It just works.

# Now let's join the hydrology table

Table_data_RF <- left_join(Table_data_RF, Table_hydro_select, 
                           by = c('Code_point_Libelle', 'Date'),
                           suffix = c('',''))

# Fantastic

## Third, the ERA5 data (tcc, u10, v10, only daily mean for now)
# events are associated by day and year (and sampling site of course)
Table_data_RF <- left_join(Table_data_RF, Table_era5_daily_mean, 
                           by = c('Year', 'Day', 'Code_point_Libelle'),
                           suffix = c('',''))

## Last, the ERA5 ssr data
# events are associated by day and year (and sampling site of course)
Table_data_RF <- left_join(Table_data_RF, Table_ssr_daily_mean, 
                           by = c('Year', 'Day', 'Code_point_Libelle'),
                           suffix = c('',''))



# Let's write this table to use it later in other scripts
# write.csv2(Table_data_RF, 'Data/Models/Table_env_RF_daily.csv',
#            fileEncoding = 'ISO-8859-1', row.names = FALSE)

# Good.

# And now the tricky part: we project the derivative on every date based on the
# calendar day.
## Note: it would be even better to integrate the random effect of the year in
# the future (done)

# First, we remove duplicates in the GAM derivatives dataset (not the one with 
# every year) because Year isn't included
gam_Dino.d_distinct <- distinct(gam_Dino.d)

Table_data_RF_uniyear <- left_join(Table_data_RF, gam_Dino.d_distinct, 
                           by = c('Code_point_Libelle', 'Day'), suffix = c('',''))

# For the one with 
Table_data_RF_multiyear <- left_join(Table_data_RF, gam_Dino.d_multiyear, 
                                     by = c('Code_point_Libelle', 'Date'), 
                                     suffix = c('','')) %>%
  # change the name of the derivative for consistency
  mutate(.derivative = deriv)

# Let's see how much data we have in each site
ggplot(Table_data_RF_multiyear) +
  geom_histogram(aes(x = Code_point_Libelle, fill = Code_point_Libelle),
                 stat = 'count') +
  theme_classic()
# That's a nice rainbow. It seems we have between 250 (Antifer, At so) and 400
# (Auger, Le Cornard) data points. This seems ok.

### WARNING ### 2025/03/17
# There seems to be a big problem with TURB values (only 505 are not NA). It
# probably has to do with FNU to NTU conversion that was not done in the data
# curating scripts. We should correct that later.

#### First random forest model: 10 sites with fewer parameters ####
### 1.1. Daily data ####
# For this first step, we will consider data at the daily level
# In a second step, we may integrate everything over a period of 2 weeks

### Preparing model data ####

# We will get rid of 1 sampling point (Teychan) because the stratification data of the
# oceanic model is bad there (problems with shallow water column in this lagoon site)
Table_data_RF_multiyear_select <- filter(Table_data_RF_multiyear,
                                          Code_point_Libelle != 'Teychan bis')

### We will split the dataset to have a training and a validation set
## First, we will exclude 5 random years that will serve as our true validation 
# set

set.seed(123)

## Validation dataset
# We keep the data for one sampling site, all years, to use as our validation
# dataset
Validation_data <- filter(Table_data_RF_multiyear_select, 
                          Code_point_Libelle == 'Men er Roue') %>%
  # Then select parameters
  select(.derivative, TEMP, SALI, CHLOROA, Stratification_Index, 
         ssr, tcc, u10, v10) %>%
  # We drop any NA value in these variables
  drop_na()

# For the data table that we will use for the random forest, we only keep the
# response variable (.derivative) and the predictor variables (TEMP, SALI,
# CHLOROA, Stratification_Index, ssr, tcc, u10, v10). Let's also try to keep
# Code_point_Libelle to see if unmonitored local factors play a big role
RF_data <- Table_data_RF_multiyear_select %>%
  # We filter out the sampling site we will use for validation
  filter(Code_point_Libelle != 'Men er Roue') %>%
  # Then select parameters
  select(.derivative, TEMP, SALI, CHLOROA, Stratification_Index, 
         ssr, tcc, u10, v10) %>%
  # We drop any NA value in these variables
  drop_na()

# Standard procedure for random forest building
splitdata_RF <- initial_split(RF_data)

RF_train <- splitdata_RF %>%
  training()

RF_test <- splitdata_RF %>%
  testing()

# We set a seed to ensure that we get consistent results each time 
# the code is run (because the split is randomly different each time)
set.seed(123)

# 'strata = .derivative' allows us to ensure that different levels of the 
# response variable (the derivative in our case) are present in each fold 
RF_train_folds <- vfold_cv(RF_train, strata = .derivative)

RF_train_folds

# This seems ok

# We tell the recipe of the model (predict '.derivative' based on everything else
RF_recipe <- recipe(.derivative ~ ., data = RF_train)

# We create a 'juiced' dataset with the recipe applied to it, for later on
RF_juiced <- prep(RF_recipe) %>%
  juice()

# OK, next


### Model building ####

# We create a model object that we can tune
tune_spec <- rand_forest(
  mtry = tune(),
  # number of trees
  trees = 500,
  min_n = tune()) %>%
  # We set the mode as regression as we want to predict a continuous variable
  set_mode('regression') %>%
  set_engine('ranger')

# Create a workflow in which we put the recipe and the model
tune_wf <- workflow() %>%
  add_recipe(RF_recipe) %>%
  add_model(tune_spec)

# We set the same seed as defined before
set.seed(123)

# This step takes some time
tune_res <- tune_grid(
  tune_wf,
  resamples = RF_train_folds,
  grid = 10)

# Plotting the results of each fold using the RMSE
tune_res %>%
  collect_metrics() %>%
  filter(.metric == 'rmse') %>%
  select(mean, min_n, mtry) %>%
  pivot_longer(min_n:mtry,
               values_to = 'value',
               names_to = 'parameter') %>%
  ggplot(aes(x = value, y = mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = 'free_x') +
  labs(x = NULL, y = 'RMSE') +
  theme_classic()

# Use a function to select the best hyper parameters (min_n and mtry) for our
# model
best_rmse <- select_best(tune_res, metric = 'rmse')

# Final random forest
final_RF <- finalize_model(
  tune_spec,
  best_rmse)

final_RF

## Model performance ####

# With our split dataset (training vs testing), we will see how the model 
# performs on independent data

final_wf <- workflow() %>%
  add_recipe(RF_recipe) %>%
  add_model(final_RF)

final_res <- final_wf %>%
  last_fit(splitdata_RF)

final_res %>%
  collect_metrics()

final_model <- final_res %>%
  extract_workflow()

final_res %>%
  unnest(cols=.predictions) %>%
ggplot() +
  geom_point(aes(x = .derivative, y = .pred)) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  labs(title = 'Model performance (training dataset)',
       x = 'True values', y = 'Predicted values') +
  theme_classic()

# Save this plot
# ggsave('Plots/RF_models/Modelperf_RF10inst_training.tiff',
#        height = 140, width = 140, units = 'mm',
#        dpi = 300, compression = 'lzw')

# Using the model on the validation dataset
MeR_prediction <- Validation_data %>%
  mutate(predict(final_model,.))

# Calculate the rmse
rmse <- rmse(MeR_prediction, truth = .derivative, estimate = .pred)
rmse_character <- as.character(round(rmse$.estimate, digits = 3))

# Plot
ggplot(MeR_prediction) +
  # Points
geom_point(aes(x = .derivative, y = .pred), color = '#2156A1') +
  # Lines
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_smooth(aes(x = .derivative, y = .pred), method = 'lm') +
  # RMSE
  annotate(geom = 'text', x = 0.055, y = 0.033, 
           label = c(paste('RMSE: ',
                           rmse_character,
                           '')),
           color = 'black',
           size = 3.25) +
  labs(title = 'Model performance
(validation dataset: Men er Roue site)',
       x = 'True values', y = 'Predicted values') +
theme_classic()

# Save this plot
# ggsave('Plots/RF_models/Modelperf_RF10inst_valid_20251014.tiff',
#        height = 140, width = 140, units = 'mm',
#        dpi = 300, compression = 'lzw')

# I think to interprete this (the rmse) we need to know the order of magnitude 
# of our response variable

ggplot(data = Table_data_RF_multiyear_select) +
  geom_histogram(aes(x = .derivative), binwidth = .01) +
  # A red line representing the RMSE value
  geom_vline(aes(xintercept = 0.022), color = 'red') +
  geom_vline(aes(xintercept = -0.022), color = 'red') +
  # Text
  annotate(geom = 'text', x = 0.25, y = 500, 
           label = c(paste('RMSE: ',
                           rmse_character,
                           '')),
           color = 'red',
           size = 3.25) +
  # Labels
  labs(title = 'RMSE compared to response
variable distribution', x = 'Count', y = 'Value of response variable') +
  scale_y_continuous() +
  # Set limits to x axis to "zoom in"
  scale_x_continuous(limits = c(-0.5,0.5)) +
  theme_classic()

# Save this plot
# ggsave('Plots/RF_models/RMSE_vs_values_RF10inst_20251014.tiff',
#        height = 140, width = 140, units = 'mm',
#        dpi = 300, compression = 'lzw')

### Investigating the most important features of the model ####

# keep the same seed
set.seed(123)

final_RF %>%
  set_engine('ranger', importance = 'permutation') %>%
  fit(.derivative~., data = RF_juiced) %>%
  vip(geom = 'point') +
  theme_classic()

# We got a nice plot of variable importance!
# But we want one that integrates the random nature of the model
# (i.e., we want error bars)

# The idea here is to run the model 20 times on 20 different seeds, 
# and to collect the variable importance data each time, so we can plot them 
# at the end.

# The for loop will include exactly the same code as before, the only thing that 
# will change every time is the seed

# We also need to create a dataframe to store the results

# Get the variables from the model
vip_info <- final_RF %>%
  set_engine('ranger', importance = 'permutation') %>%
  fit(.derivative~., data = RF_juiced)

# Convert it as tibble to get a vector of names
vip_tibble <- as_tibble(vip_info$fit$variable.importance, rownames = 'Variable')

dataplot_vip <- as_tibble_col(vip_tibble$Variable, column_name = 'Variable')
# Nice!
rm(vip_info)
rm(vip_tibble)

# A big ugly for loop to run 20 iterations of the model
# (this step takes some time, obviously)

for (i in 1:20) {
  set.seed(i)
  
  # 'strata = .derivative' allows us to ensure that different levels of the 
  # response variable (the derivative in our case) are present in each fold 
  RF_train_folds <- vfold_cv(RF_train, strata = .derivative)
  
  RF_train_folds
  
  # This seems ok
  
  # We tell the recipe of the model (predict '.derivative' based on everything else,
  # and we DON'T convert our categorical variable ('Code_point_Libelle') into a dummy 
  # variable
  RF_recipe <- recipe(.derivative ~ ., data = RF_train) # %>%
  #   step_dummy(Code_point_Libelle)
  
  # We create a 'juiced' dataset with the recipe applied to it, for later on
  RF_juiced <- prep(RF_recipe) %>%
    juice()
  
  # OK, next
  
  # We create a model object that we can tune
  tune_spec <- rand_forest(
    mtry = tune(),
    # number of trees
    trees = 500,
    min_n = tune()) %>%
    # We set the mode as regression as we want to predict a continuous variable
    set_mode('regression') %>%
    set_engine('ranger')
  
  # Create a workflow in which we put the recipe and the model
  tune_wf <- workflow() %>%
    add_recipe(RF_recipe) %>%
    add_model(tune_spec)
  
  # We set the same seed as defined before
  set.seed(i)
  
  # This step takes some time
  tune_res <- tune_grid(
    tune_wf,
    resamples = RF_train_folds,
    grid = 10)
  
  # Use a function to select the best hyper parameters (min_n and mtry) for our
  # model
  best_rmse <- select_best(tune_res, metric = 'rmse')
  
  # Final random forest
  final_RF <- finalize_model(
    tune_spec,
    best_rmse)

  # Re-run the model on same seed and store the info in a "passerelle" object
  set.seed(i)
  
  passerelle_RF <- final_RF %>%
    set_engine('ranger', importance = 'permutation') %>%
    fit(.derivative~., data = RF_juiced)

  # Append the variable importance at the end of the tibble we constructed 
  # before
  # 1st, build a tibble with the var importance and the variables
  added_tibble <- as_tibble(passerelle_RF$fit$variable.importance, 
                                rownames = 'Variable')
  # 2nd, left_join with our existing dataframe by 'Variable' with appropriate
  # suffix
  dataplot_vip <- left_join(x = dataplot_vip, y = added_tibble, 
                            by = c('Variable'), 
                            suffix = c('',as.character(i)))
  
}

# Save the result!
write.csv2(dataplot_vip, 'Data/RF_outputs/Randomforest10_vip_data_20251014.csv',
           row.names = FALSE, fileEncoding = 'ISO-8859-1')

## Variable importance plot ####

# If necessary, import data previously saved
# dataplot_vip <- read.csv2('Data/RF_outputs/Randomforest10_vip_data_20251014.csv', header = TRUE,
#                           fileEncoding = 'ISO-8859-1')

# We pivot longer the values of variable importance to get tidy data

dataplot_vip_tidy <- pivot_longer(dataplot_vip, cols = -c('Variable'), 
                                  names_to ='Seed') %>%
  group_by(Variable)

# To make a nice plot, we need to order the variables by order of importance
# So we calculate the median for each
dataplot_vip_median <- dataplot_vip_tidy %>%
  summarise(value.median = median(value),
            value.sd = sd(value), .groups = 'keep') %>%
  arrange(desc(value.median))

# Nice. Now reorder the factor
dataplot_vip_tidy <- dataplot_vip_tidy %>%
  mutate(Variable = as_factor(Variable)) %>%
  # relevel the factor in descending order of variable importance
  mutate(Variable = fct_relevel(Variable, 'TEMP', 'ssr', 'v10', 'SALI', 
                                'Stratification_Index', 'tcc', 'CHLOROA', 'u10'))

# Color palette
palette_bretagne8 <- c('#FBA823', 'red3', '#FBB665', '#11203E', '#377185',
                       '#B47E24', '#1F3700', '#FBB646')

# Plot
ggplot(data = dataplot_vip_tidy, aes(x = Variable, y = value)) +
  # A dot plot of medians
  geom_point(data = dataplot_vip_median,
             aes(x = Variable, y = value.median),
             size = 4.5, shape = 21, color = 'black', fill = 'grey40') +
  # With a dot plot of all simulations as overlay
  geom_point(position = position_jitter(width = .1),
             alpha = .4, size = 1.5, color = 'grey10') +
  # Custom label for y axis
  labs(x = NULL, y = 'Variable importance') +
  # Flip x and y axes
  coord_flip() +
  # Reverse x axis so most important variable appears at the top
  scale_x_discrete(limits = rev) +
  theme_classic()

# Save the plot
# ggsave('Plots/RF_models/VIP_RandomForest_10sites_20251014.tiff', width = 164, height = 130, units = 'mm',
#        compression = 'lzw', dpi = 300)

# According to Dr. Bede Davies, it's better to represent it as mean +-std error
# So we'll do just that and calculate the mean for each
dataplot_vip_mean <- dataplot_vip_tidy %>%
  summarise(value.mean = mean(value), .groups = 'keep') %>%
  arrange(desc(value.mean))

# Same order as the median, that's convenient.
# Calculate the standard deviation
dataplot_vip_sd <- dataplot_vip_tidy %>%
  summarise(value.sd = sd(value), .groups = 'keep')

# Group the 2 dataframes
dataplot_vip_stats <- left_join(dataplot_vip_mean, dataplot_vip_sd, 
                                by = c('Variable')) %>%
  # Nice. Now reorder the factor
  mutate(Variable = as_factor(Variable)) %>%
  # relevel the factor in descending order of variable importance
  mutate(Variable = fct_relevel(Variable, 'TEMP', 'ssr', 'v10', 'SALI', 
                                'Stratification_Index', 'tcc', 'CHLOROA', 'u10'))

# Same for the datset with the individual runs
dataplot_vip_tidy <- dataplot_vip_tidy %>%
  mutate(Variable = as_factor(Variable)) %>%
  # relevel the factor in descending order of variable importance
  mutate(Variable = fct_relevel(Variable, 'TEMP', 'ssr', 'v10', 'SALI', 
                                'Stratification_Index', 'tcc', 'CHLOROA', 'u10'))

# Color palette
palette_bretagne8 <- c('#FBA823', 'red3', '#FBB665', '#11203E', '#377185',
                       '#B47E24', '#1F3700', '#FBB646')

# Plot with mean +- sd
ggplot(data = dataplot_vip_stats, aes(x = Variable, y = value.mean)) +
  # all data points for 20 model runs
  geom_point(data = dataplot_vip_tidy, aes(x = Variable, y = value),
             position = position_jitter(width = .1),
             alpha = .4, size = 1.5, color = 'grey10') +
  
  # error bars
  geom_errorbar(aes(x = Variable, ymin = value.mean - value.sd,
                    ymax = value.mean + value.sd), linewidth = .5,
                color = 'grey10', width = .2) +
  
  # points for means
  geom_point(size = 4, color = 'black', fill = 'grey45',
             shape = 21, stroke = .5) +
  
  # Custom label for y axis
  labs(x = NULL, y = 'Variable importance score') +
  # Flip x and y axes
  coord_flip() +
  # Reverse x axis so most important variable appears at the top
  scale_x_discrete(limits = rev) +
  theme_classic()

# Save the plot
# ggsave('Plots/RF_models/VIP_RandomForest10_meansd_20251014.tiff', width = 100,
# height = 120, units = 'mm', compression = 'lzw', dpi = 300)

#### Second random forest model: 4 sites with more parameters ####

### Preparing model data ####

### Data for the second random forest model (4 sites) ###

Table_data_RF2 <- left_join(Table_data_RF_multiyear, Season_Meso_12sites,
                            by = c('Code_point_Libelle', 'ID.interne.passage'), 
                            suffix = c('',''))%>%
  # Keep only the 4 sites for which we are confident in Mesodinium data AND we
  # have pigment data
  filter(Code_point_Libelle %in% 
           c(# Baie de Seine
             'Antifer ponton pétrolier', 'Cabourg',
             # Bretagne Sud
             'Men er Roue', 'Ouest Loscolo')) %>%
  # Only years from 2016 onwards (first HPLC pigment data)
  filter(Year >= 2016) %>%
  # Keep only events with reliable values of nutrients (no NA)
  filter(is.na(NH4) == FALSE &
           is.na(PO4) == FALSE &
           is.na(SIOH) == FALSE &
           is.na(OXYGENE) == FALSE &
           is.na(NO3.NO2) == FALSE)

# 584 events remaining. Matching pigments now: this is more of a quality check,
# as HPLC pigments values were included in the hydrology table.

Table_data_RF2 <- left_join(Table_data_RF2, Table_pigments_select,
                            by = c('Code_point_Libelle', 'ID.interne.passage'),
                            suffix = c('','')) %>%
  # Removing any NAs for alloxanthin
  filter(is.na(Allo) == FALSE) %>%
  # Relevel the sampling site as a factor
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle, 
                                          'Antifer ponton pétrolier', 'Cabourg', 
                                          'Men er Roue', 'Ouest Loscolo'))

# Some NAs removed: we only have 469 events in the table now.

# color palette
palette_4sites <- c(# Baie de Seine
  'red3', 'orangered',
  # Men er Roue, Ouest Loscolo
  '#2156A1', '#5995E3')


# How much events we have in each site?
ggplot(Table_data_RF2) +
  geom_histogram(aes(x = Code_point_Libelle, fill = Code_point_Libelle),
                 stat = 'count') +
  scale_fill_discrete(type = palette_4sites) +
  theme_classic()
# We have approx. 1/3 more events in South Brittany compared with Seine Bay.
# Not ideal but eh.

# Let's draw time series of alloxanthin to check if what we have is sufficiently
# continuous

ggplot(Table_data_RF2)+
  geom_point(aes(x=Day, y=Allo, color = Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  scale_color_discrete(type = palette_4sites, guide = 'none') +
  theme_classic() +
  labs(y="Alloxanthin (µg/L)", x="Date")

# We indeed have a hole in the series in Seine Bay.

# For the data table that we will use for the random forest, we only keep the
# response variable (.derivative) and the predictor variables (TEMP, SALI,
# CHLOROA, X.14.Day_Average_SI, ssr, tcc, u10, v10).

RF_data <- Table_data_RF2 %>%
  # We get the Men er Roue sampling site out, we'll use it as our validation
  # dataset
  filter(Code_point_Libelle != 'Men er Roue') %>%
  # Keeping only the desired variables
  select(.derivative, TEMP, SALI, CHLOROA, Stratification_Index, 
         ssr, tcc, u10, v10, Mesodinium_genus, Allo,
         NH4, PO4, SIOH, OXYGENE, NO3.NO2) %>%
  # We drop any NA value in these variables (should be done already)
  drop_na()

# Keeping aside our validation dataset with the Men er Roue sampling site
RF_valid <- Table_data_RF2 %>%
  # We get the Men er Roue sampling site out, we'll use it as our validation
  # dataset
  filter(Code_point_Libelle == 'Men er Roue') %>%
  # Keeping only the desired variables
  select(.derivative, TEMP, SALI, CHLOROA, Stratification_Index, 
         ssr, tcc, u10, v10, Mesodinium_genus, Allo,
         NH4, PO4, SIOH, OXYGENE, NO3.NO2) %>%
  # We drop any NA value in these variables (should be done already)
  drop_na()


# Same number of events (469)

### We will split the dataset to have a training and a validation set
splitdata_RF <- initial_split(RF_data)

RF_train <- splitdata_RF %>%
  training()

RF_test <- splitdata_RF %>%
  testing()

# We set a seed to ensure that we get consistent results each time 
# the code is run (because the split is randomly different each time)
set.seed(234)

# 'strata = .derivative' allows us to ensure that different levels of the 
# response variable (the derivative in our case) are present in each fold 
RF_train_folds <- vfold_cv(RF_train, strata = .derivative)

RF_train_folds

# This seems ok

# We tell the recipe of the model (predict '.derivative' based on everything else,
# and we DON'T convert our categorical variable ('Code_point_Libelle') into a dummy 
# variable
RF_recipe <- recipe(.derivative ~ ., data = RF_train) # %>%
# step_dummy(Code_point_Libelle)

# We create a 'juiced' dataset with the recipe applied to it, for later on
RF_juiced <- prep(RF_recipe) %>%
  juice()

# OK, next


### Model building ####

# We create a model object that we can tune
tune_spec <- rand_forest(
  mtry = tune(),
  # number of trees
  trees = 500,
  min_n = tune()) %>%
  # We set the mode as regression as we want to predict a continuous variable
  set_mode('regression') %>%
  set_engine('ranger')

# Create a workflow in which we put the recipe and the model
tune_wf <- workflow() %>%
  add_recipe(RF_recipe) %>%
  add_model(tune_spec)

# We set the same seed as defined before
set.seed(234)

# This step takes some time
tune_res <- tune_grid(
  tune_wf,
  resamples = RF_train_folds,
  grid = 10)

# Plotting the results of each fold using the RMSE
tune_res %>%
  collect_metrics() %>%
  filter(.metric == 'rmse') %>%
  select(mean, min_n, mtry) %>%
  pivot_longer(min_n:mtry,
               values_to = 'value',
               names_to = 'parameter') %>%
  ggplot(aes(x = value, y = mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = 'free_x') +
  labs(x = NULL, y = 'RMSE') +
  theme_classic()

# Use a function to select the best hyper parameters (min_n and mtry) for our
# model
best_rmse <- select_best(tune_res, metric = 'rmse')

# Final random forest
final_RF <- finalize_model(
  tune_spec,
  best_rmse)

final_RF

## Model performance ####

# With our split dataset (training vs testing), we will see how the model 
# performs on independent data

final_wf <- workflow() %>%
  add_recipe(RF_recipe) %>%
  add_model(final_RF)

final_res <- final_wf %>%
  last_fit(splitdata_RF)

final_res %>%
  collect_metrics()

final_model <- final_res %>%
  extract_workflow()

final_res %>%
  unnest(cols=.predictions) %>%
  ggplot() +
  geom_point(aes(x = .derivative, y = .pred)) +
  theme_classic()

final_res %>%
  unnest(cols=.predictions) %>%
  ggplot() +
  geom_point(aes(x = .derivative, y = .pred)) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  labs(title = 'Model performance (training dataset)',
       x = 'True values', y = 'Predicted values') +
  theme_classic()

# Save this plot
# ggsave('Plots/RF_models/Modelperf_RF10inst_training.tiff',
#        height = 140, width = 140, units = 'mm',
#        dpi = 300, compression = 'lzw')

# Using the model on the validation dataset
MeR_prediction <- RF_valid %>%
  mutate(predict(final_model,.))

# Calculate the rmse
rmse <- rmse(MeR_prediction, truth = .derivative, estimate = .pred)
rmse_character <- as.character(round(rmse$.estimate, digits = 3))

# Plot
ggplot(MeR_prediction) +
  geom_point(aes(x = .derivative, y = .pred)) +
  annotate(geom = 'text', x = 0.04, y = -0.1, 
           label = c(paste('RMSE: ',
                           rmse_character,
                           '')),
           color = 'black',
           size = 3.25) +
  geom_smooth(aes(x = .derivative, y = .pred), method = 'lm') +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis() +
  theme_classic()

# I think to interprete this (the rmse) we need to know the order of magnitude 
# of our response variable

ggplot(data = Table_data_RF_multiyear_select) +
  geom_histogram(aes(x = .derivative), binwidth = .02) +
  # A red line representing the RMSE value
  geom_vline(aes(xintercept = 0.129), color = 'red') +
  geom_vline(aes(xintercept = -0.129), color = 'red') +
  scale_y_continuous() +
  theme_classic()

final_res %>%
  unnest(cols=.predictions) %>%
  ggplot() +
  geom_point(aes(x = .derivative, y = .pred)) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  labs(title = 'Model performance (training dataset)',
       x = 'True values', y = 'Predicted values') +
  theme_classic()

# Save this plot
# ggsave('Plots/RF_models/Modelperf_RF4inst_training.tiff',
#        height = 140, width = 140, units = 'mm',
#        dpi = 300, compression = 'lzw')

# Using the model on the validation dataset
MeR_prediction <- RF_valid %>%
  mutate(predict(final_model,.))

# Calculate the rmse
rmse <- rmse(MeR_prediction, truth = .derivative, estimate = .pred)
rmse_character <- as.character(round(rmse$.estimate, digits = 3))

# Plot
ggplot(MeR_prediction) +
  # Points
  geom_point(aes(x = .derivative, y = .pred), color = '#2156A1') +
  # Lines
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_smooth(aes(x = .derivative, y = .pred), method = 'lm') +
  # RMSE
  annotate(geom = 'text', x = 0.04, y = -0.09, 
           label = c(paste('RMSE: ',
                           rmse_character,
                           '')),
           color = 'black',
           size = 3.25) +
  labs(title = 'Model performance
(validation dataset: Men er Roue site)',
       x = 'True values', y = 'Predicted values') +
  theme_classic()

# Save this plot
# ggsave('Plots/RF_models/Modelperf_RF4inst_valid_20251105.tiff',
#        height = 140, width = 140, units = 'mm',
#        dpi = 300, compression = 'lzw')

# I think to interprete this (the rmse) we need to know the order of magnitude 
# of our response variable

ggplot(data = Table_data_RF2) +
  geom_histogram(aes(x = .derivative), binwidth = .01) +
  # A red line representing the RMSE value
  geom_vline(aes(xintercept = 0.072), color = 'red') +
  geom_vline(aes(xintercept = -0.072), color = 'red') +
  # Text
  annotate(geom = 'text', x = 0.5, y = 100, 
           label = c(paste('RMSE: ',
                           rmse_character,
                           '')),
           color = 'red',
           size = 3.25) +
  # Labels
  labs(title = 'RMSE compared to response
variable distribution', x = 'Value of response variable', y = 'Count') +
  scale_y_continuous() +
  # Set limits to x axis to "zoom in"
  scale_x_continuous(limits = c(-1.2,1.2)) +
  theme_classic()

# Save this plot
# ggsave('Plots/RF_models/RMSE_vs_values_RF4inst_20251105.tiff',
#        height = 140, width = 140, units = 'mm',
#        dpi = 300, compression = 'lzw')

### Investigating the most important features of the model ####

# keep the same seed
set.seed(234)

final_RF %>%
  set_engine('ranger', importance = 'permutation') %>%
  fit(.derivative~., data = RF_juiced) %>%
  vip(geom = 'point') +
  theme_classic()

# We got a nice plot of variable importance!
# But we want one that integrates the random nature of the model
# (i.e., we want error bars)

# The idea here is to run the model 20 times on 20 different seeds, 
# and to collect the variable importance data each time, so we can plot them 
# at the end.

# The for loop will include exactly the same code as before, the only thing that 
# will change every time is the seed

# We also need to create a dataframe to store the results

# Get the variables from the model
vip_info <- final_RF %>%
  set_engine('ranger', importance = 'permutation') %>%
  fit(.derivative~., data = RF_juiced)

# Convert it as tibble to get a vector of names
vip_tibble <- as_tibble(vip_info$fit$variable.importance, rownames = 'Variable')

dataplot_vip <- as_tibble_col(vip_tibble$Variable, column_name = 'Variable')
# Nice!
rm(vip_info)
rm(vip_tibble)

# A big ugly for loop to run 20 iterations of the model
# (this step takes some time, obviously)

for (i in 1:20) {
  set.seed(i)
  
  # 'strata = .derivative' allows us to ensure that different levels of the 
  # response variable (the derivative in our case) are present in each fold 
  RF_train_folds <- vfold_cv(RF_train, strata = .derivative)
  
  RF_train_folds
  
  # This seems ok
  
  # We tell the recipe of the model (predict '.derivative' based on everything else,
  # and we DON'T convert our categorical variable ('Code_point_Libelle') into a dummy 
  # variable
  RF_recipe <- recipe(.derivative ~ ., data = RF_train) # %>%
  #   step_dummy(Code_point_Libelle)
  
  # We create a 'juiced' dataset with the recipe applied to it, for later on
  RF_juiced <- prep(RF_recipe) %>%
    juice()
  
  # OK, next
  
  # We create a model object that we can tune
  tune_spec <- rand_forest(
    mtry = tune(),
    # number of trees
    trees = 500,
    min_n = tune()) %>%
    # We set the mode as regression as we want to predict a continuous variable
    set_mode('regression') %>%
    set_engine('ranger')
  
  # Create a workflow in which we put the recipe and the model
  tune_wf <- workflow() %>%
    add_recipe(RF_recipe) %>%
    add_model(tune_spec)
  
  # We set the same seed as defined before
  set.seed(i)
  
  # This step takes some time
  tune_res <- tune_grid(
    tune_wf,
    resamples = RF_train_folds,
    grid = 10)
  
  # Use a function to select the best hyper parameters (min_n and mtry) for our
  # model
  best_rmse <- select_best(tune_res, metric = 'rmse')
  
  # Final random forest
  final_RF <- finalize_model(
    tune_spec,
    best_rmse)
  
  # Re-run the model on same seed and store the info in a "passerelle" object
  set.seed(i)
  
  passerelle_RF <- final_RF %>%
    set_engine('ranger', importance = 'permutation') %>%
    fit(.derivative~., data = RF_juiced)
  
  # Append the variable importance at the end of the tibble we constructed 
  # before
  # 1st, build a tibble with the var importance and the variables
  added_tibble <- as_tibble(passerelle_RF$fit$variable.importance, 
                            rownames = 'Variable')
  # 2nd, left_join with our existing dataframe by 'Variable' with appropriate
  # suffix
  dataplot_vip <- left_join(x = dataplot_vip, y = added_tibble, 
                            by = c('Variable'), 
                            suffix = c('',as.character(i)))
  
}

# Save the result!
# write.csv2(dataplot_vip, 'Data/RF_outputs/Randomforest2_vip_data_20251105.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')

## Variable importance plot ####

# If necessary, import data previously saved
# dataplot_vip_import <- read.csv2('Randomforest_vip_data_20251105.csv', header = TRUE,
#                           fileEncoding = 'ISO-8859-1')

# We pivot longer the values of variable importance to get tidy data

dataplot_vip_tidy <- pivot_longer(dataplot_vip, cols = -c('Variable'), 
                                  names_to ='Seed') %>%
  group_by(Variable)

# To make a nice plot, we need to order the variables by order of importance
# So we calculate the median for each
dataplot_vip_median <- dataplot_vip_tidy %>%
  summarise(value.median = median(value),
            value.sd = sd(value), .groups = 'keep') %>%
  arrange(desc(value.median))

# Nice. Now reorder the factor
dataplot_vip_tidy <- dataplot_vip_tidy %>%
  mutate(Variable = as_factor(Variable)) %>%
  # relevel the factor in descending order of variable importance
  # Order with Code_point_Libelle (data from 2025/03/18)
  # mutate(Variable = fct_relevel(Variable, 'NH4', 'ssr', 'SIOH', 'TEMP', 'PO4', 
  #                               'NO3.NO2', 'Allo', 'u10', 'Code_point_Libelle',
  #                               'OXYGENE', 'tcc', 'CHLOROA', 'SALI',
  #                               'X14.Day_Average_SI', 'v10', 'Mesodinium_genus'  
  #                               )) %>%
  # Order without Code_point_Libelle (data from 2025/06/26)
  mutate(Variable = fct_relevel(Variable, 'NO3.NO2', 'TEMP', 'NH4', 'PO4', 'ssr',
                                'SIOH', 'SALI', 'Stratification_Index', 'CHLOROA',
                                'OXYGENE', 'Allo', 'tcc', 'Mesodinium_genus',
                                'u10', 'v10'))

# Color palette
# palette_bretagne16 <- c('#FC4D6B', 'red3', '#FF6448', '#8B064B', '#FBA823', '#DEB1CC', 
#                         '#11203E', '#377185', '#1F3700', '#2A56A1', '#9DA51E',
#                         '#B47E24', '#711412', '#FBB646', '#FBB665')

# Plot
ggplot(data = dataplot_vip_tidy, aes(x = Variable, y = value, color = Variable)) +
  # A dot plot of medians
  geom_point(data = dataplot_vip_median,
             aes(x = Variable, y = value.median, color = Variable),
             size = 4.5) +
  # With a dot plot of all simulations as overlay
  geom_point(position = position_jitter(width = .1),
             alpha = .4, size = 1.5) +
  # Custom label for y axis
  labs(x = NULL, y = 'Variable importance') +
  # Flip x and y axes
  coord_flip() +
  # Reverse x axis so most important variable appears at the top
  scale_x_discrete(limits = rev) +
  # Color palette
  scale_color_discrete(type = palette_bretagne16, guide = 'none') +
  theme_classic()

# Save the plot
# ggsave('Plots/RF_models/VIP_RandomForest2_4sites_20250318.tiff', width = 164, height = 250, units = 'mm',
#        compression = 'lzw', dpi = 300)

# According to Dr. Bede Davies, it's better to represent it as mean +-std error
# So we'll do just that and calculate the mean for each
dataplot_vip_mean <- dataplot_vip_tidy %>%
  summarise(value.mean = mean(value), .groups = 'keep') %>%
  arrange(desc(value.mean))

# Almost same order as the median, that's convenient.
# Calculate the standard deviation
dataplot_vip_sd <- dataplot_vip_tidy %>%
  summarise(value.sd = sd(value), .groups = 'keep')

# Group the 2 dataframes
dataplot_vip_stats <- left_join(dataplot_vip_mean, dataplot_vip_sd, 
                                by = c('Variable')) %>%
  # Relevel the factor to get the variables in the right order
  mutate(Variable = as_factor(Variable)) %>%
  # relevel the factor in descending order of variable importance
  mutate(Variable = fct_relevel(Variable, 'ssr', 'TEMP', 'NH4', 'NO3.NO2', 'v10',  'SIOH', 'PO4',   
                                'Allo', 'OXYGENE', 'SALI', 'tcc', 'CHLOROA', 'Stratification_Index',   
                                'u10', 'Mesodinium_genus'))

# Relevel dataplot_vip_tidy variables for the dot plot
dataplot_vip_tidy <- dataplot_vip_tidy %>%
  mutate(Variable = as_factor(Variable)) %>%
  mutate(Variable = fct_relevel(Variable, 'ssr', 'TEMP', 'NH4', 'NO3.NO2', 'v10',  'SIOH', 'PO4',   
                                'Allo', 'OXYGENE', 'SALI', 'tcc', 'CHLOROA', 'Stratification_Index',   
                                'u10', 'Mesodinium_genus'))

# Color palette with new order --> Not used
# palette_bretagne16 <- c('#FC4D6B', '#FF6448', '#DEB1CC', 'red3', '#8B064B', '#FBA823',  
#                         '#11203E', '#377185', '#2A56A1', '#9DA51E', '#1F3700', 
#                         '#B47E24', '#711412', '#FBB646', '#FBB665')

# Plot with mean +- sd
ggplot(data = dataplot_vip_stats, aes(x = Variable, y = value.mean)) +
  
  # all data points for 20 model runs
  geom_point(data = dataplot_vip_tidy, aes(x = Variable, y = value),
             position = position_jitter(width = .1),
             alpha = .4, size = 1.5, color = 'grey10') +
  
  # error bars
  geom_errorbar(aes(x = Variable, ymin = value.mean - value.sd,
                    ymax = value.mean + value.sd), linewidth = .5,
                color = 'grey10', width = .2) +
  
  # points for means
  geom_point(size = 4, color = 'black', fill = 'grey45',
             shape = 21, stroke = .5) +
  
  # Color scales
  # scale_fill_discrete(type = palette_bretagne16, guide = 'none') +
  # scale_color_discrete(type = palette_bretagne16, guide = 'none') +
  # Custom label for y axis
  labs(x = NULL, y = 'Variable importance score') +
  # Flip x and y axes
  coord_flip() +
  # Reverse x axis so most important variable appears at the top
  scale_x_discrete(limits = rev) +
  theme_classic()

# Save the plot
# ggsave('Plots/RF_models/VIP_RandomForest2_meansd_20251105.tiff', width = 164,
#        height = 180, units = 'mm', compression = 'lzw', dpi = 300)
