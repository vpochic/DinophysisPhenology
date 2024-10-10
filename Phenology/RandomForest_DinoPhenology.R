###### Random forest model for Dinophysis phenology ###
## V. POCHIC
# 2024-10-03

# The goal here is to apply a random forest model to our data on Dinophysis
# accumulation/loss rates, in order to identify the variables that are best
# predictors for the phenology

# Based on the tutorial on machine learning by Bede Davies
# https://bedeffinianrowedavies.com/statisticsmltutorials/intromachinelearning

# Careful! We need data files computed with the following scripts:
# GAM Dino unified more sites - mgcv.R (for the GAM models)
# Deriv_GAMs_DinoPhenology_REPHY.R (for the derivatives = acc. rates)
# Dino_phenology_heatmaps.R (for the environmental variables)

# We'll see in the future if we add biotic variables (alloxanthin and 
# Mesodinium cell counts)

### Packages ####
library(tidymodels)
library(tidyverse)
library(vip)
library(ranger)

### Import data ####
## Part 1: Dinophysis sampling data ####
# That's the file Season_Dino.csv with all sampling events in all sites for
# 2007-2022
Season_Dino <- read.csv2('Season_Dino.csv', header = TRUE, 
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

## Part 2: gam derivative data ####

# We import the dataframe of the derivative of the GAM.
# From the script Deriv_GAMs_DinoPhenology...
gam_Dino.d <- read.csv2('Derivatives_GAM_Dino_12sites_20241003.csv', 
                        header = TRUE, fileEncoding = 'ISO-8859-1')


## Part 3: environmental data ####
### Hydrological parameters form the REPHY dataset
Table_hydro <- read.csv2('Table1_hydro_select.csv', header = TRUE, 
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

### Stratification index from the GAMAR model
Table_stratif <- read.csv2('Stratif_index_GAMAR_12sites.csv',
                           header = TRUE, fileEncoding = 'ISO-8859-1') %>%
  # Getting it in date format
  mutate(Target_Date = ymd_hms(Target_Date))

# The stratif file is much heavier because we have ALL THE VALUES FOR ALL THE
# DATES ON ALL SITES
# By joining our files and ditching what's superfluous, we should manage


#### Joining datasets ####

# Let's join Season_Dino and stratif to see if this goes according to plan
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

# And now the tricky part: we project the derivative on every date based on the
# calendar day.
## Note: it would be even better to integrate the random effect of the year in
# the future

# First, we remove duplicates in the GAM derivatives dataset 
# (because Year isn't included)
gam_Dino.d_distinct <- distinct(gam_Dino.d)

Table_data_RF <- left_join(Table_data_RF, gam_Dino.d_distinct, 
                           by = c('Code_point_Libelle', 'Day'), suffix = c('',''))

# Let's see how much data we have in each site
ggplot(Table_data_RF) +
  geom_histogram(aes(x = Code_point_Libelle, fill = Code_point_Libelle),
                 stat = 'count') +
  theme_classic()
# That's a nice rainbow. It seems we have between 250 (Antifer, At so) and 400
# (Auger, Le Cornard) data points. This seems ok.

#### Applying the random forest model ####

### Preparing model data ####

# For the data table that we will use for the random forest, we only keep the
# response variable (.derivative) and the predictor variables (TEMP, SALI,
# CHLOROA, X.14.Day_Average_SI, Stratification_Index). Let's also try to keep
# Code_point_Libelle to see if unmonitored local factors play a big role
RF_data <- Table_data_RF %>%
select(.derivative, TEMP, SALI, CHLOROA, X14.Day_Average_SI, 
         Stratification_Index, Code_point_Libelle) %>%
  # We drop any NA value in these variables
  drop_na()

# Did this modify the balance between sites?
ggplot(RF_data) +
  geom_histogram(aes(x = Code_point_Libelle, fill = Code_point_Libelle),
                 stat = 'count') +
  theme_classic()
# Not really. just a little drop in Arcachon

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
# and we convert our categorical variable ('Code_point_Libelle') into a dummy 
# variable
RF_recipe <- recipe(.derivative ~ ., data = RF_train) %>%
  step_dummy(Code_point_Libelle)

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
  # and we convert our categorical variable ('Code_point_Libelle') into a dummy 
  # variable
  RF_recipe <- recipe(.derivative ~ ., data = RF_train) %>%
    step_dummy(Code_point_Libelle)
  
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
# write.csv2(dataplot_vip, 'Randomforest_vip_data_20241010.csv', row.names = FALSE)

## Variable importance plot ####

# We pivot longer the values of variable importance to get tidy data

dataplot_vip_tidy <- pivot_longer(dataplot_vip, cols = -c('Variable'), 
                                  names_to ='Seed') %>%
  group_by(Variable)

# To make a nice plot, we need to order the variables by order of importance
# So we calculate the median for each
dataplot_vip_median <- dataplot_vip_tidy %>%
  summarise(value.median = median(value), .groups = 'keep') %>%
  arrange(desc(value.median))

# Nice. Now reorder the factor
dataplot_vip_tidy <- dataplot_vip_tidy %>%
  mutate(Variable = as_factor(Variable)) %>%
  # relevel the factor in descending order of variable importance
  mutate(Variable = fct_relevel(Variable, 'TEMP', 'CHLOROA', 'SALI', 
                                "X14.Day_Average_SI", 
                                "Code_point_Libelle_Ouest.Loscolo",
                                "Stratification_Index",
                                "Code_point_Libelle_Arcachon...Bouée.7",
                                "Code_point_Libelle_Cabourg",
                                "Code_point_Libelle_Auger",
                                "Code_point_Libelle_Men.er.Roue",
                                "Code_point_Libelle_Le.Cornard",
                                "Code_point_Libelle_Teychan.bis",
                                "Code_point_Libelle_At.so",
                                "Code_point_Libelle_les.Hébihens",
                                "Code_point_Libelle_Point.1.Boulogne",
                                "Code_point_Libelle_Loguivy"))

# Plot
ggplot(data = dataplot_vip_tidy, aes(x = Variable, y = value)) +
  # A dot plot
  geom_point(position = position_jitter(width = .1),
             color = 'dodgerblue3', alpha = .4) +
  # With a boxplot as overlay
  geom_boxplot(fill = 'transparent', color = 'dodgerblue4', outliers = FALSE) +
  # Custom label for y axis
  labs(x = NULL, y = 'Variable importance') +
  # Flip x and y axes
  coord_flip() +
  # Reverse x axis so most important variable appears at the top
  scale_x_discrete(limits = rev) +
  theme_classic()

# Save the plot
ggsave('VIP_RandomForest_20241010.tiff', width = 225, height = 170, units = 'mm',
       compression = 'lzw', dpi = 300)
