###### Generalised Linear Models of Dinophysis change rate ###
## V. POCHIC
# 2025-06-27

# The goal here is to examine the relationship between Dinophysis change rate
# (the derivative of the phenology GAM) and several selected environmental
# parameters

# Based on the tutorial on GLMs by Bede Davies
# https://bedeffinianrowedavies.com/statisticstutorials/gaussianglms

# Careful! We need data files computed with the following scripts:
# GAM_Dino_all_sites-mgcv.R (for the GAM models)
# Deriv_GAMs_DinoPhenology_REPHY.R (for the derivatives = acc. rates)
# Dino_phenology_heatmaps.R (for the environmental variables)

### Packages ####
library(tidymodels)
library(tidyverse)
library(vip)
library(mgcv)
library(ranger)
library(DescTools)
library(viridis)
library(ggpointdensity)
library(RColorBrewer)
library(ggnewscale)
library(cmocean)

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

## Daily
# Mean
Table_era5_daily_mean <- read.csv2('Data/Models/ERA5/Outputs/era5_dataset_daily_mean.csv', 
                                   header = TRUE, fileEncoding = 'ISO-8859-1')
# Median
Table_era5_daily_median <- read.csv2('Data/Models/ERA5/Outputs/era5_dataset_daily_mean.csv', 
                                     header = TRUE, fileEncoding = 'ISO-8859-1')

## Weekly
# Mean
Table_era5_weekly_mean <- read.csv2('Data/Models/ERA5/Outputs/era5_dataset_weekly_mean.csv', 
                                    header = TRUE, fileEncoding = 'ISO-8859-1')
# Median
Table_era5_weekly_median <- read.csv2('Data/Models/ERA5/Outputs/era5_dataset_weekly_mean.csv', 
                                      header = TRUE, fileEncoding = 'ISO-8859-1')

## Fortnightly
# Mean
Table_era5_fortnight_mean <- read.csv2('Data/Models/ERA5/Outputs/era5_dataset_fortnight_mean.csv', 
                                       header = TRUE, fileEncoding = 'ISO-8859-1')
# Median
Table_era5_fortnight_median <- read.csv2('Data/Models/ERA5/Outputs/era5_dataset_fortnight_mean.csv', 
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

## Third, the ERA5 data (only daily mean for now)
# events are associated by day and year (and sampling site of course)
Table_data_RF <- left_join(Table_data_RF, Table_era5_daily_mean, 
                           by = c('Year', 'Day', 'Code_point_Libelle'),
                           suffix = c('',''))

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
  mutate(.derivative = deriv) %>%
  # Sampling site as factor and reorder
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis'))

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
#### Model building - 12 sites ####

# We will first build models based on the 12 sites at our disposal in the
# Atlantic-Channel region. Later, we will focus on only 4 sites for which we 
# have additional data. (But let's go 1 step at a time)

## Temperature ####

# The first effect we want to test is that of sea surface temperature.

# First, let's plot the data
pheno_palette12 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646')

ggplot(data = Table_data_RF_multiyear,
       aes(x = TEMP, y = .derivative, color = Code_point_Libelle)) +
  # Points
  geom_point(size = 2, alpha = .7) +
  # Add 2 "ghost points" to force the minimum y-axis range to -0.1;0.1
  geom_point(aes(x = 10, y = -.1), color = 'transparent', fill = 'transparent') +
  geom_point(aes(x = 10, y = .1), color = 'transparent', fill = 'transparent') +
  # Color scale (sampling sites)
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  # New color scale for temperature rug plot
  new_scale_color() +
  geom_rug(data = Table_data_RF_multiyear, aes(x = TEMP, y = NULL,
                                               color = TEMP)) +
  scale_color_distiller(palette = 'RdBu', direction = -1, guide = 'none') +
  # facets
  facet_wrap(facets = 'Code_point_Libelle', scales = 'free') +
  # Labels
  labs(x = 'Sea Surface Temperature (°C)', y = 'GAM derivative',) +
  # Theme
  theme_classic()

# This seems pretty fckin difficult
# Can we try to model it with a GAM? (non-linear relationship)

# We'll do this with the 10 sites we used for the random forest (Teychan and 
# Auger excluded)

Table_GAM <- Table_data_RF_multiyear %>%
  filter(is.na(TEMP) == FALSE) %>%
  filter(Code_point_Libelle != 'Teychan bis' & 
           Code_point_Libelle != 'Auger')

gam_Dino <- gam(data = Table_GAM, 
                # Only a spline for the temperature
                # 'k = -1' allows the model to fix the 'best' number of basic
                # functions (= knots)
                formula = .derivative~s(TEMP, k = -1),
                # Introducing the weights
                # weights = unif_weight,
                # Using a Poisson distribution for count data
                family = gaussian(),
                # Restricted maximum likelihood estimation (recommended method)
                method = 'REML')

summary(gam_Dino)
gam.check(gam_Dino)

# Create a new 'empty' dataset for storing model prediction (fit)
gam_Dino_newdata <- expand_grid(TEMP=seq(min(Table_GAM$TEMP),
                                         max(Table_GAM$TEMP)))

## Get the inverse link function of the model
# With this function, we can transform the prediction we make on the link
#scale to the response scale
ilink <- gam_Dino$family$linkinv

## Predict : add fit and se.fit on the **link** scale
gam_Dino_newdata <- bind_cols(gam_Dino_newdata, 
                              setNames(as_tibble(
                                predict(gam_Dino, gam_Dino_newdata, 
                                        se.fit = TRUE, type = 'link')[1:2]),
                                c('fit_link','se_link')))

## Create the 95% confidence interval (2*standard error fit) and backtransform 
# to response variable using the inverse link function
gam_Dino_newdata <- mutate(gam_Dino_newdata,
                           fit_resp  = ilink(fit_link),
                           right_upr = ilink(fit_link + (2 * se_link)),
                           right_lwr = ilink(fit_link - (2 * se_link)))
# Check the confidence interval. It should not extend below 0 (negative counts
# are impossible)
min(gam_Dino_newdata$right_lwr) # Nice :)
min(gam_Dino_newdata$fit_resp)
max(gam_Dino_newdata$right_upr) # Nice too
max(gam_Dino_newdata$fit_resp)
# The maximum value seems really small

# Saving data to make another plot in another script
# write.csv2(gam_Dino_newdata, 'gam_Dino_multiyear_data.csv', row.names = FALSE,
#            fileEncoding = 'ISO-8859-1')

### Checking the model
ModelOutputs<-data.frame(Fitted=fitted(gam_Dino),
                         Residuals=resid(gam_Dino))

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- gam_Dino$model
# then we add the values of fitted and residuals 
# (but are they in the same order as the model? -> need to check that)
qq_data <- bind_cols(qq_data, ModelOutputs)

# Plot : verify that true data matches (more or less) model fit
ggplot(qq_data)+
  geom_point(aes(x = TEMP, y = .derivative), color = 'red') +
  geom_point(aes(x = TEMP, y = Fitted), color = 'blue') +
  theme_classic() +
  labs(y = "Dinophysis derivative", x = "TEMP")

# The data is almost symmetrical around 0, the fit is therefore near a perfect 0

# Custom qq-plot
qqplot_custom <- ggplot(qq_data) +
  stat_qq(aes(sample=Residuals), color = 'red3', alpha = .7) +
  stat_qq_line(aes(sample=Residuals), color = 'black') +
  theme_classic() +
  labs(y="Sample Quantiles",x="Theoretical Quantiles",
       title = 'qq-plot: GAM of Dinophysis derivative ~ SST')

qqplot_custom

# Save the plot
# ggsave('Plots/GAMs/Drivers/qqplot_custom_GAM_TEMP.tiff', dpi = 300, height = 164, width = 164,
#                units = 'mm', compression = 'lzw')

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data)+
  geom_point(aes(x=Fitted, y=Residuals), 
             alpha = .7, color = 'red3') +
  theme_classic() +
  labs(y="Residuals",x="Fitted Values",
       title = 'Residuals vs Fitted: 
GAM of Dinophysis derivative ~ SST')

RvFplot_custom

# Save the plot
# ggsave('Plots/GAMs/Drivers/RvFplot_custom_GAM_TEMP.tiff', dpi = 300, height = 120, width = 164,
#                units = 'mm', compression = 'lzw')

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data, aes(x = Residuals)) +
  geom_histogram(binwidth = .005, fill = 'red3') +
  theme_classic() +
  labs(x='Residuals', y = 'Count',
       title = 'Histogram of residuals: 
GAM of Dinophysis derivative ~ SST')

HistRes_custom

# Save the plot
# ggsave('Plots/GAMs/Drivers/HistRes_custom_GAM_TEMP.tiff', dpi = 300, height = 164, width = 164,
#                units = 'mm', compression = 'lzw')

# Plot the model against the data
ggplot(data = gam_Dino_newdata,
       aes(x = TEMP, y = fit_resp)) +
  geom_line(color = 'red3', linewidth = .5) +
  geom_ribbon(aes(x = TEMP,
                  ymin = right_lwr, 
                  ymax = right_upr), fill = 'red3', color = 'red3', alpha = .2) +
  geom_point(data = Table_GAM, aes(x = TEMP, y = .derivative), size = .8, 
             alpha = .5, color = 'red3') +
  labs(x = 'Sea Surface Temperature (°C)', 
       y = 'Dinophysis accumulation rate (d-1)',
       title = 'GAM fit: Dinophysis derivative ~ SST') +
  theme_classic()

# Save this plot
# ggsave('Plots/GAMs/Drivers/GAM_fit_SST.tiff', height = 135, width = 164,
#        units = 'mm', compression = 'lzw')

## Surface solar radiation ####

# Let's plot the data
# Because the ssr values are very large numbers, we will look at the log10 of
# the ssr to try and build a better model
ggplot(data = Table_data_RF_multiyear,
       aes(x = log10(ssr), y = .derivative, color = Code_point_Libelle)) +
  # Points
  geom_point(size = 2, alpha = .7) +
  # Add 2 "ghost points" to force the minimum y-axis range to -0.1;0.1
  geom_point(aes(x = 5, y = -.1), color = 'transparent', fill = 'transparent') +
  geom_point(aes(x = 5, y = .1), color = 'transparent', fill = 'transparent') +
  # Color scale (sampling sites)
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  # New color scale for ssrerature rug plot
  new_scale_color() +
  geom_rug(data = Table_data_RF_multiyear, aes(x = log10(ssr), y = NULL,
                                               color = log10(ssr))) +
  scale_color_cmocean(name = 'solar', direction = 1, guide = 'none') +
  # facets
  facet_wrap(facets = 'Code_point_Libelle', scales = 'free') +
  # Labels
  labs(x = 'Surface Solar Radiation (log scale)', y = 'GAM derivative',) +
  # Theme
  theme_classic()

# That's not clear.
# Let's check with a GAM

# We'll do this with the 10 sites we used for the random forest (Teychan and 
# Auger excluded)

Table_GAM <- Table_data_RF_multiyear %>%
  filter(is.na(ssr) == FALSE) %>%
  filter(Code_point_Libelle != 'Teychan bis' & 
           Code_point_Libelle != 'Auger') %>%
  mutate(log10ssr = log10(ssr))

gam_Dino <- gam(data = Table_GAM, 
                # Only a spline for the ssrerature
                # 'k = -1' allows the model to fix the 'best' number of basic
                # functions (= knots)
                formula = .derivative~s(log10ssr, k = -1),
                # Introducing the weights
                # weights = unif_weight,
                # Using a Poisson distribution for count data
                family = gaussian(),
                # Restricted maximum likelihood estimation (recommended method)
                method = 'REML')

summary(gam_Dino)
gam.check(gam_Dino)

# Create a new 'empty' dataset for storing model prediction (fit)
gam_Dino_newdata <- expand_grid(log10ssr=seq(min(Table_GAM$log10ssr),
                                        max(Table_GAM$log10ssr), 
                                        length.out = 300))

## Get the inverse link function of the model
# With this function, we can transform the prediction we make on the link
#scale to the response scale
ilink <- gam_Dino$family$linkinv

## Predict : add fit and se.fit on the **link** scale
gam_Dino_newdata <- bind_cols(gam_Dino_newdata, 
                              setNames(as_tibble(
                                predict(gam_Dino, gam_Dino_newdata, 
                                        se.fit = TRUE, type = 'link')[1:2]),
                                c('fit_link','se_link')))

## Create the 95% confidence interval (2*standard error fit) and backtransform 
# to response variable using the inverse link function
gam_Dino_newdata <- mutate(gam_Dino_newdata,
                           fit_resp  = ilink(fit_link),
                           right_upr = ilink(fit_link + (2 * se_link)),
                           right_lwr = ilink(fit_link - (2 * se_link)))
# Check the confidence interval. It should not extend below 0 (negative counts
# are impossible)
min(gam_Dino_newdata$right_lwr) # Nice :)
min(gam_Dino_newdata$fit_resp)
max(gam_Dino_newdata$right_upr) # Nice too
max(gam_Dino_newdata$fit_resp)
# This looks really small again

### Checking the model
ModelOutputs<-data.frame(Fitted=fitted(gam_Dino),
                         Residuals=resid(gam_Dino))

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- gam_Dino$model
# then we add the values of fitted and residuals 
# (but are they in the same order as the model? -> need to check that)
qq_data <- bind_cols(qq_data, ModelOutputs)

# Plot : verify that true data matches (more or less) model fit
ggplot(qq_data)+
  geom_point(aes(x = log10ssr, y = .derivative), color = 'red') +
  geom_point(aes(x = log10ssr, y = Fitted), color = 'blue') +
  theme_classic() +
  labs(y = "Dinophysis derivative", x = "log10(ssr)")

# The data is almost symmetrical around 0, the fit is therefore near a perfect 0

# Custom qq-plot
qqplot_custom <- ggplot(qq_data) +
  stat_qq(aes(sample=Residuals), color = '#FBA823', alpha = .7) +
  stat_qq_line(aes(sample=Residuals), color = 'black') +
  theme_classic() +
  labs(y="Sample Quantiles",x="Theoretical Quantiles",
       title = 'qq-plot: GAM of Dinophysis derivative ~ log10(SSR)')

qqplot_custom

# Save the plot
# ggsave('Plots/GAMs/Drivers/qqplot_custom_GAM_log10ssr.tiff', dpi = 300, height = 164, width = 164,
#                units = 'mm', compression = 'lzw')

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data)+
  geom_point(aes(x=Fitted, y=Residuals), 
             alpha = .7, color = '#FBA823') +
  theme_classic() +
  labs(y="Residuals",x="Fitted Values",
       title = 'Residuals vs Fitted: 
GAM of Dinophysis derivative ~ log10(SSR)')

RvFplot_custom

# Save the plot
# ggsave('Plots/GAMs/Drivers/RvFplot_custom_GAM_log10ssr.tiff', dpi = 300, height = 120, width = 164,
#                units = 'mm', compression = 'lzw')

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data, aes(x = Residuals)) +
  geom_histogram(binwidth = .005, fill = '#FBA823') +
  theme_classic() +
  labs(x='Residuals', y = 'Count',
       title = 'Histogram of residuals: 
GAM of Dinophysis derivative ~ log10(SSR)')

HistRes_custom

# Save the plot
# ggsave('Plots/GAMs/Drivers/HistRes_custom_GAM_log10ssr.tiff', dpi = 300, height = 164, width = 164,
#                units = 'mm', compression = 'lzw')

# Plot the model against the data
ggplot(data = gam_Dino_newdata,
       aes(x = log10ssr, y = fit_resp)) +
  geom_line(color = '#FBA823', linewidth = .5) +
  geom_ribbon(aes(x = log10ssr,
                  ymin = right_lwr, 
                  ymax = right_upr), fill = '#FBA823', color = '#FBA823', alpha = .2) +
  geom_point(data = Table_GAM, aes(x = log10ssr, y = .derivative), size = .8, 
             alpha = .5, color = '#FBA823') +
  labs(x = 'log10(Surface Solar Radiation)', 
       y = 'Dinophysis accumulation rate (d-1)',
       title = 'GAM fit: Dinophysis derivative ~ log10(SSR)') +
  theme_classic()

# Save this plot
# ggsave('Plots/GAMs/Drivers/GAM_fit_log10SSR.tiff', height = 135, width = 164,
#        units = 'mm', compression = 'lzw')

## Stratification Index ####

# Let's plot the data
ggplot(data = Table_data_RF_multiyear,
       aes(x = Stratification_Index, y = .derivative, color = Code_point_Libelle)) +
  # Points
  geom_point(size = 2, alpha = .7) +
  # Add 2 "ghost points" to force the minimum y-axis range to -0.1;0.1
  geom_point(aes(x = 0, y = -.1), color = 'transparent', fill = 'transparent') +
  geom_point(aes(x = 0, y = .1), color = 'transparent', fill = 'transparent') +
  # Color scale (sampling sites)
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  # New color scale for Stratification_Indexerature rug plot
  new_scale_color() +
  # geom_rug(data = Table_data_RF_multiyear, aes(x = Stratification_Index, y = NULL,
  #                                              color = Stratification_Index)) +
  scale_color_distiller(palette = 'RdBu', direction = -1, guide = 'none') +
  # facets
  facet_wrap(facets = 'Code_point_Libelle', scales = 'free_y') +
  # Labels
  labs(x = 'Stratification Index', y = 'GAM derivative',) +
  # Theme
  theme_classic()

# This one's a bit more tricky.
# Let's check with a GAM

# We'll do this with the 10 sites we used for the random forest (Teychan and 
# Auger excluded)

Table_GAM <- Table_data_RF_multiyear %>%
  filter(is.na(Stratification_Index) == FALSE) %>%
  filter(Code_point_Libelle != 'Teychan bis' & 
           Code_point_Libelle != 'Auger')

gam_Dino <- gam(data = Table_GAM, 
                # Only a spline for the Stratification_Indexerature
                # 'k = -1' allows the model to fix the 'best' number of basic
                # functions (= knots)
                formula = .derivative~s(Stratification_Index, k = -1),
                # Introducing the weights
                # weights = unif_weight,
                # Using a Poisson distribution for count data
                family = gaussian(),
                # Restricted maximum likelihood estimation (recommended method)
                method = 'REML')

summary(gam_Dino)
gam.check(gam_Dino)

# Create a new 'empty' dataset for storing model prediction (fit)
gam_Dino_newdata <- expand_grid(Stratification_Index=seq(min(Table_GAM$Stratification_Index),
                                        max(Table_GAM$Stratification_Index), length.out = 300))

## Get the inverse link function of the model
# With this function, we can transform the prediction we make on the link
#scale to the response scale
ilink <- gam_Dino$family$linkinv

## Predict : add fit and se.fit on the **link** scale
gam_Dino_newdata <- bind_cols(gam_Dino_newdata, 
                              setNames(as_tibble(
                                predict(gam_Dino, gam_Dino_newdata, 
                                        se.fit = TRUE, type = 'link')[1:2]),
                                c('fit_link','se_link')))

## Create the 95% confidence interval (2*standard error fit) and backtransform 
# to response variable using the inverse link function
gam_Dino_newdata <- mutate(gam_Dino_newdata,
                           fit_resp  = ilink(fit_link),
                           right_upr = ilink(fit_link + (2 * se_link)),
                           right_lwr = ilink(fit_link - (2 * se_link)))
# Check the confidence interval. It should not extend below 0 (negative counts
# are impossible)
min(gam_Dino_newdata$right_lwr) # Nice :)
min(gam_Dino_newdata$fit_resp)
max(gam_Dino_newdata$right_upr) # Nice too
max(gam_Dino_newdata$fit_resp)
# The maximum value seems really small

### Checking the model
ModelOutputs<-data.frame(Fitted=fitted(gam_Dino),
                         Residuals=resid(gam_Dino))

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- gam_Dino$model
# then we add the values of fitted and residuals 
# (but are they in the same order as the model? -> need to check that)
qq_data <- bind_cols(qq_data, ModelOutputs)

# Plot : verify that true data matches (more or less) model fit
ggplot(qq_data)+
  geom_point(aes(x = Stratification_Index, y = .derivative), color = 'red') +
  geom_point(aes(x = Stratification_Index, y = Fitted), color = 'blue') +
  theme_classic() +
  labs(y = "Dinophysis derivative", x = "Stratification_Index")

# The data is almost symmetrical around 0, the fit is therefore near a perfect 0

# Custom qq-plot
qqplot_custom <- ggplot(qq_data) +
  stat_qq(aes(sample=Residuals), color = '#377185', alpha = .7) +
  stat_qq_line(aes(sample=Residuals), color = 'black') +
  theme_classic() +
  labs(y="Sample Quantiles",x="Theoretical Quantiles",
       title = 'qq-plot: 
GAM of Dino derivative ~ Stratification Index')

qqplot_custom

# Save the plot
# ggsave('Plots/GAMs/Drivers/qqplot_custom_GAM_Stratification_Index.tiff', dpi = 300, height = 164, width = 164,
#                units = 'mm', compression = 'lzw')

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data)+
  geom_point(aes(x=Fitted, y=Residuals), 
             alpha = .7, color = '#377185') +
  theme_classic() +
  labs(y="Residuals",x="Fitted Values",
       title = 'Residuals vs Fitted: 
GAM of Dinophysis derivative ~ Stratification Index')

RvFplot_custom

# Save the plot
# ggsave('Plots/GAMs/Drivers/RvFplot_custom_GAM_Stratification_Index.tiff', dpi = 300, height = 120, width = 164,
#                units = 'mm', compression = 'lzw')

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data, aes(x = Residuals)) +
  geom_histogram(binwidth = .005, fill = '#377185') +
  theme_classic() +
  labs(x='Residuals', y = 'Count',
       title = 'Histogram of residuals: 
GAM of Dinophysis derivative ~ Stratification_Index')

HistRes_custom

# Save the plot
# ggsave('Plots/GAMs/Drivers/HistRes_custom_GAM_Stratification_Index.tiff', dpi = 300, height = 164, width = 164,
#                units = 'mm', compression = 'lzw')

# Plot the model against the data
ggplot(data = gam_Dino_newdata,
       aes(x = Stratification_Index, y = fit_resp)) +
  geom_line(color = '#377185', linewidth = .5) +
  geom_ribbon(aes(x = Stratification_Index,
                  ymin = right_lwr, 
                  ymax = right_upr), fill = '#377185', color = '#377185', alpha = .2) +
  geom_point(data = Table_GAM, aes(x = Stratification_Index, y = .derivative), size = .8, 
             alpha = .5, color = '#377185') +
  labs(x = 'Surface Solar Radiation', 
       y = 'Dinophysis accumulation rate (d-1)',
       title = 'GAM fit: Dinophysis derivative ~ Stratification_Index') +
  theme_classic()

# Save this plot
# ggsave('Plots/GAMs/Drivers/GAM_fit_Stratification_Index.tiff', height = 135, width = 164,
#        units = 'mm', compression = 'lzw')

## NO3 + NO2 ####

# Now we're gonna check the variables of our second random forest (the one with
# only 4 sites and 7 years)

# We'll start with the most important variables (as given by the VIP analysis) :
#NO3 + NO2

# Restrict our dataset
Table_GAM <- Table_data_RF_multiyear %>%
  filter(Code_point_Libelle %in% c('Antifer ponton pétrolier', 'Cabourg',
                                   'Men er Roue', 'Ouest Loscolo')) %>%
  filter(Year >= 2016 & Year <= 2022)

# New color palette
pheno_palette4 <- c('red3', 'orangered', '#2156A1', '#5995E3')

# Let's plot the data
ggplot(data = Table_GAM,
       aes(x = NO3.NO2, y = .derivative, color = Code_point_Libelle)) +
  # Points
  geom_point(size = 2, alpha = .7) +
  # Add 2 "ghost points" to force the minimum y-axis range to -0.1;0.1
  geom_point(aes(x = 0, y = -.1), color = 'transparent', fill = 'transparent') +
  geom_point(aes(x = 0, y = .1), color = 'transparent', fill = 'transparent') +
  # Color scale (sampling sites)
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  # New color scale for NO3.NO2erature rug plot
  new_scale_color() +
  # geom_rug(data = Table_data_RF_multiyear, aes(x = NO3.NO2, y = NULL,
  #                                              color = NO3.NO2)) +
  scale_color_distiller(palette = 'RdBu', direction = -1, guide = 'none') +
  # facets
  facet_wrap(facets = 'Code_point_Libelle', scales = 'free_y') +
  # Labels
  labs(x = 'N03 + NO2 (µmol.L-1)', y = 'GAM derivative',) +
  # Theme
  theme_classic()

# Once again, close to no signal.
# Let's check with a GAM

# We'll do this with the 10 sites we used for the random forest (Teychan and 
# Auger excluded)

Table_GAM <- Table_GAM %>%
  filter(is.na(NO3.NO2) == FALSE)

gam_Dino <- gam(data = Table_GAM, 
                # Only a spline for the NO3.NO2erature
                # 'k = -1' allows the model to fix the 'best' number of basic
                # functions (= knots)
                formula = .derivative~s(NO3.NO2, k = -1),
                # Introducing the weights
                # weights = unif_weight,
                # Using a Poisson distribution for count data
                family = gaussian(),
                # Restricted maximum likelihood estimation (recommended method)
                method = 'REML')

summary(gam_Dino)
gam.check(gam_Dino)

# Create a new 'empty' dataset for storing model prediction (fit)
gam_Dino_newdata <- expand_grid(NO3.NO2=seq(min(Table_GAM$NO3.NO2),
                                                         max(Table_GAM$NO3.NO2), length.out = 300))

## Get the inverse link function of the model
# With this function, we can transform the prediction we make on the link
#scale to the response scale
ilink <- gam_Dino$family$linkinv

## Predict : add fit and se.fit on the **link** scale
gam_Dino_newdata <- bind_cols(gam_Dino_newdata, 
                              setNames(as_tibble(
                                predict(gam_Dino, gam_Dino_newdata, 
                                        se.fit = TRUE, type = 'link')[1:2]),
                                c('fit_link','se_link')))

## Create the 95% confidence interval (2*standard error fit) and backtransform 
# to response variable using the inverse link function
gam_Dino_newdata <- mutate(gam_Dino_newdata,
                           fit_resp  = ilink(fit_link),
                           right_upr = ilink(fit_link + (2 * se_link)),
                           right_lwr = ilink(fit_link - (2 * se_link)))
# Check the confidence interval. It should not extend below 0 (negative counts
# are impossible)
min(gam_Dino_newdata$right_lwr) # Nice :)
min(gam_Dino_newdata$fit_resp)
max(gam_Dino_newdata$right_upr) # Nice too
max(gam_Dino_newdata$fit_resp)
# The maximum value seems really small

### Checking the model
ModelOutputs<-data.frame(Fitted=fitted(gam_Dino),
                         Residuals=resid(gam_Dino))

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- gam_Dino$model
# then we add the values of fitted and residuals 
# (but are they in the same order as the model? -> need to check that)
qq_data <- bind_cols(qq_data, ModelOutputs)

# Plot : verify that true data matches (more or less) model fit
ggplot(qq_data)+
  geom_point(aes(x = NO3.NO2, y = .derivative), color = 'red') +
  geom_point(aes(x = NO3.NO2, y = Fitted), color = 'blue') +
  theme_classic() +
  labs(y = "Dinophysis derivative", x = "NO3.NO2")

# The data is almost symmetrical around 0, the fit is therefore near a perfect 0

# Custom qq-plot
qqplot_custom <- ggplot(qq_data) +
  stat_qq(aes(sample=Residuals), color = '#FC4D6B', alpha = .7) +
  stat_qq_line(aes(sample=Residuals), color = 'black') +
  theme_classic() +
  labs(y="Sample Quantiles",x="Theoretical Quantiles",
       title = 'qq-plot: 
GAM of Dino derivative ~ NO3+NO2')

qqplot_custom

# Save the plot
# ggsave('Plots/GAMs/Drivers/qqplot_custom_GAM_NO3.NO2.tiff', dpi = 300, height = 164, width = 164,
#                units = 'mm', compression = 'lzw')

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data)+
  geom_point(aes(x=Fitted, y=Residuals), 
             alpha = .7, color = '#FC4D6B') +
  theme_classic() +
  labs(y="Residuals",x="Fitted Values",
       title = 'Residuals vs Fitted: 
GAM of Dinophysis derivative ~ NO3+NO2')

RvFplot_custom

# Save the plot
# ggsave('Plots/GAMs/Drivers/RvFplot_custom_GAM_NO3.NO2.tiff', dpi = 300, height = 120, width = 164,
#                units = 'mm', compression = 'lzw')

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data, aes(x = Residuals)) +
  geom_histogram(binwidth = .005, fill = '#FC4D6B') +
  theme_classic() +
  labs(x='Residuals', y = 'Count',
       title = 'Histogram of residuals: 
GAM of Dinophysis derivative ~ NO3.NO2')

HistRes_custom

# Save the plot
# ggsave('Plots/GAMs/Drivers/HistRes_custom_GAM_NO3.NO2.tiff', dpi = 300, height = 164, width = 164,
#                units = 'mm', compression = 'lzw')

# Plot the model against the data
ggplot(data = gam_Dino_newdata,
       aes(x = NO3.NO2, y = fit_resp)) +
  geom_line(color = '#FC4D6B', linewidth = .5) +
  geom_ribbon(aes(x = NO3.NO2,
                  ymin = right_lwr, 
                  ymax = right_upr), fill = '#FC4D6B', color = '#FC4D6B', alpha = .2) +
  geom_point(data = Table_GAM, aes(x = NO3.NO2, y = .derivative), size = .8,
             alpha = .5, color = '#FC4D6B') +
  labs(x = 'NO3 + NO2', 
       y = 'Dinophysis accumulation rate (d-1)',
       title = 'GAM fit: Dinophysis derivative ~ NO3+NO2') +
  theme_classic()

# Save this plot
# ggsave('Plots/GAMs/Drivers/GAM_fit_NO3.NO2.tiff', height = 135, width = 164,
#        units = 'mm', compression = 'lzw')
