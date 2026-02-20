#### Script: Visualisations of multiyear data on Dinophysis phenology ###
# Author: V. POCHIC
# Last modif: 2026/02/18

## General information ####

## Description: 
# This script produces visualisations of multiyear data/trends in Dinophysis
# phenology based on REPHY data. It also contains an analysis of trends in
# Dinophysis maxima with Generalised Linear Models (GLMs)

## Files required:
# Data/REPHY_outputs/Season_Dino.csv
# Data/REPHY_outputs/Table_hydro.csv
## --> REPHY data to be plotted
# Data/GAM_outputs/Dinophysis/All_sites/Dino_GAM_response_pred.csv
## --> GAM fit of Dinophysis phenology (function of DOY)
# Data/Models/ERA5/Outputs/...
# Data/Models/GAMAR/Outputs/Stratif_index_GAMAR_12sites.csv
## --> Data from oceanographic/atmospheric physical models

## Outputs:
## Files:
# Data/GLM_outputs/Table_S1_summary_glm_trends_beta.csv
## --> Summary of the GLM for multiannual trends of timing of Dinophysis maxima

## Figures:
# Figure 4 - Timing of Dinophysis maxima over the study period (+ GLM fit)
# Figure 5 - Timing of Dinophysis maxima + heatmap of environmental variables
# Figure S6 - Timing of Dinophysis maxima over latitude
# Figure S7 - Timing of Dinophysis maxima + heatmap of ssr

#### Required packages ####

library(tidyverse)
library(tidymodels)
library(ggplot2)
library(mgcv)
library(gratia)
library(gcplyr)
library(ggnewscale)
library(cmocean)
library(RColorBrewer)
library(ggpubr)

####-----------------------------------------------------------------------####
#### Evolution of the timing of maximum count over the study period ####

# One question is: did the phenology of Dinophysis evolve over the study
# period (2007-2022)?
# The GAM we designed is not adequate to answer this question. A **very basic** 
# way to investigate it is to look at the date of the maximum for each site, for
# each year, and check if there is a trend.

#### Import data ####

# Import seasonalised Dinophysis data
Season_Dino <- read.csv2('Data/REPHY_outputs/Season_Dino.csv', 
                         fileEncoding = "ISO-8859-1")

# Import fit of Dinophysis GAM
Dino_response_pred <- read.csv2(
  'Data/GAM_outputs/Dinophysis/All_sites/Dino_GAM_response_pred_multiyear.csv',
  header = TRUE, fileEncoding = 'ISO-8859-1')

Season_Dino_nozeros <- Season_Dino %>%
  filter(true_count != 0)

#### Maxima ####

# Ok so now we look at maxima of Dinophysis counts for each site.

## However there's a catch: the GAM showed us that many sites have 2 peaks in 
# their phenology (all except Seine Bay sites), let's split the year in 2 to get 
# both peaks. Our splitting point will be the minimum of the GAM between the 2 
# peaks.

# We're going to identify the turning point in the GAM: the Day for which 
# fit(Day)<fit(Day-1) AND fit(Day)<fit(Day+1)

# Because the phenology of the GAM should be the same for all years, we should
# get the same date (at each site) whatever the year we consider. But we're 
# going to check that.

GAM_pit <- Dino_response_pred %>%
  # Group by site and Year
  group_by(Code_point_Libelle, Year) %>%
  # Filter rows for which the value of fit is inferior to both the previous
  # (lag) and following (lead) value
  filter(lag(fit) > fit & lead(fit) > fit)

# Here we do the opposite: we're looking for the peaks
GAM_peak <- Dino_response_pred %>%
  # Group by site and Year
  group_by(Code_point_Libelle, Year) %>%
  # Filter rows for which the value of fit is inferior to both the previous
  # (lag) and following (lead) value
  filter(lag(fit) < fit & lead(fit) < fit)

# Cleaning up: filter out Seine Bay sites (only 1 peak) and sampling sites where
# there is almost no Dinophysis
GAM_pit_select <- filter(GAM_pit, Code_point_Libelle %in% 
                         c('Men er Roue', 'Ouest Loscolo',
                         'Le Cornard', 'Auger',
                         'Arcachon - Bouée 7', 'Teychan bis',
                         'Bouzigues (a)', 'Parc Leucate 2', 'Sète mer')) %>%
# We see there are 2 pits for many sites: 1 corresponds to the pit between the 
# 2 peaks, and 1 to the "renewal" after the winter. We will select only the 
# former (Day > 70 and < 300). This works for all sites except for Parc Leucate.
# We create a special condition for Parc Leucate (pit at day 55)
  filter(ifelse(Code_point_Libelle != 'Parc Leucate 2', Day > 70 & Day < 300,
                Day == 55)) %>%
  # Reorder the site as a factor
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Bouzigues (a)', 'Parc Leucate 2',
                                          'Sète mer'))

# For peaks, filter out sites with almost no Dino
GAM_peak_select <- filter(GAM_peak, Code_point_Libelle %in% 
                           c('Antifer ponton pétrolier', 'Cabourg',
                             'Men er Roue', 'Ouest Loscolo',
                             'Le Cornard', 'Auger',
                             'Arcachon - Bouée 7', 'Teychan bis',
                             'Bouzigues (a)', 'Parc Leucate 2',
                             'Sète mer')) %>%
  # There's a "rogue" maximum in Teychan in winter due to noise : get rid of it
  filter(Day != 359) %>%
  # Reorder the site as a factor
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Bouzigues (a)', 'Parc Leucate 2',
                                          'Sète mer'))

# Next step: split the maxima in 2 periods (before and after the pit)

# Getting the Dinophysis maxima in REPHY data
Maxima_Dino <- Season_Dino_nozeros %>%
  # group
  group_by(Code_point_Libelle, Year) %>%
  # We'll need to create a period variable (1/2) for each site
  # This piece of code is not pretty but it works
  mutate(period = 
           #Men er Roue
           ifelse(Code_point_Libelle == 'Men er Roue' & Day <= 223,
                  1,
                  ifelse(Code_point_Libelle == 'Men er Roue',
                         2,
           # Ouest Loscolo
           ifelse(Code_point_Libelle == 'Ouest Loscolo' & Day <= 232,
                  1,
                  ifelse(Code_point_Libelle == 'Ouest Loscolo',
                         2,
           # Le Cornard
           ifelse(Code_point_Libelle == 'Le Cornard' & Day <= 210,
                  1,
                  ifelse(Code_point_Libelle == 'Le Cornard',
                         2,
           # Auger
           ifelse(Code_point_Libelle == 'Auger' & Day <= 209,
                  1,
                  ifelse(Code_point_Libelle == 'Auger',
                         2,
           # Arcachon Bouée 7
           ifelse(Code_point_Libelle == 'Arcachon - Bouée 7' & Day <= 195,
                  1,
                  ifelse(Code_point_Libelle == 'Arcachon - Bouée 7',
                         2,
           # Teychan bis
           ifelse(Code_point_Libelle == 'Teychan bis' & Day <= 201,
                  1,
                  ifelse(Code_point_Libelle == 'Teychan bis',
                         2,
           # Parc Leucate 2
           ifelse(Code_point_Libelle == 'Parc Leucate 2' & Day <= 234 
                  # Need to add this condition because of winter peak continuing
                  # into the next year
                  & Day >= 55,
                  1,
                  ifelse(Code_point_Libelle == 'Parc Leucate 2',
                         2,
          # Bouzigues (a)
           ifelse(Code_point_Libelle == 'Bouzigues (a)' & Day <= 232,
                  1,
                  ifelse(Code_point_Libelle == 'Bouzigues (a)',
                         2,
          # Sète mer
           ifelse(Code_point_Libelle == 'Sète mer' & Day <= 209,
                  1,
                  ifelse(Code_point_Libelle == 'Sète mer',
                         2,
          # Diana centre
           ifelse(Code_point_Libelle == 'Diana centre' & Day <= 164,
                  1,
                  ifelse(Code_point_Libelle == 'Diana centre',
                         2,
           # All other sites
           1)))))))))))))))))))))

# Let's compute the number of observations of Dinophysis by period
Maxima_Dino_n <- Maxima_Dino %>%
  # groups
  group_by(Code_point_Libelle, Year, period) %>%
  # summarise to get n
  summarise(n_max = n(), max_Dino = max(true_count), .groups = 'keep') %>%
  # Create a "shape" variable for plotting
  mutate(shape = ifelse(period == 1 & n_max > 1 & max_Dino > 1, 1, 
                        ifelse(period == 2 & n_max > 1 & max_Dino > 1, 2,
                               ifelse(period == 1 & (n_max == 1 || max_Dino == 1), 3,
                                      4))))

# Now extract maxima. Careful, this step has to be done AFTER the previous one.
Maxima_Dino <- Maxima_Dino %>%
  # Add period as a grouping variable
  group_by(Code_point_Libelle, Year, period) %>%
  # extract rows with yearly maxima of 'true_count' (slice() retains groups)
  slice(which.max(true_count)) %>%
# Keep only Dinophysis cell densities >= 100 cells/L (more than 1 cell counted)
filter(true_count >= 1)

# Join both datasets to get n() info on maxima
Maxima_Dino <- left_join(Maxima_Dino, Maxima_Dino_n) %>%
  # Code_point_Libelle as factor
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

# Seems perfect.

# GAM fit
GAM_peak_select <- GAM_peak_select %>%
  # group
  group_by(Code_point_Libelle) %>%
  # We'll need to create a period variable (1/2) for each site, based on the
  # pits of the GAM
  # This piece of code is not pretty but it works
  mutate(period = 
           #Men er Roue
           ifelse(Code_point_Libelle == 'Men er Roue' & Day <= 223,
                  1,
                  ifelse(Code_point_Libelle == 'Men er Roue',
                         2,
                         # Ouest Loscolo
                         ifelse(Code_point_Libelle == 'Ouest Loscolo' & Day <= 232,
                                1,
                                ifelse(Code_point_Libelle == 'Ouest Loscolo',
                                       2,
                                       # Le Cornard
                                       ifelse(Code_point_Libelle == 'Le Cornard' & Day <= 210,
                                              1,
                                              ifelse(Code_point_Libelle == 'Le Cornard',
                                                     2,
                                                     # Auger
                                                     ifelse(Code_point_Libelle == 'Auger' & Day <= 209,
                                                            1,
                                                            ifelse(Code_point_Libelle == 'Auger',
                                                                   2,
                                                                   # Arcachon Bouée 7
                                                                   ifelse(Code_point_Libelle == 'Arcachon - Bouée 7' & Day <= 195,
                                                                          1,
                                                                          ifelse(Code_point_Libelle == 'Arcachon - Bouée 7',
                                                                                 2,
                                                                                 # Teychan bis
                                                                                 ifelse(Code_point_Libelle == 'Teychan bis' & Day <= 201,
                                                                                        1,
                                                                                        ifelse(Code_point_Libelle == 'Teychan bis',
                                                                                               2,
                                                                                               # All other sites
                                                                                               1))))))))))))) %>%
  # Add period as a grouping variable
  group_by(Code_point_Libelle, Year, period) %>%
  # Code_point_Libelle as factor
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis'))

# Getting some stats for this second thing
Maxima_Dino_stats <- Maxima_Dino %>%
  group_by(Code_point_Libelle, period) %>%
  summarise(median_daymax = median(Day), mean_daymax = mean(Day), 
            stdev_daymax = sd(Day), min_daymax = min(Day), max_daymax = max(Day),
            median.lat = median(Latitude),
            .groups = 'keep')

# Let's plot it

# Beautiful color palette
pheno_palette16 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646',
                     '#642C3A', '#DEB1CC', '#FC4D6B', '#791D40')


# plot
ggplot(Maxima_Dino, 
       aes(x = Year, y = Day, color = Code_point_Libelle, 
           shape = as.factor(period))) +
  # then the points
  geom_point(size = 4, alpha = .8) +
  # new color scale for the points
  scale_color_discrete(type = pheno_palette16, guide = 'none') +
  # aesthetics
  facet_wrap(facets = 'Code_point_Libelle', nrow = 4) +
  scale_y_continuous(limits = c(1, 365), 
                     breaks = c(1, 100, 200, 300, 365)) +
  scale_shape_discrete(guide = 'none') +
  # scale_shape_manual(values = c(1, 3), guide = 'none') +
  # labels
  labs(y = 'Day of maximum Dinophysis count') +
  theme_classic()

# We have some problems ar Parc Leucate 2: the maxima are cut in 2 groups
# because they occur around the 1st of January.
# To fix this, we will need to reorder the seasonality so that maxima in winter 
# don't get separated artificially. For this we will manipulate the data to
# create a year that's more than 365 days (yeah, a bit weird.)

#### Statistical evaluation of trends: GLM ####

# We want to test statistically if there are linear trends in the timing of 
# Dinophysis maxima

# Very simple: our response variable is day of Dino max (DOY), our predictor is
# the Year

# We'll do GLMs, with the help of the tutorial by Dr. Bede Davies
# https://bedeffinianrowedavies.com/statisticstutorials/poissonglms

# We still have an issue with the winter peak at Parc Leucate 2. We'll 
# circumvent it with a trick: create a 'NewDay' variable that lengthens the Year
# So that very early January dates occur after December dates.
Maxima_Dino_trends <- Maxima_Dino %>%
  # Only interesting sampling sites
  filter(Code_point_Libelle %in% 
    c('Antifer ponton pétrolier', 'Cabourg',
    'Men er Roue', 'Ouest Loscolo',
    'Le Cornard', 'Auger',
    'Arcachon - Bouée 7', 'Teychan bis',
    'Parc Leucate 2', 'Bouzigues (a)', 'Sète mer',
    'Diana centre')) %>%
  # Here is the 'NewDay' trick:
  mutate(NewDay = ifelse(Day <= 35, Day + 365, Day)) %>%
  # create a y-trend variable to use as the predictor
  # that's equal to NewDay in the problematic situation we identified,
  # (Parc Leucate 2) or equal to Day for everything else
  mutate(y_trend = ifelse(Code_point_Libelle == 'Parc Leucate 2' &
                            period == 2, NewDay, Day)) %>%
  # And create a category variable that is the sampling site + the period
  mutate(category = paste(Code_point_Libelle, period, sep = "_")) %>%
  # Relevel as a factor (category = Sampling site + period)
  mutate(category = as_factor(category)) %>%
  mutate(category = fct_relevel(category,
                                        'Antifer ponton pétrolier_1', 'Cabourg_1',
                                        'Men er Roue_1', 'Men er Roue_2',
                                        'Ouest Loscolo_1', 'Ouest Loscolo_2',
                                        'Le Cornard_1', 'Le Cornard_2', 
                                        'Auger_1', 'Auger_2',
                                        'Arcachon - Bouée 7_1', 'Arcachon - Bouée 7_2', 
                                        'Teychan bis_1', 'Teychan bis_2',
                                        'Parc Leucate 2_1', 'Parc Leucate 2_2',
                                        'Bouzigues (a)_1', 'Bouzigues (a)_2',
                                        'Sète mer_1', 'Sète mer_2',
                                        'Diana centre_1', 'Diana centre_2'))

# Let's transform our response variable so it follows a beta distribution
# (bounded between 0 and 1). This allows us to fit a GLM properly
Maxima_Dino_trends <- Maxima_Dino_trends %>%
  # For this we divide y_trend (our slightly modified DOY) by the maximum value
  # of this variable
  mutate(y_trend_beta = y_trend/max(Maxima_Dino_trends$y_trend))

# Now, we'll fit a GLM for every series of maxima (1 for each period in each 
# sampling site)
# With the beta-transformed variable (we'll use mgcv's gam() function)
glm_trends_beta <- gam(y_trend_beta~Year*category,
                       data = Maxima_Dino_trends,
                       # Specifying a beta distribution of the variable
                       family = betar(link="logit"))

summary(glm_trends_beta)

### Diagnostic plots ###

### Checking the model
ModelOutputs<-data.frame(Fitted=fitted(glm_trends_beta),
                         Residuals=resid(glm_trends_beta))

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- glm_trends_beta$model
# then we add the values of fitted and residuals
qq_data <- bind_cols(qq_data, ModelOutputs)

# Reorder data
qq_data_reordered <- qq_data %>%
  mutate(category = fct_relevel(category,
                                          'Antifer ponton pétrolier_1', 'Cabourg_1',
                                          'Men er Roue_1', 'Men er Roue_2',
                                          'Ouest Loscolo_1', 'Ouest Loscolo_2',
                                          'Le Cornard_1', 'Le Cornard_2',
                                          'Auger_1', 'Auger_2',
                                          'Arcachon - Bouée 7_1', 'Arcachon - Bouée 7_2',
                                          'Teychan bis_1', 'Teychan bis_2',
                                          'Parc Leucate 2_1', 'Parc Leucate 2_2',
                                          'Bouzigues (a)_1', 'Bouzigues (a)_2',
                                          'Sète mer_1', 'Sète mer_2', 
                                          'Diana centre_1', 'Diana centre_2')) %>%
  # we ungroup the data frame to produce only 1 qq-plot
  ungroup()

# And (qq-)plot
qqplot_custom <- ggplot(qq_data_reordered) +
  stat_qq(aes(sample=Residuals, color = category), alpha = .7) +
  stat_qq_line(aes(sample=Residuals, color = category)) +
  facet_wrap(facets = c('category')) +
  scale_color_discrete(guide = 'none') +
  theme_classic() +
  labs(y="Sample Quantiles",x="Theoretical Quantiles")

qqplot_custom

# Ok!

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data_reordered)+
  geom_point(aes(x=Fitted,y=Residuals, color = category), 
             alpha = .7) +
  facet_wrap(facets = c('category'), scales = 'free') +
  scale_color_discrete(guide = 'none') +
  theme_classic() +
  labs(y="Residuals",x="Fitted Values")

RvFplot_custom

# Okish

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data_reordered, aes(x = Residuals, 
                                                fill = category))+
  geom_histogram(binwidth = 1)+
  facet_wrap(facets = c('category'), scales = 'free') +
  scale_fill_discrete(guide = 'none') +
  theme_classic() +
  labs(x='Residuals', y = 'Count')

HistRes_custom

# For some sites it's not that great...
# For others it's very good

# Exporting the summary of the GLM as a table
summary_glm_trends_beta <- tidy(glm_trends_beta, parametric = TRUE, smooth = FALSE)

# write.csv2(summary_glm_trends_beta, 'Data/GLM_outputs/Table_S1_summary_glm_trends_beta.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')

# Let's predict with the GLM
NewData_1 <- expand_grid(# Year variable
                         Year=seq(min(Maxima_Dino_trends$Year), 
                                  max(Maxima_Dino_trends$Year)),
                         # We also add the site data
                         category = unique(Maxima_Dino_trends$category))

Pred <- predict(glm_trends_beta, NewData_1, se.fit=TRUE, type="response")

NewData<-NewData_1 %>% 
  mutate(response=Pred$fit,
         se.fit=Pred$se.fit,
         Upr=response+(se.fit*1.96),
         Lwr=response-(se.fit*1.96))

# Join datasets of prediction and true data
Data_trends <- left_join(Maxima_Dino_trends, NewData)

# Color palette for 12 sampling sites
pheno_palette12 <- c('red3', 'orangered', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646',
                     '#642C3A', '#DEB1CC', '#FC4D6B', '#791D40')

# Plot (we multiply the predictions made on the beta-transformed variable
# by a factor so that it fits the true data)
ggplot(Data_trends)+
  # confidence interval of model
  geom_ribbon(aes(x = Year,
                  ymax = Upr*max(Maxima_Dino_trends$y_trend),
                  ymin = Lwr*max(Maxima_Dino_trends$y_trend),
                  fill=Code_point_Libelle, 
                  group = category),
              alpha=0.35) +
  # fit of model
  geom_line(aes(x = Year,
                y = response*max(Maxima_Dino_trends$y_trend),
                color=Code_point_Libelle,
                group = category)) +
  # true data
  geom_point(aes(x = Year, y = y_trend, color = Code_point_Libelle, 
                 shape = as_factor(shape))) +
  # color scales
  scale_color_discrete(palette = pheno_palette12, guide = 'none') +
  scale_fill_discrete(palette = pheno_palette12, guide = 'none') +
  # shape scale
  scale_shape_manual(values = c(16, 17, 1, 2), guide = 'none') +
  # y-axis scale
  scale_y_continuous(limits = c(1, 410), 
                     breaks = c(1, 100, 200, 300, 365)) +
  # Labels
  labs(x = "Year", 
       y = "Day of maximum Dinophysis count") +
  facet_wrap(facets = 'Code_point_Libelle') +
  theme_classic() +
  theme(strip.text = element_text(
    size = 7, color = 'black'))

# Good sh*t
# ggsave('Plots/REPHY/Fig4_Trends_maxima_Dino.tiff', dpi = 300, height = 150, width = 164,
#                units = 'mm', compression = 'lzw')

####-----------------------------------------------------------------------####
#### Environmental data ####

# In this section, we want to plot Dinophysis maxima in relation to several
# environmental variables, to see if there are patterns between those and the
# timing of the maximum.

#### Import data ####

## Here we will call a table of environmental data to plot a heatmap of
# environmental variables underneath our Dinophysis data.

#--- Hydrology data from REPHY dataset ---#

Table_hydro <- read.csv2('Data/REPHY_outputs/Table_hydro.csv', 
                               header = TRUE,
                               fileEncoding = 'ISO-8859-1')

# However, we will not plot the data as raw values. Instead, we want to get a
# median by fortnight, and plot that.
# How do we do that? First: summarise by fortnight.
Table_hydro_fortnight <- Table_hydro %>%
  # Restrict to sampling sites and period of interest
  # (We don't consider Mediterranean sampling sites)
  filter(Code_point_Libelle %in% c('Point 1 Boulogne', 'At so',
                                   'Antifer ponton pétrolier', 'Cabourg',
                                   'les Hébihens', 'Loguivy',
                                   'Men er Roue', 'Ouest Loscolo',
                                   'Le Cornard', 'Auger',
                                   'Arcachon - Bouée 7', 'Teychan bis')) %>%
  filter(Year >= 2007) %>%
  # Select only relevant variables
  select(c('Code_point_Libelle', 'lon', 'lat', 'Year', 'Month', 'Date', 
           'ID.interne.passage', 'CHLOROA', 'SALI', 'TEMP')) %>%
  # Create a DOY, week, and fortnight variables
  mutate(Date = ymd(Date)) %>%
  mutate(Day = yday(Date)) %>%
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2)) %>%
  # Summarise by fortnight! We calculate the median by fortnight for [chl a],
  # salinity and temperature.
  group_by(Code_point_Libelle, lon, lat, Fortnight) %>%
  summarise(CHLOROA.med = median(CHLOROA, na.rm = TRUE), 
            SALI.med = median(SALI, na.rm = TRUE),
            TEMP.med = median(TEMP, na.rm = TRUE), .groups = 'keep')
  

# Now we need to expand these fortnightly medians into a daily table, for 
# plotting.

# We create a complete table of DOY:
Table_daily <- expand_grid(Day=seq(1, 365),
     # We also add the sampling site
     Code_point_Libelle = unique(Table_hydro_fortnight$Code_point_Libelle)) %>%
  # And now, we compute the fortnight
  mutate(Fortnight = ceiling(Day/14))

# Aaaaaand, we left_join both tables
Table_hydro_daily <- left_join(Table_daily, Table_hydro_fortnight) %>%
  # The following lines put the table in a more human-readable format
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis')) %>%
  arrange(Code_point_Libelle) %>%
  # And we get rid of fortnight 27 because we don't have enough data
  filter(Fortnight < 27)

# Great!

#--- Stratification data from GAMAR model ---#

## Now, we will import the stratification data from the GAMAR model
Table_stratif <- read.csv2('Data/Models/GAMAR/Outputs/Stratif_index_GAMAR_12sites.csv', 
          header = TRUE, fileEncoding = 'ISO-8859-1')

# Same sh*t as before: need to summarise this into fortnightly medians. Proceed
# exactly as before
Table_stratif_fortnight <- Table_stratif %>%
  # Select only relevant variables
  select(c('Code_point_Libelle', 'Target_Date', 'Stratification_Index')) %>%
  # Create a DOY, week, and fortnight variables
  mutate(Date = ymd_hms(Target_Date)) %>%
  mutate(Day = yday(Date)) %>%
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2)) %>%
  # Summarise by fortnight! We calculate the median by fortnight for the
  # Stratification Index
  group_by(Code_point_Libelle, Fortnight) %>%
  summarise(SI.med = median(Stratification_Index, na.rm = TRUE),
            .groups = 'keep')

# Left_join with the daily table previously defined
Table_stratif_daily <- left_join(Table_daily, Table_stratif_fortnight) %>%
  # The following lines put the table in a more human-readable format
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis')) %>%
  arrange(Code_point_Libelle) %>%
  # And we get rid of fortnight 27 because we don't have enough data
  filter(Fortnight < 27)

# Excellent

#--- Surface solar radiation (ssr) data from ERA5 model ---#

## This one is already summarised by fortnight, which is great, but we need to
# take the median over several years
Table_era5_ssr_fortnight <- read.csv2(
  'Data/Models/ERA5/Outputs/era5_dataset_ssr_fortnight_median.csv', 
  header = TRUE, fileEncoding = 'ISO-8859-1') %>%
  # Take the median over several years
  group_by(Fortnight, Code_point_Libelle) %>%
  summarise(ssr.med = median(ssr, na.rm = TRUE), 
            .groups = 'keep')

Table_ssr_daily <- left_join(Table_daily, Table_era5_ssr_fortnight,
                                 by = c('Fortnight', 'Code_point_Libelle')) %>%
  # The following lines put the table in a more human-readable format
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis')) %>%
  arrange(Code_point_Libelle) %>%
  # And we get rid of fortnight 27 because we don't have enough data
  filter(Fortnight < 27)

# Perfect!!!

#### Plotting ####

#--- Selection of sampling sites ---#

# We'll create a sub-selection of Channel/Atlantic sites, without Teychan bis
# for which stratification data from the GAMAR model are too far from ground
# truth values (because it is located in a semi-enclosed shallow bay)
Table_hydro_daily_select <- filter(Table_hydro_daily,
                                   Code_point_Libelle %in%
                                     c('Point 1 Boulogne', 'At so',
                                       'les Hébihens', 'Loguivy',
                                       'Antifer ponton pétrolier', 'Cabourg',
                                       'Men er Roue', 'Ouest Loscolo',
                                       'Le Cornard', 'Auger',
                                       'Arcachon - Bouée 7')) %>%
  # Recode the Code_point_Libelle
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7'))

# Stratification
Table_stratif_daily_select <- filter(Table_stratif_daily,
                                 Code_point_Libelle %in%
                                   c('Point 1 Boulogne', 'At so',
                                     'les Hébihens', 'Loguivy',
                                     'Antifer ponton pétrolier', 'Cabourg',
                                     'Men er Roue', 'Ouest Loscolo',
                                     'Le Cornard', 'Auger',
                                     'Arcachon - Bouée 7')) %>%
  # Recode the Code_point_Libelle
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7'))

# SSR
Table_ssr_daily_select <- filter(Table_ssr_daily,
                             Code_point_Libelle %in%
                               c('Point 1 Boulogne', 'At so',
                                 'les Hébihens', 'Loguivy',
                                 'Antifer ponton pétrolier', 'Cabourg',
                                 'Men er Roue', 'Ouest Loscolo',
                                 'Le Cornard', 'Auger',
                                 'Arcachon - Bouée 7')) %>%
  # Recode the Code_point_Libelle
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7'))

# Dinophysis maxima
Maxima_Dino_select <- filter(Maxima_Dino,
                         Code_point_Libelle %in%
                           c('Point 1 Boulogne', 'At so',
                             'les Hébihens', 'Loguivy',
                             'Antifer ponton pétrolier', 'Cabourg',
                             'Men er Roue', 'Ouest Loscolo',
                             'Le Cornard', 'Auger',
                             'Arcachon - Bouée 7')) %>%
  # Recode the Code_point_Libelle
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7'))

# GAM peak select table (without 1 site: Teychan)
GAM_peak_plot_select <- select(GAM_peak_select, c(Code_point_Libelle, Day, 
                                                  period, Year)) %>%
  group_by(Code_point_Libelle, period) %>%
  # this summarise function effectively retains only one row, as the Day of
  # maximum (Day) is the same for each year - this is a feature of the GAM
  summarise(Day = mean(Day), .groups = 'keep') %>%
  filter(Code_point_Libelle %in% 
                             c('Antifer ponton pétrolier', 'Cabourg',
                               'Men er Roue', 'Ouest Loscolo',
                               'Le Cornard', 'Auger',
                               'Arcachon - Bouée 7')) %>%
  # shape variable to plot differently depending on the period
  mutate(shape = ifelse(period == 1, 5, 6))

# Very good

#--- Plots ---#




# New color scales
# For all points
pheno_palette11 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003', '#F7B41D')

# For GAM peaks
pheno_palette7 <- c('red3', 'orangered', '#2156A1', '#5995E3', 
                    '#1F3700', '#649003', '#F7B41D')

## Plotting the heatmap, maxima and GAM peaks
# 5 versions of the plot : temperature/chl a/salinity/stratification/ssr heatmaps
# names for the variable that contains the plot: plot_temp, plot_chloroa, 
# plot_sali, plot_stratif, plot_ssr

# !!! Don't forget to change the name !!!

plot_chloroa <- ggplot() +
  geom_tile(data = 
            Table_hydro_daily_select,
            # Table_stratif_daily_select,
            # Table_ssr_daily_select,
            aes(x = Day, y = Code_point_Libelle, 
                fill = 
                  # ssr.med
                  # TEMP.med
                  CHLOROA.med
                  # SALI.med
                  # SI.med
            ),
            alpha = .8) +
  # color palette for fill
  # Stratification version
  # scale_fill_cmocean(name = 'deep', direction = 1) +
  # Temperature version
  # scale_fill_distiller(palette = 'RdBu', direction = -1) +
  # Chl a version
  scale_fill_cmocean(name = 'algae') +
  # Salinity version
  # scale_fill_cmocean(name = 'haline') +
  # ssr version
  # scale_fill_cmocean(name = 'solar') +
  
  # Add labels here so the name of the 'fill' legend is correct
  labs(x = 'Day of the year', y = NULL, fill = 
         # 'Median SST (°C)',
         # 'Mean Surface Solar Radiation (J/m2)',
         'Median [chl a] (mg/m3)',
         # 'Median Salinity (PSU)',
         # 'Median Stratification Index (-)',
       title = NULL) +
  
  # New fill scale for points
  new_scale_fill() +
  
  # All maxima as smaller translucid points
  geom_point(data = Maxima_Dino_select, aes(x = Day, y = Code_point_Libelle, 
                                            color = Code_point_Libelle,
                                        shape = as_factor(shape)), 
             size = 2.5, alpha = .8) +
  # Color scales (Maxima)
  scale_color_discrete(type = pheno_palette11, guide = 'none') +
  # New color scale
  new_scale_fill() +
  
  # points for GAM peak
  geom_point(data = GAM_peak_plot_select, aes(y = Code_point_Libelle, x = Day, 
                 fill = Code_point_Libelle, shape = as_factor(shape)), 
             size = 4, color = 'grey10',stroke = .5) +
  # Color scales (GAM peak)
  scale_fill_discrete(type = pheno_palette7, guide = 'none') +
  
  # Scale for shape (all points) +
  scale_shape_manual(values = c(16, 17, 1, 2, 21, 24), guide = 'none') +
  
  
  # Axis stuff
  # reverse y axis to get sites in the desired order
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(1,365), breaks = c(1, 100, 200, 300, 365)) +
  
  # Theme
  theme(plot.title = element_text(size = 11), 
        # Axis
        axis.title.x = element_text(size=8), 
        axis.title.y =element_text(size=8),
        axis.text = element_text(size=6, color = 'black'),
        # Alternative text axis for paper figure
        axis.text.x = element_text(size=6, color = 'black'),
        axis.text.y = element_blank(),
        # Axis lines
        axis.line.x = element_line(linewidth = .2, color = 'black'),
        axis.line.y = element_line(linewidth = .2, color = 'black'),
        # Legend
        legend.background = element_rect(linewidth = .25, color = 'grey10'),
        legend.title = element_text(size = 6, color = 'grey5'),
        legend.frame = element_rect(linewidth = .25, color = 'grey10'),
        legend.ticks = element_line(linewidth = .15, color = 'grey25'),
        legend.text = element_text(size = 6, color = 'grey5'),
        legend.position = 'bottom',
        # Panel
        panel.grid = element_blank(),
        panel.background = element_blank(),
        # Facet labels
        strip.background = element_rect(fill = 'transparent',
                                        linewidth = 1,
                                        color = 'grey10'),
        strip.text = element_text(color = 'grey5', size = 7.5))

# Call plot (optional)
# plot_temp

# Great! Let's save that (for plot_ssr only)
# ggsave('Plots/REPHY/FigS7_Dino_gam_max_ssr_RF_plot.tiff', dpi = 600, height = 150, width = 100,
#        units = 'mm', compression = 'lzw')

#### Arranging 4 heatmap plots on the same figure ####

# For paper figure: plot 4 heatmaps side by side
plot_4heatmaps <- ggarrange(plot_temp, plot_chloroa, plot_sali, plot_stratif,
                            nrow = 2, ncol = 2,
                            align = 'hv', legend = 'bottom')

plot_4heatmaps

# Saving
# ggsave('Plots/REPHY/Fig5_Dino_4heatmaps_plot.tiff', dpi = 600, height = 160,
#        width = 134, units = 'mm', compression = 'lzw')

#### Plotting maxima by latitude ####

# Last but not least, plot the timing of the maximum against latitude of the
# sampling sites.

# Here we take only Atlantic sites but we can plot Teychan bis
# So color palette for 8 sampling sites for GAM peaks
pheno_palette8 <- c('red3', 'orangered', '#2156A1', '#5995E3', 
                    '#1F3700', '#649003', '#F7B41D', '#FBB646')

# We need to get latitude into our GAM_peak table
# GAM peak table (with Teychan)
GAM_peak_plot <- select(GAM_peak_select, c(Code_point_Libelle, Day, 
                                                  period, Year)) %>%
  group_by(Code_point_Libelle, period) %>%
  # this summarise function effectively retains only one row, as the Day of
  # maximum (Day) is the same for each year - this is a feature of the GAM
  summarise(Day = mean(Day), .groups = 'keep') %>%
  filter(Code_point_Libelle %in% 
           c('Antifer ponton pétrolier', 'Cabourg',
             'Men er Roue', 'Ouest Loscolo',
             'Le Cornard', 'Auger',
             'Arcachon - Bouée 7', 'Teychan bis')) %>%
  # shape variable to plot differently depending on the period
  mutate(shape = ifelse(period == 1, 5, 6))

# Then join with Maxima_Dino table to get latitude
Maxima_Dino_lat <- Maxima_Dino %>%
  ungroup() %>%
  select(c('Code_point_Libelle', 'period', 'Latitude')) %>%
  # summarise (all rows are equal)
  group_by(Code_point_Libelle, period) %>%
  summarise(Latitude = mean(Latitude), .groups = 'keep')
  

GAM_peak_plot <- left_join(GAM_peak_plot, Maxima_Dino_lat, 
                           by = c('Code_point_Libelle', 'period'))

# It's not the prettiest code ever but it works!

#--- Plot ---#
ggplot(data  = subset(Maxima_Dino, 
                      Code_point_Libelle %in% c('Point 1 Boulogne','At so',
                                                'Antifer ponton pétrolier', 'Cabourg',
                                                'les Hébihens', 'Loguivy',
                                                'Men er Roue', 'Ouest Loscolo',
                                                'Le Cornard', 'Auger',
                                                'Arcachon - Bouée 7', 'Teychan bis')), 
       aes(x = Day, y = Latitude, color = Code_point_Libelle, 
           shape = as.factor(shape))) +
  # First the points for maxima of each year
  geom_point(size = 2.5, alpha = .8) +
  
  # Then the bigger points for the GAM
  geom_point(data = GAM_peak_plot, aes(y = Latitude, x = Day, 
                 fill = Code_point_Libelle, shape = as_factor(shape)), 
             size = 4, color = 'grey10', stroke = .5) +
  
  # x scale
  scale_x_continuous(limits = c(1, 365), 
                     breaks = c(1, 100, 200, 300, 365)) +
  
  # Color scale for all
  scale_color_discrete(type = pheno_palette16, guide = 'none') +
  # fill scale for all
  scale_fill_discrete(type = pheno_palette8, guide = 'none') +
  # Scale for shape
  scale_shape_manual(values = c(16, 17, 1, 2, 21, 24), guide = 'none') +
  # labels
  labs(y = 'Latitude (degrees)',
       x = 'Day of maximum Dinophysis count (DOY)') +
  theme_classic()

# Save that good stuff
# ggsave('Plots/REPHY/FigS6_Dino_maxima_by_latitude.tiff', dpi = 300, height = 150, width = 100,
#                units = 'mm', compression = 'lzw')

####-------------------------- End of script ------------------------------####