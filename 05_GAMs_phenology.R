#### Script: Generalised Additive Models (GAMs) of Dinophysis, Mesodinium ###
### and alloxanthin based on REPHY data ###
# Author: V. POCHIC
# Last modif: 2026/02/20

## General information ####

## Description: 
# This script computes Generalised additive models from plankton count data
# (Dinophysis, Mesodinium) and pigment concentration (alloxanthin) from the
# REPHY dataset. It also computes the derivative of the GAM - a variation rate
# of Dinophysis as a function of the day of the year (DOY).

## Files required:
# Data/REPHY_outputs/Season_Dino.csv
# Data/REPHY_outputs/Season_Meso.csv
# Data/REPHY_outputs/Season_Allo.csv
## REPHY data to serve as input for the GAMs

## Outputs:
## Files:
# Data/GAM_outputs/Dinophysis/All_sites/Dino_GAM_response_pred.csv
# Data/GAM_outputs/Dinophysis/4_sites/Dino_GAM_4sites_response_pred.csv
# Data/GAM_outputs/Mesodinium/Meso_GAM_response_pred.csv
# Data/GAM_outputs/Alloxanthin/Allo_GAM_response_pred.csv
## --> All GAM fits for all sites Dino GAM + trophic chain GAMs on less sampling 
## --> sites
# Data/GAM_outputs/Dinophysis/All_sites/Dino_GAM_response_pred_multiyear.csv
# Data/GAM_outputs/Dinophysis/4_sites/Dino_GAM_response_pred_multiyear_4sites.csv
## --> These 2 files are similar to the first 2 but are not grouped by DOY, so 
## --> it keeps the variability due to the Year as random effect.
# Data/GAM_outputs/Dinophysis/All_sites/Dino_GAM_derivative.csv
# Data/GAM_outputs/Dinophysis/4_sites/Dino_GAM_derivative_4sites.csv
## --> Derivative values (variation rate) from the Dino GAMs.

## Figures:
# Figure 3 - Phenology GAM of Dinophysis
# Figure S3 - Dinophysis GAM fit on different years at sampling site Ouest Loscolo
# Figures S13-S18 - diagnostic plots for GAMs
# Figure 6 - Phenology GAMs of alloxanthin/Mesodinium/Dinophysis at sampling
# sites Cabourg and Ouest Loscolo
# Figure S11 - Phenology GAMs of alloxanthin/Mesodinium/Dinophysis at sampling
# sites Antifer and Men er Roue

#### Required packages ####

library(tidyverse)
library(ggplot2)
library(mgcv)
library(gratia)
library(gcplyr)
library(ggpubr)

####-----------------------------------------------------------------------####
#### Import data ####

# All seasonalised data:
Season_Dino <- read.csv2('Data/REPHY_outputs/Season_Dino.csv', header = TRUE,
                         fileEncoding = 'ISO-8859-1')
Season_Meso <- read.csv2('Data/REPHY_outputs/Season_Meso.csv', header = TRUE,
                         fileEncoding = 'ISO-8859-1')
Season_Allo <- read.csv2('Data/REPHY_outputs/Season_Allo.csv', header = TRUE,
                         fileEncoding = 'ISO-8859-1')

# Filter Season_Dino to compute a smaller dataset to follow the other GAMs of 
# the trophic chain (only 4 sampling sites, time period 2016-2022):
Season_Dino_4sites <- filter(Season_Dino,
                             Code_point_Libelle %in% 
                               c('Antifer ponton pétrolier', 'Cabourg',
                                 'Men er Roue', 'Ouest Loscolo')
                             &
                               Year >= 2016)



####-----------------------------------------------------------------------####
#### On running GAM models with this script ####

# Because some steps involve heavy calculations, I recommend running the models
# 1 by 1, saving the results and clearing the workspace between each model.

####-----------------------------------------------------------------------####
#### Dinophysis GAM - all sites ####

#--- Model and basic model checks ---#

# Formulate the GAM

# We will use the Year as a random effect and for this we need it as a factor
Season_Dino_factor <- Season_Dino %>%
  mutate(Year = as_factor(Year)) %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle))


gam_Dino <- gam(data = Season_Dino_factor, 
                # Only a spline for the day of the year
                # The use of a cyclic basis spline helps to make ends meet at the
                # first and last days of the year
                # 'k = -1' allows the model to fix the 'best' number of basic
                # functions (= knots)
                formula = true_count~s(Day, bs = 'cc', k = -1,
                                       # separate each site
                                       by = Code_point_Libelle) 
                # We add the Year as a random effect. This will help to assess and
                # smooth the effects of interannual variability in the phenology
                + s(Year, bs = 're', k = -1),
                # Using a Poisson distribution for count data
                family = poisson(),
                # Restricted maximum likelihood estimation (recommended method)
                method = 'REML')

summary(gam_Dino)

# Create a new 'empty' dataset for storing model prediction (fit)
gam_Dino_newdata <- expand_grid(Day=seq(1, 365),
                                # We add a Year vector as it has become a factor
                                # of the model
                                Year=seq(min(Season_Dino$Year), 
                                        max(Season_Dino$Year)),
                                # We also add the site data
                                Code_point_Libelle = unique(Season_Dino$Code_point_Libelle))

## Get the inverse link function of the model
# With this function, we can transform the prediction we make on the link
#scale to the response scale
ilink <- gam_Dino$family$linkinv

## Predict : add fit and se.fit on the **link** scale
gam_Dino_newdata <- bind_cols(gam_Dino_newdata, 
                              setNames(as_tibble(
                                predict(gam_Dino, gam_Dino_newdata, 
                                        se.fit = TRUE, type = 'link',
                                        re.form = ~ 1|Year)[1:2]),
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
# But we can see that the confidence interval is ridiculously small around the
# min and max of the model fit --> we'll fix that later down the script

# Quick plot

# Get sampling sites in the correct order
gam_Dino_newdata <- gam_Dino_newdata %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

Season_Dino_factor <- Season_Dino_factor %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

# Plot:
ggplot(gam_Dino_newdata, aes(x = Day, y = fit_resp, color = Year))+
  geom_line(linewidth = 1) +
  geom_point(data = Season_Dino, aes(x = Day, y = true_count, color = Year), 
             shape = 21, alpha = .3) +
  facet_wrap(facets = 'Code_point_Libelle', scales = 'free_y') +
  theme_classic()

# Ok nice. So we see that 1) each sampling site has a unique model fit, that
# matches more or less with the data, which is great. 2) the 'Year' random 
# effect is uniform among all sites, that is normal given how we formulated the
# model.

#--- Checking the model with diagnostic plots ---#
ModelOutputs<-data.frame(Fitted=fitted(gam_Dino),
                         Residuals=resid(gam_Dino))

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- gam_Dino$model
# then we add the values of fitted and residuals
qq_data <- bind_cols(qq_data, ModelOutputs)

# Reorder sampling sites
qq_data <- qq_data %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

# Plot: verify that true data matches (more or less) model fit
ggplot(qq_data)+
  geom_point(aes(x = Day, y = true_count), color = 'red') +
  geom_point(aes(x = Day, y = Fitted), color = 'blue') +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  theme_classic() +
  labs(y = "Dinophysis count", x = "Calendar day")

# It matches! great!

# Now for the "official" qq plot, with color of the points depending on site
# Color palette
# Blues for Northern Brittany: '#0A1635', '#2B4561'
# Terra cotta for Pas de Calais: 'sienna4', 'tan3'
pheno_palette16 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646',
                     '#642C3A', '#DEB1CC', '#FC4D6B', '#791D40')

# We need to reorder the factor 'Code_point_Libelle' so the sites appear in the
# order we want
qq_data_reordered <- qq_data %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre')) %>%
  # we ungroup the data frame to produce only 1 qq-plot
  ungroup()

# And (qq-)plot
qqplot_custom <- ggplot(qq_data_reordered) +
  stat_qq(aes(sample=Residuals, color = Code_point_Libelle), alpha = .7) +
  stat_qq_line(aes(sample=Residuals, color = Code_point_Libelle)) +
  facet_wrap(facets = c('Code_point_Libelle')) +
  scale_color_discrete(type = pheno_palette16, guide = 'none') +
  theme_classic() +
  labs(y="Sample Quantiles",x="Theoretical Quantiles")

qqplot_custom

## Save the plot
# ggsave('Plots/GAMs/Dinophysis/Diagnostic_plots/qqplot_Dino_allsites_GAM.tiff',
# dpi = 300, height = 164, width = 164, units = 'mm', compression = 'lzw')

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data_reordered)+
  geom_point(aes(x=Fitted,y=Residuals, color  =Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free') +
  scale_color_discrete(type = pheno_palette16, guide = 'none') +
  theme_classic() +
  labs(y="Residuals",x="Fitted Values")

RvFplot_custom

## Save the plot
# ggsave('Plots/GAMs/Dinophysis/Diagnostic_plots/RvFplot_Dino_allsites_GAM.tiff',
# dpi = 300, height = 164, width = 164, units = 'mm', compression = 'lzw')

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data_reordered, aes(x = Residuals, 
                                                fill = Code_point_Libelle))+
  geom_histogram(binwidth = 1)+
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free') +
  scale_fill_discrete(type = pheno_palette16, guide = 'none') +
  theme_classic() +
  labs(x='Residuals', y = 'Count')

HistRes_custom

## Save the plot
# ggsave('Plots/GAMs/Dinophysis/Diagnostic_plots/HistRes_Dino_allsites_GAM.tiff',
# dpi = 300, height = 164, width = 164, units = 'mm', compression = 'lzw')

#### Confidence intervals ####

# We can try to obtain better confidence intervals.
# This section is almost entirely taken from a post by Gavin Simpson on his
# blog 'From the bottom of the heap'
# Reference : 
# https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/

# We want to define better confidence intervals and plot some examples of
# model fits taken from the bayesian posterior distribution of the model
# to illustrate model variability

# First
# A function for generating random values from a multivariate normal
rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m*n), m, n))
}


# We extract a few things from the fitted gam

# Bayesian covariance matrix (unconditional = TRUE means we adjust for the 
# smoothing parameters being estimated rather than known values)
Vb <- vcov(gam_Dino)
# New data
newd <- with(Season_Dino, data.frame(Day = rep(seq(1, 365, length = 365),
                                               length(unique(Year))*
                                                 length(unique(Code_point_Libelle))
                                               )))

newd <- group_by(newd, Day) %>%
  arrange(Day, by_group = TRUE)

# Add the Year vector
newd$Year <- rep(unique(Season_Dino$Year), 
                 365*length(unique(Season_Dino$Code_point_Libelle)))

newd <- group_by(newd, Year, Day) %>%
  arrange(Year, by_group = TRUE)

# And the site vector
newd$Code_point_Libelle <- rep(unique(Season_Dino$Code_point_Libelle), 
                               365*length(unique(Season_Dino$Year)))

newd <- ungroup(newd) %>%
  group_by(Code_point_Libelle) %>%
  arrange(Code_point_Libelle)

# Prediction by the model on the ***link*** scale, on the new data
pred <- predict(gam_Dino, newd, se.fit = TRUE, type = 'link')
# Isolate standard error of the predicted fit
se.fit <- pred$se.fit

# Set the pseudo-random seed to make results reproducible (?)
set.seed(42)
# specify the number of simulations to generate
N <- 5000

# N draws from a multivariate normal distributed matrix with mean 0 (?)
BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)

# Remove
rm(Vb)

# Calculating a function now (?)
Cg <- predict(gam_Dino, newd, type = "lpmatrix")
simDev <- Cg %*% t(BUdiff)

# Remove on the way to free memory space
rm(Cg)
rm(BUdiff)

### CAREFUL ! The 'absdev' and 'masd' steps can take a lot of time, especially
# when N is large. Don't panic, it eventually finishes running after a few tens
# of minutes.

# Finding the absolute values of the standardized deviation from the true model (?)
absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
rm(simDev)

# Computing the maximum of the absolute standard deviation
masd <- apply(absDev, 2L, max)
rm(absDev)

# Finding the critical value used to scale the standard errors to yield 
# the simultaneous interval
crit <- quantile(masd, prob = 0.95, type = 8)

# Now that crit is calculated, remove the big masd table that takes all the 
# memory space
rm(masd)

# It now runs much quicker. Nice.

# Adjusting the variables in the prediction to take
pred <- transform(cbind(data.frame(pred), newd),
                  uprP = fit + (2 * se.fit),
                  lwrP = fit - (2 * se.fit),
                  # The simultaneous CI is based on the critical value
                  # calculated just above
                  uprS = fit + (crit * se.fit),
                  lwrS = fit - (crit * se.fit))

# Calculating response_pred -> the gam fitted on every year on the repsonse
# scale
response_pred <- pred %>%
  mutate(across(c('fit', 'lwrS', 'uprS', 'lwrP', 'uprP'), ~ ilink(.)))

# We save response_pred now because it contains the variability linked to the
# random effect of the Year:
# write.csv2(response_pred,
# 'Data/GAM_outputs/Dinophysis/All_sites/Dino_GAM_response_pred_multiyear.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')

# Constructing a CI based on maximum and minimum simultaneous interval for
# plotting
pred_plot <- pred %>%
  group_by(Day, Code_point_Libelle) %>%
  summarise(
    # median
    median.fit = median(fit),
    # Confidence intervals
    # Here we take minimum and maximum to encompass all possible years (17), as
    # the prediction varies among years
    # Simultaneous confidence interval (CI)
    lwrS = min(lwrS), uprS = max(uprS),
    # Point-wise CI
    lwrP = min(lwrP), uprP = max(uprP), 
    .groups = 'keep')

# Transform values so they are expressed in the response variable units
# (by applying the inverse link function)
response_pred_plot <- as.data.frame(pred_plot) %>%
  mutate(across(c('median.fit', 'lwrS', 'uprS', 'lwrP', 'uprP'), ~ ilink(.)))

# And save that shall we
# write.csv2(response_pred_plot, 
#            'Data/GAM_outputs/Dinophysis/All_sites/Dino_GAM_response_pred.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')

# You can check that the backtransformation by the inverse link function
# worked by comparing the pred_plot and response_pred_plot tables.
# The response_pred_plot values should be higher.

# Checking the confidence intervals
# minimum and maximum values
min(response_pred_plot$median.fit)
min(response_pred_plot$lwrS) # Not < 0, nice :)
min(response_pred_plot$lwrP) # Not < 0, nice :)
max(response_pred_plot$median.fit)
max(response_pred_plot$uprS)
max(response_pred_plot$uprP)
# We can see that the simultaneous interval is slightly wider than the
# point-wise one. This is normal.

#### Plot aesthetics ####

# Defining a beautiful color palette
pheno_palette16 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646',
                     '#642C3A', '#DEB1CC', '#FC4D6B', '#791D40')

# Import 'response_pred' if necessary
response_pred_plot <- read.csv2(
  'Data/GAM_outputs/Dinophysis/All_sites/Dino_GAM_response_pred.csv',
                                header = TRUE, fileEncoding = 'ISO-8859-1')

# And Season_Dino
Season_Dino <- read.csv2('Data/REPHY_outputs/Season_Dino.csv',
                         header = TRUE, fileEncoding = 'ISO-8859-1')

# Reordering the factor 'Code_point_Libelle' so that the sites appear in the
# plot in the desired order
response_pred_plot <- response_pred_plot %>%
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

# And do the exact same thing in the Season_Dino dataset
Season_Dino <- Season_Dino %>%
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

#### Cropping out data points for plotting ####
# Crop out highest data points of each site to have better view of the model
Season_Dino_crop <- Season_Dino %>%
  # Point 1 Boulogne
  # filter(ifelse(true_count > 2 & Code_point_Libelle == 'Point 1 Boulogne',
  #               # if condition met, drop the line
  #               FALSE,
  #               # else, keep the line
  #               TRUE)) %>%
  # Antifer
  filter(ifelse(true_count > 50 & Code_point_Libelle == 'Antifer ponton pétrolier',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Cabourg
  filter(ifelse(true_count > 25 & Code_point_Libelle == 'Cabourg',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Men er Roue
  filter(ifelse(true_count > 8 & Code_point_Libelle == 'Men er Roue',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Ouest Loscolo
  filter(ifelse(true_count > 30 & Code_point_Libelle == 'Ouest Loscolo',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Le Cornard
  filter(ifelse(true_count > 5 & Code_point_Libelle == 'Le Cornard',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Auger
  filter(ifelse(true_count > 10 & Code_point_Libelle == 'Auger',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Arcachon - Bouée 7
  filter(ifelse(true_count > 13 & Code_point_Libelle == 'Arcachon - Bouée 7',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Teychan
  filter(ifelse(true_count > 10 & Code_point_Libelle == 'Teychan bis',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Parc Leucate
  filter(ifelse(true_count > 17 & Code_point_Libelle == 'Parc Leucate 2',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Bouzigues
  filter(ifelse(true_count > 8 & Code_point_Libelle == 'Bouzigues (a)',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Sète mer
  filter(ifelse(true_count > 5 & Code_point_Libelle == 'Sète mer',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Diana
  filter(ifelse(true_count > 8 & Code_point_Libelle == 'Diana centre',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE))


#### Plot the GAM ####

# We'll multiply by 100 to get a y-value in cells per L, which is more
# readable than a cell count.
pred_plot_cells_per_L <- response_pred_plot %>%
  mutate(across(c('median.fit', 'uprS', 'lwrS'), ~ .*100
  ))

# New plot in cells per L

ggplot(pred_plot_cells_per_L, aes(x = Day, y = median.fit, 
                               color = Code_point_Libelle,
                               fill = Code_point_Libelle)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
  # geom_path(data = stackFits_response, mapping = aes(y = values_response, x = Day, group= ind),
  #           alpha = 0.1, colour = "grey20") +
  geom_path(lwd = 1) +
  # scale_y_continuous(limits = c(0,40)) +
  geom_point(data = Season_Dino_crop, aes(x = Day, y = Dinophysis_genus), 
             size = .8, alpha = .5) +
  labs(y = c(expression(paste("Dinophysis cells.L"^'-1'))),
       x = "Day of the year",
       title = NULL
  ) +
  # Add a "ghost point" to force the minimum y-axis range to 5
  geom_point(aes(x = 1, y = 250), color = 'transparent', fill = 'transparent',
             size = .8, alpha = .5) +
  
  # Facet wrap
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette16, guide = 'none') +
  scale_fill_discrete(type = pheno_palette16, guide = 'none') +
  theme_classic() +
  theme(strip.text.x = element_text(size = 7.2, colour = "black"))

# Saving plot
# ggsave('Plots/GAMs/Dinophysis/Fig3_gam_Dino_cells_per_L.tiff',
#        dpi = 300, height = 135, width = 164,
#        units = 'mm', compression = 'lzw')

#### Plotting the GAM on all different years ####

## The logic here is to see how well the GAM fits the data year by year
# Can we see some years where the phenology significantly deviates from
# what the GAM predicts?

# We'll do this for one sampling site: Ouest Loscolo

#--- Ouest Loscolo ---#

response_pred_plot <- read.csv2(
  'Data/GAM_outputs/Dinophysis/All_sites/Dino_GAM_response_pred_multiyear.csv', 
                                header = TRUE, fileEncoding = 'ISO-8859-1')

rpp_OL <- filter(response_pred_plot, Code_point_Libelle == 'Ouest Loscolo')
Season_Dino_OL <- filter(Season_Dino, Code_point_Libelle == 'Ouest Loscolo')

ggplot(Season_Dino_OL, aes(x = Day, y = true_count)) +
  geom_point(color = 'dodgerblue3', alpha = .8) +
  geom_line(data = rpp_OL, aes(x = Day, y = fit), color = 'dodgerblue3') +
  geom_ribbon(data = rpp_OL, aes(x = Day, y = fit,
                                 ymin = lwrS, ymax = uprS), 
              color = 'dodgerblue3', alpha = 0.2) +
  facet_wrap(facets = c('Year'), scales = 'free_y') +
  labs(x = 'Day of the year', y = 'Dinophysis cells counted') +
  theme_classic()

# The GAM deviates significantly from the true data in some years (notably 2010).
# This is because the random effect for the Year is calculated on all sampling
# sites. In other words, the random effect determines 'good' and 'bad' years
# for Dinophysis based on the entire dataset. This is something that we could
# improve by reformulating the model, at the expense of some computing time.

# Save this plot
# ggsave('Plots/GAMs/Dinophysis/FigS3_Dino_GAM_OL_several_years.tiff',
# height = 200, width = 300, dpi = 300, unit = 'mm', compression = 'lzw')

#### First derivative of GAM fit ####

# We want to get the Dinophysis variation rate as a function of the DOY. For this,
# we'll simply calculate the 1st derivative of the GAM fit.

## We would like to have the derivative calculated separately for each year!
# This will be especially useful for our random forest model.
# So, we're going to do it the old-fashioned way with the function calc_deriv()
# of the gcplyr package

# First, we get the gam fit data for all years
gam_Dino_multiyear <- read.csv2(
  'Data/GAM_outputs/Dinophysis/All_sites/Dino_GAM_response_pred_multiyear.csv', 
                                header = TRUE, fileEncoding = 'ISO-8859-1')

# We need to add a vector of x values that correspond to sequential days
# from day 1 of 2007 to day 365 of 2022
gam_Dino_multiyear <- gam_Dino_multiyear %>%
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
  # And voila, we have our vector of sequential numbers:
  mutate(Date_number = as.numeric(Date))

# Now we just have to apply the function calc_deriv to calculate the derivative
gam_Dino_multiyear_deriv <- gam_Dino_multiyear %>%
  mutate(deriv = calc_deriv(y = fit, x = Date_number, return = 'derivative',
                            # Different derivative for each site
                            subset_by = Code_point_Libelle)) %>%
  mutate(deriv_upperS = calc_deriv(y = uprS, x = Date_number, return = 'derivative',
                                   # Different derivative for each site
                                   subset_by = Code_point_Libelle)) %>%
  mutate(deriv_lowerS = calc_deriv(y = lwrS, x = Date_number, return = 'derivative',
                                   # Different derivative for each site
                                   subset_by = Code_point_Libelle))

# We still need to do a few things
# First, we need to discard all derivatives calculated on the 1st of january 
# and 31st of december because of discontinuities between years
gam_Dino_multiyear_deriv <- gam_Dino_multiyear_deriv %>%
  filter(Day != 1) %>%
  filter(Day != 365) %>%
  # Reorder the sites
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

# Let's plot that to check if everything's ok

# We will plot on only 1 year, so we group by site and Day and compute a median, 
# + the min and max of the simultaneous confidence interval
gam_Dino_multiyear_deriv_plot <- gam_Dino_multiyear_deriv %>%
  group_by(Code_point_Libelle, Day) %>%
  summarise(median.deriv = median(deriv), upper = max(deriv_upperS), 
            lower = min(deriv_lowerS), .groups = 'keep')

# Color palette definition if needed
pheno_palette16 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646',
                     '#642C3A', '#DEB1CC', '#FC4D6B', '#791D40')

# Plot command
ggplot(gam_Dino_multiyear_deriv_plot, aes(x = Day, y = median.deriv, 
                                                 color = Code_point_Libelle,
                                                 fill = Code_point_Libelle)) +
  # Confidence interval
  geom_ribbon(aes(ymin = lower, ymax = upper, color = Code_point_Libelle,
                  fill = Code_point_Libelle), alpha = 0.2) +
  # Derivative fit
  geom_path(lwd = 1) +
  # Draw a line at 0 to separate accumulation from loss
  geom_line(aes(x = Day, y = 0), color = 'grey10', linewidth = .7) +
  labs(y = "Dinophysis change rate (cells/10mL/d-1)",
       x = "Day of the year",
       title = "1st derivative of Dinophysis GAM"
  ) +
  # Add 2 "ghost points" to force the limits of the y-axis to [-0.1 ; 0.1]
  geom_point(aes(x = 1, y = 0.1), color = 'transparent', fill = 'transparent') +
  geom_point(aes(x = 1, y = -0.1), color = 'transparent', fill = 'transparent') +
  # Separate by site
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y',
             # 4 rows to highlight the latitudinal change in phenology
             nrow = 4) +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette16, guide = 'none') +
  scale_fill_discrete(type = pheno_palette16, guide = 'none') +
  theme_classic()

# Everything seems fine?

# Save that table of 1st derivative values (we save the complete table, not the
# summarised one)
# write.csv2(gam_Dino_multiyear_deriv,
# 'Data/GAM_outputs/Dinophysis/All_sites/Dino_GAM_derivative.csv',
#                       row.names = FALSE, fileEncoding = 'ISO-8859-1')


####-----------------------------------------------------------------------####
#### Dinophysis GAM - 4 sites ####

# This is the same model than previously, but run on less data: only sampling
# sites and years for which we have Mesodinium and alloxanthin data.

# /!\ CAUTION /!\ For convenience, we use the same notation than in the previous
# section. So be sure you saved everything you needed, because tables/plots will
# be overwritten. Names of saved files (incl. figures) are of course different.

#--- Model and basic model checks ---#

# Formulate the GAM
  
  # Factorise some variables for the GAM
  Season_Dino_factor <- Season_Dino_4sites %>%
  # Transform variables as factors
  mutate(Year = as_factor(Year)) %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle))


gam_Dino <- gam(data = Season_Dino_factor, 
                # Only a spline for the day of the year
                # The use of a cyclic basis spline helps to make ends meet at the
                # first and last days of the year
                # 'k = -1' allows the model to fix the 'best' number of basic
                # functions (= knots)
                formula = true_count~s(Day, bs = 'cc', k = -1,
                                       # separate each site
                                       by = Code_point_Libelle) 
                # We add the Year as a random effect. This will help to assess and
                # smooth the effects of interannual variability in the phenology
                + s(Year, bs = 're', k = -1),
                # Using a Poisson distribution for count data
                family = poisson(),
                # Restricted maximum likelihood estimation (recommended method)
                method = 'REML')

summary(gam_Dino)

# Create a new 'empty' dataset for storing model prediction (fit)
gam_Dino_newdata <- expand_grid(Day=seq(1, 365),
                                # We add a Year vector as it has become a factor
                                # of the model
                                Year=seq(min(Season_Dino_4sites$Year), 
                                    max(Season_Dino_4sites$Year)),
                                # We also add the site data
                                Code_point_Libelle = unique(Season_Dino_4sites$Code_point_Libelle))

## Get the inverse link function of the model
# With this function, we can transform the prediction we make on the link
#scale to the response scale
ilink <- gam_Dino$family$linkinv

## Predict : add fit and se.fit on the **link** scale
gam_Dino_newdata <- bind_cols(gam_Dino_newdata, 
                              setNames(as_tibble(
                                predict(gam_Dino, gam_Dino_newdata, 
                                        se.fit = TRUE, type = 'link',
                                        re.form = ~ 1|Year)[1:2]),
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
# But we can see that the confidence interval is ridiculously small around the
# min and max of the model fit --> we'll fix that later down the script

# Quick plot

# Get sampling sites in the correct order
gam_Dino_newdata <- gam_Dino_newdata %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

Season_Dino_factor <- Season_Dino_factor %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# Plot:
ggplot(gam_Dino_newdata, aes(x = Day, y = fit_resp, color = as_factor(Year)))+
  geom_line(linewidth = 1) +
  geom_point(data = Season_Dino_4sites, aes(x = Day, y = true_count, 
                                     color = as_factor(Year)), 
             shape = 21, alpha = .6) +
  facet_wrap(facets = 'Code_point_Libelle', scales = 'free_y') +
  theme_classic()

# Ok nice. So we see that 1) each sampling site has a unique model fit, that
# matches more or less with the data, which is great. 2) the 'Year' random 
# effect is uniform among all sites, that is normal given how we formulated the
# model.

#--- Checking the model with diagnostic plots ---#
ModelOutputs<-data.frame(Fitted=fitted(gam_Dino),
                         Residuals=resid(gam_Dino))

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- gam_Dino$model
# then we add the values of fitted and residuals
qq_data <- bind_cols(qq_data, ModelOutputs)

# Reorder sampling sites
qq_data <- qq_data %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# Plot: verify that true data matches (more or less) model fit
ggplot(qq_data)+
  geom_point(aes(x = Day, y = true_count), color = 'red') +
  geom_point(aes(x = Day, y = Fitted), color = 'blue') +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  theme_classic() +
  labs(y = "Dinophysis count", x = "Calendar day")

# It matches! great!

# Now for the "official" qq plot, with color of the points depending on site
# Color palette with only 4 colors
pheno_palette4 <- c('red3', 'orangered', '#2156A1', '#5995E3')

# We need to reorder the factor 'Code_point_Libelle' so the sites appear in the
# order we want
qq_data_reordered <- qq_data %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo')) %>%
  # we ungroup the data frame to produce only 1 qq-plot
  ungroup()

# And (qq-)plot
qqplot_custom <- ggplot(qq_data_reordered) +
  stat_qq(aes(sample=Residuals, color = Code_point_Libelle), alpha = .7) +
  stat_qq_line(aes(sample=Residuals, color = Code_point_Libelle)) +
  facet_wrap(facets = c('Code_point_Libelle')) +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  labs(y="Sample Quantiles",x="Theoretical Quantiles")

qqplot_custom

# Save the plot
# ggsave('Plots/GAMs/Dinophysis/Diagnostic_plots/qqplot_Dino_4sites_GAM.tiff',
# dpi = 300, height = 164, width = 164, units = 'mm', compression = 'lzw')

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data_reordered)+
  geom_point(aes(x=Fitted,y=Residuals, color  =Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free') +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  labs(y="Residuals",x="Fitted Values")

RvFplot_custom

# Save the plot
# ggsave('Plots/GAMs/Dinophysis/Diagnostic_plots/RvFplot_Dino_4sites_GAM.tiff',
# dpi = 300, height = 164, width = 164, units = 'mm', compression = 'lzw')

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data_reordered, aes(x = Residuals, 
                                                fill = Code_point_Libelle))+
  geom_histogram(binwidth = 1)+
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free') +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  labs(x='Residuals', y = 'Count')

HistRes_custom

# Save the plot
# ggsave('Plots/GAMs/Dinophysis/Diagnostic_plots/HistRes_Dino_4sites_GAM.tiff',
# dpi = 300, height = 164, width = 164, units = 'mm', compression = 'lzw')

#### Confidence intervals ####

# We can try to obtain better confidence intervals.
# This section is almost entirely taken from a post by Gavin Simpson on his
# blog 'From the bottom of the heap'
# Reference : 
# https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/

# We want to define better confidence intervals and plot some examples of
# model fits taken from the bayesian posterior distribution of the model
# to illustrate model variability

# First
# A function for generating random values from a multivariate normal
rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m*n), m, n))
}


# We extract a few things from the fitted gam

# Bayesian covariance matrix (unconditional = TRUE means we adjust for the 
# smoothing parameters being estimated rather than known values)
Vb <- vcov(gam_Dino)
# New data
newd <- with(Season_Dino_4sites, data.frame(Day = rep(seq(1, 365, length = 365),
                                               length(unique(Year))*
                                                 length(unique(Code_point_Libelle))
)))

newd <- group_by(newd, Day) %>%
  arrange(Day, by_group = TRUE)

# Add the Year vector
newd$Year <- rep(unique(Season_Dino_4sites$Year), 
                 365*length(unique(Season_Dino_4sites$Code_point_Libelle)))

newd <- group_by(newd, Year, Day) %>%
  arrange(Year, by_group = TRUE)

# And the site vector
newd$Code_point_Libelle <- rep(unique(Season_Dino_4sites$Code_point_Libelle), 
                               365*length(unique(Season_Dino_4sites$Year)))

newd <- ungroup(newd) %>%
  group_by(Code_point_Libelle) %>%
  arrange(Code_point_Libelle)

# Prediction by the model on the ***link*** scale, on the new data
pred <- predict(gam_Dino, newd, se.fit = TRUE, type = 'link')
# Isolate standard error of the predicted fit
se.fit <- pred$se.fit

# Set the pseudo-random seed to make results reproducible (?)
set.seed(42)
# specify the number of simulations to generate
N <- 5000

# N draws from a multivariate normal distributed matrix with mean 0 (?)
BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)

# Remove
rm(Vb)

# Calculating a function now (?)
Cg <- predict(gam_Dino, newd, type = "lpmatrix")
simDev <- Cg %*% t(BUdiff)

# Remove on the way to free memory space
rm(Cg)
rm(BUdiff)

### CAREFUL ! The 'absdev' and 'masd' steps can take a lot of time, especially
# when N is large. Don't panic, it eventually finishes running after a few tens
# of minutes.

# Finding the absolute values of the standardized deviation from the true model (?)
absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
rm(simDev)

# Computing the maximum of the absolute standard deviation
masd <- apply(absDev, 2L, max)
rm(absDev)

# Finding the critical value used to scale the standard errors to yield 
# the simultaneous interval
crit <- quantile(masd, prob = 0.95, type = 8)

# Now that crit is calculated, remove the big masd table that takes all the 
# memory space
rm(masd)

# It now runs much quicker. Nice.

# Adjusting the variables in the prediction to take
pred <- transform(cbind(data.frame(pred), newd),
                  uprP = fit + (2 * se.fit),
                  lwrP = fit - (2 * se.fit),
                  # The simultaneous CI is based on the critical value
                  # calculated just above
                  uprS = fit + (crit * se.fit),
                  lwrS = fit - (crit * se.fit))

# Calculating response_pred -> the gam fitted on every year on the repsonse
# scale
response_pred <- pred %>%
  mutate(across(c('fit', 'lwrS', 'uprS', 'lwrP', 'uprP'), ~ ilink(.)))

# We save response_pred now because it contains the variability linked to the
# random effect of the Year:
write.csv2(response_pred,
'Data/GAM_outputs/Dinophysis/4_sites/Dino_GAM_response_pred_multiyear_4sites.csv',
           row.names = FALSE, fileEncoding = 'ISO-8859-1')

# Constructing a CI based on maximum and minimum simultaneous interval for
# plotting
pred_plot <- pred %>%
  group_by(Day, Code_point_Libelle) %>%
  summarise(
    # median
    median.fit = median(fit),
    # Confidence intervals
    # Here we take minimum and maximum to encompass all possible years (17), as
    # the prediction varies among years
    # Simultaneous confidence interval (CI)
    lwrS = min(lwrS), uprS = max(uprS),
    # Point-wise CI
    lwrP = min(lwrP), uprP = max(uprP), 
    .groups = 'keep')

# Transform values so they are expressed in the response variable units
# (by applying the inverse link function)
response_pred_plot <- as.data.frame(pred_plot) %>%
  mutate(across(c('median.fit', 'lwrS', 'uprS', 'lwrP', 'uprP'), ~ ilink(.)))

# And save that shall we
write.csv2(response_pred_plot,
'Data/GAM_outputs/Dinophysis/4_sites/Dino_GAM_response_pred_4sites.csv',
           row.names = FALSE, fileEncoding = 'ISO-8859-1')

# You can check that the backtransformation by the inverse link function
# worked by comparing the pred_plot and response_pred_plot tables.
# The response_pred_plot values should be higher.

# Checking the confidence intervals
# minimum and maximum values
min(response_pred_plot$median.fit)
min(response_pred_plot$lwrS) # Not < 0, nice :)
min(response_pred_plot$lwrP) # Not < 0, nice :)
max(response_pred_plot$median.fit)
max(response_pred_plot$uprS)
max(response_pred_plot$uprP)
# We can see that the simultaneous interval is slightly wider than the
# point-wise one. This is normal.

#### Plot aesthetics ####

# Defining a beautiful color palette for only 4 sampling sites
pheno_palette4 <- c('red3', 'orangered', '#2156A1', '#5995E3')

# Import 'response_pred' if necessary
response_pred_plot <- read.csv2(
  'Data/GAM_outputs/Dinophysis/4_sites/Dino_GAM_response_pred_4sites.csv',
                                header = TRUE, fileEncoding = 'ISO-8859-1')

# And Season_Dino
Season_Dino <- read.csv2('Data/REPHY_outputs/Season_Dino.csv',
                         header = TRUE, fileEncoding = 'ISO-8859-1') %>%
  # Filter for sampling sites and period
  filter(Code_point_Libelle %in% c('Antifer ponton pétrolier', 'Cabourg',
                                   'Men er Roue', 'Ouest Loscolo')) %>%
  filter(Year >= 2016)

# Reordering the factor 'Code_point_Libelle' so that the sites appear in the
# plot in the desired order
response_pred_plot <- response_pred_plot %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# And do the exact same thing in the Season_Dino dataset
Season_Dino <- Season_Dino %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

#### Cropping out data points for plotting ####
# Crop out highest data points of each site to have better view of the model
Season_Dino_crop <- Season_Dino %>%
  # Antifer
  filter(ifelse(true_count > 50 & Code_point_Libelle == 'Antifer ponton pétrolier',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Cabourg
  filter(ifelse(true_count > 25 & Code_point_Libelle == 'Cabourg',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Men er Roue
  filter(ifelse(true_count > 8 & Code_point_Libelle == 'Men er Roue',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Ouest Loscolo
  filter(ifelse(true_count > 30 & Code_point_Libelle == 'Ouest Loscolo',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE))


#### Plot the GAM ####

# We'll multiply by 100 to get a y-value in cells per L, which is more
# readable than a cell count.
pred_plot_cells_per_L <- response_pred_plot %>%
  mutate(across(c('median.fit', 'uprS', 'lwrS'), ~ .*100
  ))

# New plot in cells per L

ggplot(pred_plot_cells_per_L, aes(x = Day, y = median.fit, 
                                  color = Code_point_Libelle,
                                  fill = Code_point_Libelle)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
  geom_path(lwd = 1) +
  geom_point(data = Season_Dino_crop, aes(x = Day, y = Dinophysis_genus), 
             size = .8, alpha = .5) +
  labs(y = c(expression(paste("Dinophysis cells.L"^'-1'))),
       x = "Day of the year",
       title = NULL
  ) +
  # Add a "ghost point" to force the minimum y-axis range to 5
  geom_point(aes(x = 1, y = 250), color = 'transparent', fill = 'transparent',
             size = .8, alpha = .5) +
  
  # Facet wrap
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  theme(strip.text.x = element_text(size = 7.2, colour = "black"))

# Looks nice?

#### First derivative of GAM fit ####

# We want to get the Dinophysis variation rate as a function of the DOY. For this,
# we'll simply calculate the 1st derivative of the GAM fit.

## We would like to have the derivative calculated separately for each year!
# This will be especially useful for our random forest model.
# So, we're going to do it the old-fashioned way with the function calc_deriv()
# of the gcplyr package

# First, we get the gam fit data for all years
gam_Dino_multiyear <- read.csv2(
  'Data/GAM_outputs/Dinophysis/4_sites/Dino_GAM_response_pred_multiyear_4sites.csv', 
  header = TRUE, fileEncoding = 'ISO-8859-1')

# We need to add a vector of x values that correspond to sequential days
# from day 1 of 2016 to day 365 of 2022
gam_Dino_multiyear <- gam_Dino_multiyear %>%
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
  # And voila, we have our vector of sequential numbers:
  mutate(Date_number = as.numeric(Date))

# Now we just have to apply the function calc_deriv to calculate the derivative
gam_Dino_multiyear_deriv <- gam_Dino_multiyear %>%
  mutate(deriv = calc_deriv(y = fit, x = Date_number, return = 'derivative',
                            # Different derivative for each site
                            subset_by = Code_point_Libelle)) %>%
  mutate(deriv_upperS = calc_deriv(y = uprS, x = Date_number, return = 'derivative',
                                   # Different derivative for each site
                                   subset_by = Code_point_Libelle)) %>%
  mutate(deriv_lowerS = calc_deriv(y = lwrS, x = Date_number, return = 'derivative',
                                   # Different derivative for each site
                                   subset_by = Code_point_Libelle))

# We still need to do a few things
# First, we need to discard all derivatives calculated on the 1st of january 
# and 31st of december because of discontinuities between years
gam_Dino_multiyear_deriv <- gam_Dino_multiyear_deriv %>%
  filter(Day != 1) %>%
  filter(Day != 365) %>%
  # Reorder the sites
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# Let's plot that to check if everything's ok

# We will plot on only 1 year, so we group by site and Day and compute a median, 
# + the min and max of the simultaneous confidence interval
gam_Dino_multiyear_deriv_plot <- gam_Dino_multiyear_deriv %>%
  group_by(Code_point_Libelle, Day) %>%
  summarise(median.deriv = median(deriv), upper = max(deriv_upperS), 
            lower = min(deriv_lowerS), .groups = 'keep')

# Color palette definition if needed
pheno_palette16 <- c('red3', 'orangered', '#2156A1', '#5995E3')

# Plot command
ggplot(gam_Dino_multiyear_deriv_plot, aes(x = Day, y = median.deriv, 
                                          color = Code_point_Libelle,
                                          fill = Code_point_Libelle)) +
  # Confidence interval
  geom_ribbon(aes(ymin = lower, ymax = upper, color = Code_point_Libelle,
                  fill = Code_point_Libelle), alpha = 0.2) +
  # Derivative fit
  geom_path(lwd = 1) +
  # Draw a line at 0 to separate accumulation from loss
  geom_line(aes(x = Day, y = 0), color = 'grey10', linewidth = .7) +
  labs(y = "Dinophysis change rate (cells/10mL/d-1)",
       x = "Day of the year",
       title = "1st derivative of Dinophysis GAM"
  ) +
  # Add 2 "ghost points" to force the limits of the y-axis to [-0.1 ; 0.1]
  geom_point(aes(x = 1, y = 0.1), color = 'transparent', fill = 'transparent') +
  geom_point(aes(x = 1, y = -0.1), color = 'transparent', fill = 'transparent') +
  # Separate by site
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y',
             # 4 rows to highlight the latitudinal change in phenology
             nrow = 4) +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic()

# Everything seems fine?

# Save that table of 1st derivative values (we save the complete table, not the
# summarised one)
# write.csv2(gam_Dino_multiyear_deriv,
# 'Data/GAM_outputs/Dinophysis/4_sites/Dino_GAM_derivative_4sites.csv',
#                       row.names = FALSE, fileEncoding = 'ISO-8859-1')


####-----------------------------------------------------------------------####
#### Mesodinium GAM - 4 sites ####

# Now we do the same model but for Mesodinium. Only on the 4 sampling sites and 
# 2016-2022 period for which we have trustworthy data.

#--- Model and basic model checks ---#

# Formulate the GAM

# Only 4 selected sampling sites and from 2016 on --> Normally this step is 
# redundant because Season_Meso.csv was created on the relevant data only.
Season_Meso <- Season_Meso %>%
filter(Code_point_Libelle %in% c('Antifer ponton pétrolier', 'Cabourg',
                                 'Men er Roue', 'Ouest Loscolo')) %>%
  filter(Year >= 2016) %>%
  # and replace counts that are non-integers (less than 1 Mesodinium cell
  # counted) with 0
  mutate(true_count = 
           ifelse(true_count < 1, 0, true_count)) %>%
  # If the count is too high, discard it, as we did for Dino
  filter(true_count < 500)

# We will use the Year as a random effect and for this we need it as a factor
Season_Meso_factor <- Season_Meso %>%
  mutate(Year = as_factor(Year)) %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle))


gam_Meso <- gam(data = Season_Meso_factor, 
                # Only a spline for the day of the year
                # The use of a cyclic basis spline helps to make ends meet at the
                # first and last days of the year
                # 'k = -1' allows the model to fix the 'best' number of basic
                # functions (= knots)
                formula = true_count~s(Day, bs = 'cc', k = -1,
                                       # separate each site
                                       by = Code_point_Libelle) 
                # We add the Year as a random effect. This will help to assess and
                # smooth the effects of interannual variability in the phenology
                + s(Year, bs = 're', k = -1),
                # Using a Poisson distribution for count data
                family = poisson(),
                # Restricted maximum likelihood estimation (recommended method)
                method = 'REML')

summary(gam_Meso)

# Create a new 'empty' dataset for storing model prediction (fit)
gam_Meso_newdata <- expand_grid(Day=seq(1, 365),
                                # We add a Year vector as it has become a factor
                                # of the model
                                Year=seq(min(Season_Meso$Year), 
                                         max(Season_Meso$Year)),
                                # We also add the site data
                                Code_point_Libelle = unique(Season_Meso$Code_point_Libelle))

## Get the inverse link function of the model
# With this function, we can transform the prediction we make on the link
#scale to the response scale
ilink <- gam_Meso$family$linkinv

## Predict : add fit and se.fit on the **link** scale
gam_Meso_newdata <- bind_cols(gam_Meso_newdata, 
                              setNames(as_tibble(
                                predict(gam_Meso, gam_Meso_newdata, 
                                        se.fit = TRUE, type = 'link',
                                        re.form = ~ 1|Year)[1:2]),
                                c('fit_link','se_link')))

## Create the 95% confidence interval (2*standard error fit) and backtransform 
# to response variable using the inverse link function
gam_Meso_newdata <- mutate(gam_Meso_newdata,
                           fit_resp  = ilink(fit_link),
                           right_upr = ilink(fit_link + (2 * se_link)),
                           right_lwr = ilink(fit_link - (2 * se_link)))
# Check the confidence interval. It should not extend below 0 (negative counts
# are impossible)
min(gam_Meso_newdata$right_lwr) # Nice :)
min(gam_Meso_newdata$fit_resp)
max(gam_Meso_newdata$right_upr) # Nice too
max(gam_Meso_newdata$fit_resp)
# But we can see that the confidence interval is ridiculously small around the
# min and max of the model fit --> we'll fix that later down the script

# Quick plot

# Get sampling sites in the correct order
gam_Meso_newdata <- gam_Meso_newdata %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

Season_Meso_factor <- Season_Meso_factor %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# Plot:
ggplot(gam_Meso_newdata, aes(x = Day, y = fit_resp, 
                             color = as_factor(Year)))+
  geom_line(linewidth = 1) +
  geom_point(data = Season_Meso, aes(x = Day, y = true_count, 
                                     color = as_factor(Year)), 
             shape = 21, alpha = .6) +
  facet_wrap(facets = 'Code_point_Libelle', scales = 'free_y') +
  theme_classic()

# Ok nice. So we see that 1) each sampling site has a unique model fit, that
# matches more or less with the data, which is great. 2) the 'Year' random 
# effect is uniform among all sites, that is normal given how we formulated the
# model.

#--- Checking the model with diagnostic plots ---#
ModelOutputs<-data.frame(Fitted=fitted(gam_Meso),
                         Residuals=resid(gam_Meso))

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- gam_Meso$model
# then we add the values of fitted and residuals
qq_data <- bind_cols(qq_data, ModelOutputs)

# Reorder sampling sites
qq_data <- qq_data %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# Plot: verify that true data matches (more or less) model fit
ggplot(qq_data)+
  geom_point(aes(x = Day, y = true_count), color = 'red') +
  geom_point(aes(x = Day, y = Fitted), color = 'blue') +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  theme_classic() +
  labs(y = "Mesodinium count", x = "Calendar day")

# It' matches! great!'s not ideal to be honest...

# Now for the "official" qq plot, with color of the points depending on site
# Color palette with only 4 colors
pheno_palette4 <- c('red3', 'orangered', '#2156A1', '#5995E3')

# We need to reorder the factor 'Code_point_Libelle' so the sites appear in the
# order we want
qq_data_reordered <- qq_data %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo')) %>%
  # we ungroup the data frame to produce only 1 qq-plot
  ungroup()

# And (qq-)plot
qqplot_custom <- ggplot(qq_data_reordered) +
  stat_qq(aes(sample=Residuals, color = Code_point_Libelle), alpha = .7) +
  stat_qq_line(aes(sample=Residuals, color = Code_point_Libelle)) +
  facet_wrap(facets = c('Code_point_Libelle')) +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  labs(y="Sample Quantiles",x="Theoretical Quantiles")

qqplot_custom

# Save the plot
ggsave('Plots/GAMs/Mesodinium/Diagnostic_plots/qqplot_Meso_GAM.tiff',
dpi = 300, height = 164, width = 164, units = 'mm', compression = 'lzw')

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data_reordered)+
  geom_point(aes(x=Fitted,y=Residuals, color  =Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free') +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  labs(y="Residuals",x="Fitted Values")

RvFplot_custom

# Save the plot
ggsave('Plots/GAMs/Mesodinium/Diagnostic_plots/RvFplot_Meso_GAM.tiff',
dpi = 300, height = 164, width = 164, units = 'mm', compression = 'lzw')

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data_reordered, aes(x = Residuals, 
                                                fill = Code_point_Libelle))+
  geom_histogram(binwidth = 1)+
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free') +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  labs(x='Residuals', y = 'Count')

HistRes_custom

# Save the plot
ggsave('Plots/GAMs/Mesodinium/Diagnostic_plots/HistRes_Meso_GAM.tiff',
dpi = 300, height = 164, width = 164, units = 'mm', compression = 'lzw')

#### Confidence intervals ####

# We can try to obtain better confidence intervals.
# This section is almost entirely taken from a post by Gavin Simpson on his
# blog 'From the bottom of the heap'
# Reference : 
# https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/

# We want to define better confidence intervals and plot some examples of
# model fits taken from the bayesian posterior distribution of the model
# to illustrate model variability

# First
# A function for generating random values from a multivariate normal
rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m*n), m, n))
}


# We extract a few things from the fitted gam

# Bayesian covariance matrix (unconditional = TRUE means we adjust for the 
# smoothing parameters being estimated rather than known values)
Vb <- vcov(gam_Meso)
# New data
newd <- with(Season_Meso, data.frame(Day = rep(seq(1, 365, length = 365),
                                               length(unique(Year))*
                                                 length(unique(Code_point_Libelle))
)))

newd <- group_by(newd, Day) %>%
  arrange(Day, by_group = TRUE)

# Add the Year vector
newd$Year <- rep(unique(Season_Meso$Year), 
                 365*length(unique(Season_Meso$Code_point_Libelle)))

newd <- group_by(newd, Year, Day) %>%
  arrange(Year, by_group = TRUE)

# And the site vector
newd$Code_point_Libelle <- rep(unique(Season_Meso$Code_point_Libelle), 
                               365*length(unique(Season_Meso$Year)))

newd <- ungroup(newd) %>%
  group_by(Code_point_Libelle) %>%
  arrange(Code_point_Libelle)

# Prediction by the model on the ***link*** scale, on the new data
pred <- predict(gam_Meso, newd, se.fit = TRUE, type = 'link')
# Isolate standard error of the predicted fit
se.fit <- pred$se.fit

# Set the pseudo-random seed to make results reproducible (?)
set.seed(42)
# specify the number of simulations to generate
N <- 5000

# N draws from a multivariate normal distributed matrix with mean 0 (?)
BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)

# Remove
rm(Vb)

# Calculating a function now (?)
Cg <- predict(gam_Meso, newd, type = "lpmatrix")
simDev <- Cg %*% t(BUdiff)

# Remove on the way to free memory space
rm(Cg)
rm(BUdiff)

### CAREFUL ! The 'absdev' and 'masd' steps can take a lot of time, especially
# when N is large. Don't panic, it eventually finishes running after a few tens
# of minutes.

# Finding the absolute values of the standardized deviation from the true model (?)
absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
rm(simDev)

# Computing the maximum of the absolute standard deviation
masd <- apply(absDev, 2L, max)
rm(absDev)

# Finding the critical value used to scale the standard errors to yield 
# the simultaneous interval
crit <- quantile(masd, prob = 0.95, type = 8)

# Now that crit is calculated, remove the big masd table that takes all the 
# memory space
rm(masd)

# It now runs much quicker. Nice.

# Adjusting the variables in the prediction to take
pred <- transform(cbind(data.frame(pred), newd),
                  uprP = fit + (2 * se.fit),
                  lwrP = fit - (2 * se.fit),
                  # The simultaneous CI is based on the critical value
                  # calculated just above
                  uprS = fit + (crit * se.fit),
                  lwrS = fit - (crit * se.fit))

# Calculating response_pred -> the gam fitted on every year on the repsonse
# scale
response_pred <- pred %>%
  mutate(across(c('fit', 'lwrS', 'uprS', 'lwrP', 'uprP'), ~ ilink(.)))

# We save response_pred now because it contains the variability linked to the
# random effect of the Year:
write.csv2(response_pred,
'Data/GAM_outputs/Mesodinium/Meso_GAM_response_pred_multiyear.csv',
           row.names = FALSE, fileEncoding = 'ISO-8859-1')

# Constructing a CI based on maximum and minimum simultaneous interval for
# plotting
pred_plot <- pred %>%
  group_by(Day, Code_point_Libelle) %>%
  summarise(
    # median
    median.fit = median(fit),
    # Confidence intervals
    # Here we take minimum and maximum to encompass all possible years (17), as
    # the prediction varies among years
    # Simultaneous confidence interval (CI)
    lwrS = min(lwrS), uprS = max(uprS),
    # Point-wise CI
    lwrP = min(lwrP), uprP = max(uprP), 
    .groups = 'keep')

# Transform values so they are expressed in the response variable units
# (by applying the inverse link function)
response_pred_plot <- as.data.frame(pred_plot) %>%
  mutate(across(c('median.fit', 'lwrS', 'uprS', 'lwrP', 'uprP'), ~ ilink(.)))

# And save that shall we
write.csv2(response_pred_plot,
'Data/GAM_outputs/Mesodinium/Meso_GAM_response_pred.csv',
           row.names = FALSE, fileEncoding = 'ISO-8859-1')

# You can check that the backtransformation by the inverse link function
# worked by comparing the pred_plot and response_pred_plot tables.
# The response_pred_plot values should be higher.

# Checking the confidence intervals
# minimum and maximum values
min(response_pred_plot$median.fit)
min(response_pred_plot$lwrS) # Not < 0, nice :)
min(response_pred_plot$lwrP) # Not < 0, nice :)
max(response_pred_plot$median.fit)
max(response_pred_plot$uprS)
max(response_pred_plot$uprP)
# We can see that the simultaneous interval is slightly wider than the
# point-wise one. This is normal.

#### Plot aesthetics ####

# Defining a beautiful color palette for only 4 sampling sites
pheno_palette4 <- c('red3', 'orangered', '#2156A1', '#5995E3')

# Import 'response_pred' if necessary
response_pred_plot <- read.csv2(
  'Data/GAM_outputs/Mesodinium/Meso_GAM_response_pred.csv',
  header = TRUE, fileEncoding = 'ISO-8859-1')

# And Season_Meso
Season_Meso <- read.csv2('Data/REPHY_outputs/Season_Meso.csv',
                         header = TRUE, fileEncoding = 'ISO-8859-1') %>%
  # Filter for sampling sites and period
  filter(Code_point_Libelle %in% c('Antifer ponton pétrolier', 'Cabourg',
                                   'Men er Roue', 'Ouest Loscolo')) %>%
  filter(Year >= 2016)

# Reordering the factor 'Code_point_Libelle' so that the sites appear in the
# plot in the desired order
response_pred_plot <- response_pred_plot %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# And do the exact same thing in the Season_Meso dataset
Season_Meso <- Season_Meso %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

#### Cropping out data points for plotting ####
# Crop out highest data points of each site to have better view of the model
Season_Meso_crop <- Season_Meso %>%
  # Antifer
  filter(ifelse(true_count > 25 & Code_point_Libelle == 'Antifer ponton pétrolier',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Cabourg
  filter(ifelse(true_count > 90 & Code_point_Libelle == 'Cabourg',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Men er Roue
  filter(ifelse(true_count > 200 & Code_point_Libelle == 'Men er Roue',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Ouest Loscolo
  filter(ifelse(true_count > 200 & Code_point_Libelle == 'Ouest Loscolo',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE))


#### Plot the GAM ####

# We'll multiply by 100 to get a y-value in cells per L, which is more
# readable than a cell count.
pred_plot_cells_per_L <- response_pred_plot %>%
  mutate(across(c('median.fit', 'uprS', 'lwrS'), ~ .*100
  ))

# New plot in cells per L

ggplot(pred_plot_cells_per_L, aes(x = Day, y = median.fit, 
                                  color = Code_point_Libelle,
                                  fill = Code_point_Libelle)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
  geom_path(lwd = 1) +
  geom_point(data = Season_Meso_crop, aes(x = Day, y = Mesodinium_genus), 
             size = .8, alpha = .5) +
  labs(y = c(expression(paste("Mesodinium cells.L"^'-1'))),
       x = "Day of the year",
       title = NULL
  ) +
  # Add a "ghost point" to force the minimum y-axis range to 5
  geom_point(aes(x = 1, y = 250), color = 'transparent', fill = 'transparent',
             size = .8, alpha = .5) +
  
  # Facet wrap
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  theme(strip.text.x = element_text(size = 7.2, colour = "black"))

# Looks nice?

####-----------------------------------------------------------------------####
#### Alloxanthin GAM - 4 sites ####

# Now we do the same model but for the alloxanthin concentration. 
# Only on the 4 sampling sites and 2016-2022 period for which we have data.

# Because the alloxanthin concentration data doesn't follow a Poisson 
# distribution but a Gamma distribution (close to a log-normal distribution)

#--- Model and basic model checks ---#

# Formulate the GAM

# Clean dataset
Season_Allo <- Season_Allo %>%
  # Only 4 selected sampling sites and from 2016 on --> Normally this step is 
  # redundant because Season_Allo.csv was created on the relevant data only.
  filter(Code_point_Libelle %in% c('Antifer ponton pétrolier', 'Cabourg',
                                   'Men er Roue', 'Ouest Loscolo')) %>%
  filter(Year >= 2016) %>%
  # Filter 4 dates for which the alloxanthin concentration is = 0 (unlikely)
  filter(Allo != 0)

# Factorise Year and sampling site
Season_Allo_factor <- Season_Allo %>%
  mutate(Year = as_factor(Year)) %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle))

gam_Allo <- gam(data = Season_Allo_factor, 
                # Only a spline for the day of the year
                # The use of a cyclic basis spline helps to make ends meet at the
                # first and last days of the year
                # 'k = -1' allows the model to fix the 'best' number of basic
                # functions (= knots)
                formula = Allo~s(Day, bs = 'cc', k = -1,
                                       # separate each site
                                       by = Code_point_Libelle) 
                # We add the Year as a random effect. This will help to assess and
                # smooth the effects of interannual variability in the phenology
                + s(Year, bs = 're', k = -1),
                # Using a Gamma distribution (close to a log-normal distribution)
                family = Gamma(),
                # Restricted maximum likelihood estimation (recommended method)
                method = 'REML')

summary(gam_Allo)

# Create a new 'empty' dataset for storing model prediction (fit)
gam_Allo_newdata <- expand_grid(Day=seq(1, 365),
                                # We add a Year vector as it has become a factor
                                # of the model
                                Year=seq(min(Season_Allo$Year), 
                                         max(Season_Allo$Year)),
                                # We also add the site data
                                Code_point_Libelle = unique(Season_Allo$Code_point_Libelle))

## Get the inverse link function of the model
# With this function, we can transform the prediction we make on the link
#scale to the response scale
ilink <- gam_Allo$family$linkinv

## Predict : add fit and se.fit on the **link** scale
gam_Allo_newdata <- bind_cols(gam_Allo_newdata, 
                              setNames(as_tibble(
                                predict(gam_Allo, gam_Allo_newdata, 
                                        se.fit = TRUE, type = 'link',
                                        re.form = ~ 1|Year)[1:2]),
                                c('fit_link','se_link')))

## Create the 95% confidence interval (2*standard error fit) and backtransform 
# to response variable using the inverse link function
gam_Allo_newdata <- mutate(gam_Allo_newdata,
                           fit_resp  = ilink(fit_link),
                           right_upr = ilink(fit_link + (2 * se_link)),
                           right_lwr = ilink(fit_link - (2 * se_link)))
# Check the confidence interval. It should not extend below 0 (negative counts
# are impossible)
min(gam_Allo_newdata$right_lwr) # Nice :)
min(gam_Allo_newdata$fit_resp)
max(gam_Allo_newdata$right_upr) # Nice too
max(gam_Allo_newdata$fit_resp)
# But we can see that the confidence interval is ridiculously small around the
# min and max of the model fit --> we'll fix that later down the script

# Quick plot

# Get sampling sites in the correct order
gam_Allo_newdata <- gam_Allo_newdata %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

Season_Allo_factor <- Season_Allo_factor %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# Plot:
ggplot(gam_Allo_newdata, aes(x = Day, y = fit_resp, 
                             color = as_factor(Year)))+
  geom_line(linewidth = 1) +
  geom_point(data = Season_Allo, aes(x = Day, y = Allo, 
                                     color = as_factor(Year)), 
             shape = 21, alpha = .6) +
  facet_wrap(facets = 'Code_point_Libelle', scales = 'free_y') +
  theme_classic()

# Ok nice. So we see that 1) each sampling site has a unique model fit, that
# matches more or less with the data, which is great. 2) the 'Year' random 
# effect is uniform among all sites, that is normal given how we formulated the
# model.

#--- Checking the model with diagnostic plots ---#
ModelOutputs<-data.frame(Fitted=fitted(gam_Allo),
                         Residuals=resid(gam_Allo))

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- gam_Allo$model
# then we add the values of fitted and residuals
qq_data <- bind_cols(qq_data, ModelOutputs)

# Reorder sampling sites
qq_data <- qq_data %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# Plot: verify that true data matches (more or less) model fit
ggplot(qq_data)+
  geom_point(aes(x = Day, y = Allo), color = 'red') +
  geom_point(aes(x = Day, y = Fitted), color = 'blue') +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  theme_classic() +
  labs(y = "Alloxanthin count", x = "Calendar day")

# It matches! great!

# Now for the "official" qq plot, with color of the points depending on site
# Color palette with only 4 colors
pheno_palette4 <- c('red3', 'orangered', '#2156A1', '#5995E3')

# We need to reorder the factor 'Code_point_Libelle' so the sites appear in the
# order we want
qq_data_reordered <- qq_data %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo')) %>%
  # we ungroup the data frame to produce only 1 qq-plot
  ungroup()

# And (qq-)plot
qqplot_custom <- ggplot(qq_data_reordered) +
  stat_qq(aes(sample=Residuals, color = Code_point_Libelle), alpha = .7) +
  stat_qq_line(aes(sample=Residuals, color = Code_point_Libelle)) +
  facet_wrap(facets = c('Code_point_Libelle')) +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  labs(y="Sample Quantiles",x="Theoretical Quantiles")

qqplot_custom

# Save the plot
# ggsave('Plots/GAMs/Alloxanthin/Diagnostic_plots/qqplot_Allo_GAM.tiff',
# dpi = 300, height = 164, width = 164, units = 'mm', compression = 'lzw')

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data_reordered)+
  geom_point(aes(x=Fitted,y=Residuals, color  =Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free') +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  labs(y="Residuals",x="Fitted Values")

RvFplot_custom

# Save the plot
# ggsave('Plots/GAMs/Alloxanthin/Diagnostic_plots/RvFplot_Allo_GAM.tiff',
# dpi = 300, height = 164, width = 164, units = 'mm', compression = 'lzw')

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data_reordered, aes(x = Residuals, 
                                                fill = Code_point_Libelle))+
  geom_histogram(binwidth = .2)+
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free') +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  labs(x='Residuals', y = 'Count')

HistRes_custom

# Save the plot
# ggsave('Plots/GAMs/Alloxanthin/Diagnostic_plots/HistRes_Allo_GAM.tiff',
# dpi = 300, height = 164, width = 164, units = 'mm', compression = 'lzw')

#### Confidence intervals ####

# We can try to obtain better confidence intervals.
# This section is almost entirely taken from a post by Gavin Simpson on his
# blog 'From the bottom of the heap'
# Reference : 
# https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/

# We want to define better confidence intervals and plot some examples of
# model fits taken from the bayesian posterior distribution of the model
# to illustrate model variability

# First
# A function for generating random values from a multivariate normal
rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m*n), m, n))
}


# We extract a few things from the fitted gam

# Bayesian covariance matrix (unconditional = TRUE means we adjust for the 
# smoothing parameters being estimated rather than known values)
Vb <- vcov(gam_Allo)
# New data
newd <- with(Season_Allo, data.frame(Day = rep(seq(1, 365, length = 365),
                                               length(unique(Year))*
                                                 length(unique(Code_point_Libelle))
)))

newd <- group_by(newd, Day) %>%
  arrange(Day, by_group = TRUE)

# Add the Year vector
newd$Year <- rep(unique(Season_Allo$Year), 
                 365*length(unique(Season_Allo$Code_point_Libelle)))

newd <- group_by(newd, Year, Day) %>%
  arrange(Year, by_group = TRUE)

# And the site vector
newd$Code_point_Libelle <- rep(unique(Season_Allo$Code_point_Libelle), 
                               365*length(unique(Season_Allo$Year)))

newd <- ungroup(newd) %>%
  group_by(Code_point_Libelle) %>%
  arrange(Code_point_Libelle)

# Prediction by the model on the ***link*** scale, on the new data
pred <- predict(gam_Allo, newd, se.fit = TRUE, type = 'link')
# Isolate standard error of the predicted fit
se.fit <- pred$se.fit

# Set the pseudo-random seed to make results reproducible (?)
set.seed(42)
# specify the number of simulations to generate
N <- 5000

# N draws from a multivariate normal distributed matrix with mean 0 (?)
BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)

# Remove
rm(Vb)

# Calculating a function now (?)
Cg <- predict(gam_Allo, newd, type = "lpmatrix")
simDev <- Cg %*% t(BUdiff)

# Remove on the way to free memory space
rm(Cg)
rm(BUdiff)

### CAREFUL ! The 'absdev' and 'masd' steps can take a lot of time, especially
# when N is large. Don't panic, it eventually finishes running after a few tens
# of minutes.

# Finding the absolute values of the standardized deviation from the true model (?)
absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
rm(simDev)

# Computing the maximum of the absolute standard deviation
masd <- apply(absDev, 2L, max)
rm(absDev)

# Finding the critical value used to scale the standard errors to yield 
# the simultaneous interval
crit <- quantile(masd, prob = 0.95, type = 8)

# Now that crit is calculated, remove the big masd table that takes all the 
# memory space
rm(masd)

# It now runs much quicker. Nice.

# Adjusting the variables in the prediction to take
pred <- transform(cbind(data.frame(pred), newd),
                  uprP = fit + (2 * se.fit),
                  lwrP = fit - (2 * se.fit),
                  # The simultaneous CI is based on the critical value
                  # calculated just above
                  uprS = fit + (crit * se.fit),
                  lwrS = fit - (crit * se.fit))
# The simultaneous interval goes crazy here, with negative values. We'll have to
# stick with the point-wise interval

# Calculating response_pred -> the gam fitted on every year on the repsonse
# scale
response_pred <- pred %>%
  mutate(across(c('fit', 'lwrS', 'uprS', 'lwrP', 'uprP'), ~ ilink(.)))

# We save response_pred now because it contains the variability linked to the
# random effect of the Year:
# write.csv2(response_pred,
# 'Data/GAM_outputs/Alloxanthin/Allo_GAM_response_pred_multiyear.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')

# Constructing a CI based on maximum and minimum simultaneous interval for
# plotting (or point-wise cI in this case, see below)
pred_plot <- pred %>%
  group_by(Day, Code_point_Libelle) %>%
  summarise(
    # median
    median.fit = median(fit),
    # Confidence intervals
    # Here we take minimum and maximum to encompass all possible years (16), as
    # the prediction varies among years
    # Simultaneous confidence interval (CI)
    lwrS = min(lwrS), uprS = max(uprS),
    # Point-wise CI
    lwrP = min(lwrP), uprP = max(uprP), 
    .groups = 'keep')

# Transform values so they are expressed in the response variable units
# (by applying the inverse link function)
response_pred_plot <- as.data.frame(pred_plot) %>%
  mutate(across(c('median.fit', 'lwrS', 'uprS', 'lwrP', 'uprP'), ~ ilink(.)))

# And save that shall we
# write.csv2(response_pred_plot,
# 'Data/GAM_outputs/Alloxanthin/Allo_GAM_response_pred.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')

# Checking the confidence intervals
# minimum and maximum values
min(response_pred_plot$median.fit)
min(response_pred_plot$lwrS) # < 0, not nice :(
min(response_pred_plot$lwrP) # Not < 0, nice :)
max(response_pred_plot$median.fit)
max(response_pred_plot$uprS)
max(response_pred_plot$uprP)
# We can see that the simultaneous interval goes bonkers here. 
# We'll use the point-wise interval.

#### Plot aesthetics ####

# Defining a beautiful color palette for only 4 sampling sites
pheno_palette4 <- c('red3', 'orangered', '#2156A1', '#5995E3')

# Import 'response_pred' if necessary
response_pred_plot <- read.csv2(
  'Data/GAM_outputs/Alloxanthin/Allo_GAM_response_pred.csv',
  header = TRUE, fileEncoding = 'ISO-8859-1')

# And Season_Allo
Season_Allo <- read.csv2('Data/REPHY_outputs/Season_Allo.csv',
                         header = TRUE, fileEncoding = 'ISO-8859-1') %>%
  # Filter for sampling sites and period --> this shouldn't do anything but eh
  filter(Code_point_Libelle %in% c('Antifer ponton pétrolier', 'Cabourg',
                                   'Men er Roue', 'Ouest Loscolo')) %>%
  filter(Year >= 2016)

# Reordering the factor 'Code_point_Libelle' so that the sites appear in the
# plot in the desired order
response_pred_plot <- response_pred_plot %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# And do the exact same thing in the Season_Allo dataset
Season_Allo <- Season_Allo %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

#### Plot the GAM ####

# Plot
ggplot(response_pred_plot, aes(x = Day, y = median.fit, 
                                  color = Code_point_Libelle,
                                  fill = Code_point_Libelle)) +
  geom_ribbon(aes(ymin = lwrP, ymax = uprP), alpha = 0.2) +
  geom_path(lwd = 1) +
  geom_point(data = Season_Allo, aes(x = Day, y = Allo), 
             size = .8, alpha = .5) +
  labs(y = c(expression(paste("[Alloxanthin] mg.L"^'-1'))),
       x = "Day of the year",
       title = NULL
  ) +
  # Facet wrap
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  theme(strip.text.x = element_text(size = 7.2, colour = "black"))

# Looks nice?
####-----------------------------------------------------------------------####
#### Plotting Dino, Meso and Allo GAM together --> Figures 6 and S11 ####

# Now that we have all the GAM fits saved in data files, we can use them to plot
# phenology GAMs of the 3 levels of the trophic chain in the same sampling site

# Of course, this will be restricted to the period and sampling sites for which
# we had data for the Mesodinium and alloxanthin GAMS:
# - 2016-2022
# - South Brittany and Seine Bay sampling sites (4 in total)

#### Import data ####

#--- Observation data ---####

# Get the true observation data for Dinophysis, Mesodinium and alloxanthin
Season_Dino <- read.csv2('Data/REPHY_outputs/Season_Dino.csv', header = TRUE, 
                         fileEncoding = 'ISO-8859-1') %>%
  # Filter data to fit relevant temporal extent.
  filter(Year >= 2016)

Season_Meso <- read.csv2('Data/REPHY_outputs/Season_Meso.csv', header = TRUE, 
                         fileEncoding = 'ISO-8859-1') %>%
  filter(true_count < 500)

Season_Allo <- read.csv2('Data/REPHY_outputs/Season_Allo.csv', header = TRUE, 
                         fileEncoding = 'ISO-8859-1')

# Separate sampling sites (we will plot sampling sites 1 by 1:
# Dinophysis
Season_Dino_Antifer <- filter(Season_Dino, 
                              Code_point_Libelle == 'Antifer ponton pétrolier')
Season_Dino_Cabourg <- filter(Season_Dino, Code_point_Libelle == 'Cabourg')
Season_Dino_MeR <- filter(Season_Dino, Code_point_Libelle == 'Men er Roue')
Season_Dino_OL <- filter(Season_Dino, Code_point_Libelle == 'Ouest Loscolo')

# Mesodinium
Season_Meso_Antifer <- filter(Season_Meso, 
                              Code_point_Libelle == 'Antifer ponton pétrolier')
Season_Meso_Cabourg <- filter(Season_Meso, Code_point_Libelle == 'Cabourg')
Season_Meso_MeR <- filter(Season_Meso, Code_point_Libelle == 'Men er Roue')
Season_Meso_OL <- filter(Season_Meso, Code_point_Libelle == 'Ouest Loscolo')

# Alloxanthin
Season_Allo_Antifer <- filter(Season_Allo, 
                              Code_point_Libelle == 'Antifer ponton pétrolier')
Season_Allo_Cabourg <- filter(Season_Allo, Code_point_Libelle == 'Cabourg')
Season_Allo_MeR <- filter(Season_Allo, Code_point_Libelle == 'Men er Roue')
Season_Allo_OL <- filter(Season_Allo, Code_point_Libelle == 'Ouest Loscolo')

#--- GAM data ---####

## Dinophysis GAM data ##
gam_Dino <- read.csv2(
  'Data/GAM_outputs/Dinophysis/4_sites/Dino_GAM_response_pred_multiyear_4sites.csv',
                      header = TRUE, fileEncoding = 'ISO-8859-1')

# We need to summarise on 1 year to get a median fit and boundaries for the CI
gam_Dino <- gam_Dino %>%
  group_by(Code_point_Libelle, Day) %>%
  summarise(median.fit = median(fit), uprS = max(uprS), lwrS = min(lwrS),
            .groups = 'keep')

# Separate sampling sites
gam_Dino_Antifer <- filter(gam_Dino, 
                           Code_point_Libelle == 'Antifer ponton pétrolier')
gam_Dino_Cabourg <- filter(gam_Dino, Code_point_Libelle == 'Cabourg')
gam_Dino_MeR <- filter(gam_Dino, Code_point_Libelle == 'Men er Roue')
gam_Dino_OL <- filter(gam_Dino, Code_point_Libelle == 'Ouest Loscolo')

## Mesodinium GAM data ##
gam_Meso <- read.csv2(
  'Data/GAM_outputs/Mesodinium/Meso_GAM_response_pred_multiyear.csv', 
  header = TRUE, fileEncoding = 'ISO-8859-1')

# Same than for GAM Dino: summarise
gam_Meso <- gam_Meso %>%
  group_by(Code_point_Libelle, Day) %>%
  summarise(median.fit = median(fit), uprS = max(uprS), lwrS = min(lwrS),
            .groups = 'keep')

# Separate sites
gam_Meso_Antifer <- filter(gam_Meso, 
                           Code_point_Libelle == 'Antifer ponton pétrolier')
gam_Meso_Cabourg <- filter(gam_Meso, Code_point_Libelle == 'Cabourg')
gam_Meso_MeR <- filter(gam_Meso, Code_point_Libelle == 'Men er Roue')
gam_Meso_OL <- filter(gam_Meso, Code_point_Libelle == 'Ouest Loscolo')

## Alloxanthin GAM data ##
gam_Allo <- read.csv2(
  'Data/GAM_outputs/Alloxanthin/Allo_GAM_response_pred.csv',
                           header = TRUE, fileEncoding = 'ISO-8859-1')

# Separate sites
gam_Allo_Antifer <- filter(gam_Allo, 
                           Code_point_Libelle == 'Antifer ponton pétrolier')
gam_Allo_Cabourg <- filter(gam_Allo, Code_point_Libelle == 'Cabourg')
gam_Allo_MeR <- filter(gam_Allo, Code_point_Libelle == 'Men er Roue')
gam_Allo_OL <- filter(gam_Allo, Code_point_Libelle == 'Ouest Loscolo')

#--- Satellite data ---####
# Ok nice. Now we'll try to add dots for satellite observations of Mesodinium
# red tides in the regions of interest.

# Loire-Vilaine (= South Brittany)
Satellite_survey_LV <- read.csv2(
  'Data/Satellite/Bloom_dates/Blooms_list_satellite_Vilaine-Gironde.csv',
                                 header = TRUE, fileEncoding = 'ISO-8859-1')

# Let's focus on Mesodinium in the Loire-Atlantique region
Satellite_survey_Meso_LV <- filter(Satellite_survey_LV,
                                   (grepl('Mesodinium', Main.species) |
                                      grepl('Mesodinium', Comment) |
                                      grepl('Dark red', Comment)) &
                                     Region == 'Loire Atlantique' &
                                     Satellite == 'S2') %>%
  # Add the date information that we will need for plotting
  mutate(Date = ymd(Date)) %>%
  mutate(Year = year(Date)) %>%
  mutate(Day = yday(Date)) %>%
  # Finally, add a false "true_count" value for plotting alongside the Dinophysis
  # phenology
  mutate(true_count = 50) %>%
  select(-c(7:11))

# Seine Bay
Satellite_survey_SB <- read.csv2(
  'Data/Satellite/Bloom_dates/Blooms_list_satellite_Seine_Bay_Meso.csv',
                                 header = TRUE, fileEncoding = 'ISO-8859-1') %>%
  # Add the date information that we will need for plotting
  mutate(Date = ymd(Date)) %>%
  mutate(Year = year(Date)) %>%
  mutate(Day = yday(Date)) %>%
  # Finally, add a false "true_count" value for plotting alongside the Dinophysis
  # phenology
  mutate(true_count = 50)

### Join both datasets
Satellite_survey_Meso <- bind_rows(Satellite_survey_Meso_LV, Satellite_survey_SB) %>%
  # Attribute a Code_point_Libelle value to facilitate plotting
  mutate(Code_point_Libelle = ifelse(Region == 'Loire Atlantique', 'Ouest Loscolo',
                                     'Cabourg'))

# Good

#### Plots ####

#Ok the objective here is to produce inidivdual plots for each GAM in each site,
# and combine them afterwards with ggarrange()

#--- Antifer ---####

# Alloxanthin
plot_Allo_Antifer <- ggplot() +
  # plot the data
  geom_point(data = Season_Allo_Antifer, aes(x = Day, y = Allo), size = 2,
             alpha = .3, color = 'red3') +
  # plot the GAM
  geom_ribbon(data = gam_Allo_Antifer, aes(x = Day, ymin = lwrP, 
                                           ymax = uprP),
              linewidth = .35, alpha = .2,
              color = 'red3', fill = 'red3') +
  geom_line(data = gam_Allo_Antifer, aes(x = Day, y = median.fit),
            linewidth = .7, color = 'red3') +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(# title = 'A. Antifer',
    subtitle = '', x = NULL, 
    y = NULL) +
  theme_classic()

plot_Allo_Antifer

# Mesodinium
plot_Meso_Antifer <- ggplot() +
  # plot the data as points
  geom_point(data = Season_Meso_Antifer, aes(x = Day, y = true_count*100), size = 2,
             alpha = .3, color = 'red3') +
  # plot the GAM
  geom_ribbon(data = gam_Meso_Antifer, aes(x = Day, ymin = lwrS*100, 
                                           ymax = uprS*100),
              linewidth = .35, alpha = .2,
              color = 'red3', fill = 'red3') +
  geom_line(data = gam_Meso_Antifer, aes(x = Day, y = median.fit*100),
            linewidth = .7, color = 'red3') +
  # Satellite observations of Mesodinium blooms
  geom_point(data = subset(Satellite_survey_Meso, Code_point_Libelle == 'Cabourg'),
             aes(x = Day, y = 2000),
             color = 'firebrick4', shape = 8, size = 3.5, stroke = .45) +
  # Text
  labs(subtitle = '', x = NULL, 
       y = NULL) +
  theme_classic()

plot_Meso_Antifer

# Dinophysis
plot_Dino_Antifer <- ggplot() +
  # plot the data as points
  geom_point(data = Season_Dino_Antifer, aes(x = Day, y = true_count*100), size = 2,
             alpha = .3, color = 'red3') +
  # plot the GAM
  geom_ribbon(data = gam_Dino_Antifer, aes(x = Day, ymin = lwrS*100, 
                                           ymax = uprS*100),
              linewidth = .35, alpha = .2,
              color = 'red3', fill = 'red3') +
  geom_line(data = gam_Dino_Antifer, aes(x = Day, y = median.fit*100),
            linewidth = .7, color = 'red3') +
  # Text
  labs(subtitle = '', x = 'Day of the year', y = NULL) +
  theme_classic()

plot_Dino_Antifer

# Arrange the plot
succession_Antifer <- ggarrange(plot_Allo_Antifer, plot_Meso_Antifer, plot_Dino_Antifer, nrow = 3,
                                align = 'v')

# Save the plot
# ggsave('Plots/GAMs/Successions/Succession_plot_Antifer_pub.tiff',
#        height = 130, width = 82, dpi = 300, unit = 'mm', compression = 'lzw')


#--- Men er Roue ---####

# Alloxanthin
plot_Allo_MeR <- ggplot() +
  # plot the data
  geom_point(data = Season_Allo_MeR, aes(x = Day, y = Allo), size = 2,
             alpha = .3, color = '#2156A1') +
  # plot the GAM
  geom_ribbon(data = gam_Allo_MeR, aes(x = Day, ymin = lwrP, 
                                       ymax = uprP),
              linewidth = .35, alpha = .2,
              color = '#2156A1', fill = '#2156A1') +
  geom_line(data = gam_Allo_MeR, aes(x = Day, y = median.fit),
            linewidth = .7, color = '#2156A1') +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(# title = 'B. Ouest LoscMeRo',
    subtitle = '',
    x = NULL, # subtitle = 'Alloxanthin (microgram.L-1)', 
    y = NULL) +
  theme_classic()

plot_Allo_MeR

# Mesodinium
plot_Meso_MeR <- ggplot() +
  # plot the data as points
  geom_point(data = Season_Meso_MeR, aes(x = Day, y = true_count*100), size = 2,
             alpha = .3, color = '#2156A1') +
  # plot the GAM
  geom_ribbon(data = gam_Meso_MeR, aes(x = Day, ymin = lwrS*100, 
                                       ymax = uprS*100),
              linewidth = .35, alpha = .2,
              color = '#2156A1', fill = '#2156A1') +
  geom_line(data = gam_Meso_MeR, aes(x = Day, y = median.fit*100),
            linewidth = .7, color = '#2156A1') +
  # Satellite observations of Mesodinium blooms
  geom_point(data = subset(Satellite_survey_Meso, Code_point_Libelle == 'Ouest Loscolo'),
             aes(x = Day, y = 25000),
             color = 'firebrick4', shape = 8, size = 3.5, stroke = .45) +
  # Text
  labs(subtitle = '',
       x = NULL, # subtitle = 'Mesodinium (cells counted)',
       y = NULL) +
  theme_classic()

plot_Meso_MeR

# Dinophysis
plot_Dino_MeR <- ggplot() +
  # plot the data as points
  geom_point(data = Season_Dino_MeR, aes(x = Day, y = true_count*100), size = 2,
             alpha = .3, color = '#2156A1') +
  # plot the GAM
  geom_ribbon(data = gam_Dino_MeR, aes(x = Day, ymin = lwrS*100, 
                                       ymax = uprS*100),
              linewidth = .35, alpha = .2,
              color = '#2156A1', fill = '#2156A1') +
  geom_line(data = gam_Dino_MeR, aes(x = Day, y = median.fit*100),
            linewidth = .7, color = '#2156A1') +
  # Text
  labs(subtitle = '', 
       x = 'Day of the year', y = NULL) + # subtitle = 'Dinophysis (cells counted)', 
  theme_classic()

plot_Dino_MeR

# Arrange the plot
succession_MeR <- ggarrange(plot_Allo_MeR, plot_Meso_MeR, plot_Dino_MeR, nrow = 3,
                            align = 'v')

# Save the plot
# ggsave('Plots/GAMs/Successions/Succession_plot_MeR_pub.tiff',
#        height = 130, width = 82, dpi = 300, unit = 'mm', compression = 'lzw')

#--- Ouest Loscolo ---####

# Alloxanthin
plot_Allo_OL <- ggplot() +
  # plot the data
  geom_point(data = Season_Allo_OL, aes(x = Day, y = Allo), size = 2,
             alpha = .3, color = '#5995E3') +
  # plot the GAM
  geom_ribbon(data = gam_Allo_OL, aes(x = Day, ymin = lwrP, 
                                      ymax = uprP),
              linewidth = .35, alpha = .2,
              color = '#5995E3', fill = '#5995E3') +
  geom_line(data = gam_Allo_OL, aes(x = Day, y = median.fit),
            linewidth = .7, color = '#5995E3') +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(# title = 'B. Ouest Loscolo',
    subtitle = '',
    x = NULL, # subtitle = 'Alloxanthin (microgram.L-1)', 
    y = NULL) +
  theme_classic()

plot_Allo_OL

# Mesodinium
plot_Meso_OL <- ggplot() +
  # plot the data as points
  geom_point(data = Season_Meso_OL, aes(x = Day, y = true_count*100), size = 2,
             alpha = .3, color = '#5995E3') +
  # plot the GAM
  geom_ribbon(data = gam_Meso_OL, aes(x = Day, ymin = lwrS*100, 
                                      ymax = uprS*100),
              linewidth = .35, alpha = .2,
              color = '#5995E3', fill = '#5995E3') +
  geom_line(data = gam_Meso_OL, aes(x = Day, y = median.fit*100),
            linewidth = .7, color = '#5995E3') +
  # Satellite observations of Mesodinium blooms
  geom_point(data = subset(Satellite_survey_Meso, Code_point_Libelle == 'Ouest Loscolo'),
             aes(x = Day, y = true_count*300),
             color = 'firebrick4', shape = 8, size = 3.5, stroke = .45) +
  # Text
  labs(subtitle = '',
       x = NULL, # subtitle = 'Mesodinium (cells counted)',
       y = NULL) +
  theme_classic()

plot_Meso_OL

# Dinophysis
plot_Dino_OL <- ggplot() +
  # plot the data as points
  geom_point(data = Season_Dino_OL, aes(x = Day, y = true_count*100), size = 2,
             alpha = .3, color = '#5995E3') +
  # plot the GAM
  geom_ribbon(data = gam_Dino_OL, aes(x = Day, ymin = lwrS*100, 
                                      ymax = uprS*100),
              linewidth = .35, alpha = .2,
              color = '#5995E3', fill = '#5995E3') +
  geom_line(data = gam_Dino_OL, aes(x = Day, y = median.fit*100),
            linewidth = .7, color = '#5995E3') +
  # Text
  labs(subtitle = '', 
       x = 'Day of the year', y = NULL) + # subtitle = 'Dinophysis (cells counted)', 
  theme_classic()

plot_Dino_OL

# Arrange the plot
succession_OL <- ggarrange(plot_Allo_OL, plot_Meso_OL, plot_Dino_OL, nrow = 3,
                           align = 'v')

# Save the plot
# ggsave('Plots/GAMs/Successions/Succession_plot_OL_pub.tiff',
#        height = 130, width = 82, dpi = 300, unit = 'mm', compression = 'lzw')

#--- Cabourg ---####

# Alloxanthin
plot_Allo_Cabourg <- ggplot() +
  # plot the data
  geom_point(data = Season_Allo_Cabourg, aes(x = Day, y = Allo), size = 2,
             alpha = .3, color = 'orangered') +
  # plot the GAM
  geom_ribbon(data = gam_Allo_Cabourg, aes(x = Day, ymin = lwrP, 
                                           ymax = uprP),
              linewidth = .35, alpha = .2,
              color = 'orangered', fill = 'orangered') +
  geom_line(data = gam_Allo_Cabourg, aes(x = Day, y = median.fit),
            linewidth = .7, color = 'orangered') +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(# title = 'A. Cabourg',
    subtitle = '', x = NULL, 
    y = NULL) +
  theme_classic()

plot_Allo_Cabourg

# Mesodinium
plot_Meso_Cabourg <- ggplot() +
  # plot the data as points
  geom_point(data = Season_Meso_Cabourg, aes(x = Day, y = true_count*100), size = 2,
             alpha = .3, color = 'orangered') +
  # plot the GAM
  geom_ribbon(data = gam_Meso_Cabourg, aes(x = Day, ymin = lwrS*100, 
                                           ymax = uprS*100),
              linewidth = .35, alpha = .2,
              color = 'orangered', fill = 'orangered') +
  geom_line(data = gam_Meso_Cabourg, aes(x = Day, y = median.fit*100),
            linewidth = .7, color = 'orangered') +
  # Satellite observations of Mesodinium blooms
  geom_point(data = subset(Satellite_survey_Meso, Code_point_Libelle == 'Cabourg'),
             aes(x = Day, y = true_count*300),
             color = 'firebrick4', shape = 8, size = 3.5, stroke = .45) +
  # Text
  labs(subtitle = '', x = NULL, 
       y = NULL) +
  theme_classic()

plot_Meso_Cabourg

# Dinophysis
plot_Dino_Cabourg <- ggplot() +
  # plot the data as points
  geom_point(data = Season_Dino_Cabourg, aes(x = Day, y = true_count*100), size = 2,
             alpha = .3, color = 'orangered') +
  # plot the GAM
  geom_ribbon(data = gam_Dino_Cabourg, aes(x = Day, ymin = lwrS*100, 
                                           ymax = uprS*100),
              linewidth = .35, alpha = .2,
              color = 'orangered', fill = 'orangered') +
  geom_line(data = gam_Dino_Cabourg, aes(x = Day, y = median.fit*100),
            linewidth = .7, color = 'orangered') +
  # Text
  labs(subtitle = '', x = 'Day of the year', y = NULL) +
  theme_classic()

plot_Dino_Cabourg

# Arrange the plot
succession_Cab <- ggarrange(plot_Allo_Cabourg, plot_Meso_Cabourg, plot_Dino_Cabourg, nrow = 3,
                            align = 'v')

# Save the plot
# ggsave('Plots/GAMs/Successions/Succession_plot_Cabourg_pub.tiff',
#        height = 130, width = 82, dpi = 300, unit = 'mm', compression = 'lzw')


## Combining plots ####
#--- Cabourg + Ouest Loscolo (Figure 6) ---####

ggarrange(succession_Cab, succession_OL, nrow = 1,
          align = 'h')

# Save that pretty stuff
# ggsave('Plots/GAMs/Successions/Fig6_Succession_plot_Cab_OL.tiff',
#        width = 154, height = 170, units = 'mm', dpi = 300, compression = 'lzw')

### Little hidden bonus: getting peak dates for the 3 GAMs
# Cabourg

max_allo_cab <- gam_Allo_Cabourg %>%
  filter(median.fit == max(gam_Allo_Cabourg$median.fit))

max_meso_cab <- gam_Meso_Cabourg %>%
  filter(median.fit == max(gam_Meso_Cabourg$median.fit))

max_dino_cab <- gam_Dino_Cabourg %>%
  filter(median.fit == max(gam_Dino_Cabourg$median.fit))

# Ouest Loscolo
# Because this one is "double-peak" style, we need to first isolate the first or
# second peak, depending on the one we want. Here we isolate the first peak, just
# need to reverse '<=' to '>=' to get the second.

gam_Allo_OL_1 <- gam_Allo_OL %>%
  filter(Day <= 200)
gam_Meso_OL_1 <- gam_Meso_OL %>%
  filter(Day <= 200)
gam_Dino_OL_1 <- gam_Dino_OL %>%
  filter(Day <= 225)

max_allo_OL <- gam_Allo_OL %>%
  filter(median.fit == max(gam_Allo_OL_1$median.fit))

max_meso_OL <- gam_Meso_OL %>%
  filter(median.fit == max(gam_Meso_OL_1$median.fit))

max_dino_OL <- gam_Dino_OL %>%
  filter(median.fit == max(gam_Dino_OL_1$median.fit))

#--- Antifer + Men er Roue (Figure S11) ---####

ggarrange(succession_Antifer, succession_MeR, nrow = 1,
          align = 'h')

# ggsave('Plots/GAMs/Successions/FigS11_Succession_plot_Antifer_MeR.tiff',
#        width = 154, height = 170, units = 'mm', dpi = 300, compression = 'lzw')

####-------------------------- End of script ------------------------------####