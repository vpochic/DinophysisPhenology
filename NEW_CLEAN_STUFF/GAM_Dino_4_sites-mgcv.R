###### GAM for Dinophysis phenology with mgcv - only 4 sites ###
## V. POCHIC
# 2025-11-07

# The purpose of this script is to create a GAM for Dinophysis on the 4 sites
# (Southern Brittany and Seine Bay), and on the years 2016-2022, to be used in
# Random Forest/GLMs along the other GAMs

#### Packages and functions ####
library(tidyverse)
library(ggplot2)
library(mgcv)
library(gratia)
library(viridis)
library(cmocean)

#### Import data ####

# We will import the Season_Dino that is defined in the "big GAM" script for
# Dinophysis ('GAM_Dino_all_sites-mgcv.R')

Season_Dino <- read.csv2('Data/REPHY_outputs/Season_Dino_20250604.csv',
                         header = TRUE, fileEncoding = 'ISO-8859-1')

Season_Dino <- Season_Dino %>%
  filter(Year >= 2016 & Code_point_Libelle %in% c(
    'Antifer ponton pétrolier', 'Cabourg', 'Men er Roue', 'Ouest Loscolo'
  )) %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# Alright

########## Running the GAM ###########
#### Model and basic model checks ####

# Formulate the GAM

# We will use the Year as a random effect and for this we need it as a factor
Season_Dino_factor <- Season_Dino %>%
  mutate(Year = as_factor(Year))


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
                # Introducing the weights
                # weights = unif_weight,
                # Using a Poisson distribution for count data
                family = poisson(),
                # Restricted maximum likelihood estimation (recommended method)
                method = 'REML')

summary(gam_Dino)
gam.check(gam_Dino)

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
# But we can see that the confidence interval is quite small around the
# min and max of the model fit.

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
  geom_point(aes(x = Day, y = true_count), color = 'red') +
  geom_point(aes(x = Day, y = Fitted), color = 'blue') +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  theme_classic() +
  labs(y = "Dinophysis count", x = "Calendar day")

# It matches! great! (A bit less great for Men er Roue. But this will do)

# Now for the "official" qq plot, with color of the points depending on site
# Color palette
# Blues for Northern Brittany: '#0A1635', '#2B4561'
# Terra cotta for Pas de Calais: 'sienna4', 'tan3'
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
  facet_wrap(facets = c('Code_point_Libelle'), nrow = 1) +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  labs(title = NULL,
       y="Sample Quantiles",x="Theoretical Quantiles")

qqplot_custom

# Save the plot
# ggsave('Plots/GAMs/Dinophysis/qqplot_custom_Dino_4sites.tiff', dpi = 300,
#        height = 55, width = 164, units = 'mm', compression = 'lzw')

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data_reordered)+
  geom_point(aes(x=Fitted,y=Residuals, color = Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free', nrow = 1) +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  labs(title = NULL,
       y="Residuals",x="Fitted Values")

RvFplot_custom

# Save the plot
# ggsave('Plots/GAMs/Dinophysis/Rvfplot_custom_Dino_4sites.tiff', dpi = 300, height = 55, width = 164,
#                units = 'mm', compression = 'lzw')

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data_reordered, aes(x = Residuals, 
                                                fill = Code_point_Libelle))+
  geom_histogram(binwidth = 1)+
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free', nrow = 1) +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  labs(title = NULL,
       x='Residuals', y = 'Count')

HistRes_custom

# Save the plot
# ggsave('Plots/GAMs/Dinophysis/HistRes_custom_Dino_4sites.tiff', dpi = 300,
#        height = 55, width = 164, units = 'mm', compression = 'lzw')

#### Confidence intervals ####

# We can try to obtain better confidence intervals.
# This section is almost entirely taken from a post by Gavin Simpson on his
# blog 'From the bottom of the heap' (big thank you Gavin Simpson!)
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

# Save response_pred
# write.csv2(response_pred, 'Data/GAM_outputs/Dinophysis/4_sites/response_pred_GAMDino_4sites_20250604.csv',
# row.names = FALSE, fileEncoding = 'ISO-8859-1')

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

# plotting
ggplot(response_pred_plot, aes(x = Day)) +
  geom_line(aes(y = median.fit), lwd = 1, color = 'black') +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2, fill = "red") +
  geom_ribbon(aes(ymin = lwrP, ymax = uprP), alpha = 0.2, fill = "red") +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  labs(y = "Dinophysis (cells.L-1)",
       x = "Day of the year")


###### Plot aesthetics #########

# Defining a beautiful color palette
pheno_palette4 <- c('red3', 'orangered', '#2156A1', '#5995E3')

# Reordering the factor 'Code_point_Libelle' so that the sites appear in the
# plot in the desired order
response_pred_plot <- response_pred_plot %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# And do the exact same thing in the Season_Dino dataset
Season_Dino <- Season_Dino %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

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


##### Plot this #####
ggplot(response_pred_plot, aes(x = Day, y = median.fit, 
                               color = Code_point_Libelle,
                               fill = Code_point_Libelle)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
  geom_path(lwd = 1) +
  # scale_y_continuous(limits = c(0,40)) +
  geom_point(data = Season_Dino_crop, aes(x = Day, y = true_count), 
             size = .8, alpha = .5) +
  # Add a "ghost point" to force the minimum y-axis range to 5
  geom_point(aes(x = 1, y = 5), color = 'transparent', fill = 'transparent',
             size = .8, alpha = .5) +
  labs(y = "Dinophysis cells observed",
       x = "Julian day",
       title = "Poisson GAM of Dinophysis phenology"
  ) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic()

# Saving plot
# ggsave('Plots/GAMs/gam_Dino_4sites_crop.tiff', dpi = 300, height = 175, width = 250,
#        units = 'mm', compression = 'lzw')

# This is quite nice, but we might want to plot it as Dinophysis cells.L-1

pred_plot_cells_per_L <- response_pred_plot %>%
  mutate(across(-c('Code_point_Libelle', 'Day'),
                ~ .*100
  ))


# New plot in cells per L

ggplot(pred_plot_cells_per_L, aes(x = Day, y = median.fit, 
                                  color = Code_point_Libelle,
                                  fill = Code_point_Libelle)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
  geom_path(lwd = 1) +
  geom_point(data = Season_Dino_crop, aes(x = Day, y = Dinophysis_genus), 
             size = .8, alpha = .5) +
  # Add a "ghost point" to force the minimum y-axis range to 5
  geom_point(aes(x = 1, y = 500), color = 'transparent', fill = 'transparent',
             size = .8, alpha = .5) +
  # Labels
  labs(y = c(expression(paste("Dinophysis cells.L"^'-1'))),
       x = "Julian day",
       title = "Poisson GAM of Dinophysis phenology") +
  # facets
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic()

# Saving plot
# ggsave('Plots/GAMs/Dinophysis/gam_Dino_4sites_crop_cells_per_L.tiff', dpi = 300, height = 175, width = 250,
#        units = 'mm', compression = 'lzw')
