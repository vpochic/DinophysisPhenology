### GAMs of pigment concentrations ###
# Part of the Dinophysis Phenology project
# V. POCHIC 2025-06-05

### Required packages ####

library(tidyverse)
library(mgcv)
library(gratia)
library(ggpubr)
library(grid)

### Import data ####

Table_pigments <- read.csv2('Data/REPHY_outputs/Table_pigments_2007-2022.csv', header = TRUE,
                            fileEncoding = 'ISO-8859-1')

# Select points of interest for Dinophysis phenology
Table_pigments_select <- filter(Table_pigments, Code_point_Libelle %in%
                                  c('Antifer ponton pétrolier', 'Cabourg',
                                    'Men er Roue', 'Ouest Loscolo')) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# Separate the dataset
# Chlorophyll a
Table_CHLOROA <- Table_pigments_select %>%
  # Remove NAs
  filter(is.na(CHLOROA) == FALSE)

# Alloxanthin
Table_Allo <- Table_pigments_select %>%
  # Remove NAs
  filter(is.na(Allo) == FALSE)

# Alloxanthin/CHLOROA ratio
Table_Allo_ratio <- Table_pigments_select %>%
  # Remove NAs
  filter(is.na(Allo) == FALSE & is.na(CHLOROA) == FALSE & CHLOROA > 0) %>%
  mutate(Allo_ratio = Allo/CHLOROA)

# Fucoxanthin (and derivatives)
Table_Fuco <- Table_pigments_select %>%
  # Remove NAs
  filter(is.na(Fuco) == FALSE
         & is.na(But.fuco) == FALSE 
         & is.na(Hex.fuco) == FALSE) %>%
  mutate(Fuco_deriv = Fuco + But.fuco + Hex.fuco)

# Fucoxanthin/CHLOROA ratio
Table_Fuco_ratio <- Table_Fuco %>%
  # Remove NAs
  filter(is.na(CHLOROA) == FALSE & CHLOROA > 0) %>%
  mutate(Fuco_ratio = Fuco/CHLOROA) %>%
  mutate(Fuco_deriv_ratio = Fuco_deriv/CHLOROA)

# Peridinin
Table_Peri <- Table_pigments_select %>%
  # Remove NAs
  filter(is.na(Peri) == FALSE)

# Alloxanthin/CHLOROA ratio
Table_Peri_ratio <- Table_pigments_select %>%
  # Remove NAs
  filter(is.na(Peri) == FALSE & is.na(CHLOROA) == FALSE & CHLOROA > 0) %>%
  mutate(Peri_ratio = Peri/CHLOROA)

### GAM of chlorophyll a ####

Table_CHLOROA_factor <- Table_CHLOROA %>%
  mutate(Year = as_factor(Year))

gam_CHLOROA <- gam(data = Table_CHLOROA_factor, 
                   # Only a spline for the day of the year
                   # The use of a cyclic basis spline helps to make ends meet at the
                   # first and last days of the year
                   # 'k = -1' allows the model to fix the 'best' number of basic
                   # functions (= knots)
                   formula = CHLOROA~s(Day, bs = 'cc', k = -1,
                                       # separate each site
                                          by = Code_point_Libelle) 
                   # We add the Year as a random effect. This will help to assess and
                   # smooth the effects of interannual variability in the phenology
                   + s(Year, bs = 're', k = -1),
                   # Using a Gaussian distribution
                   family = scat(),
                   # Restricted maximum likelihood estimation (recommended method)
                   method = 'REML')

summary(gam_CHLOROA)
gam.check(gam_CHLOROA)

# Create a new 'empty' dataset for storing model prediction (fit)
gam_CHLOROA_newdata <- expand_grid(Day=seq(1, 365),
                                   # We add a Year vector as it has become a factor
                                   # of the model
                                   Year=seq(min(Table_CHLOROA$Year), 
                                            max(Table_CHLOROA$Year)),
                                   # We also add the site data
                                   Code_point_Libelle = unique(Table_CHLOROA$Code_point_Libelle))

## Get the inverse link function of the model
# With this function, we can transform the prediction we make on the link
#scale to the response scale
ilink <- gam_CHLOROA$family$linkinv

## Predict : add fit and se.fit on the **link** scale
gam_CHLOROA_newdata <- bind_cols(gam_CHLOROA_newdata, 
                                 setNames(as_tibble(
                                   predict(gam_CHLOROA, gam_CHLOROA_newdata, 
                                           se.fit = TRUE, type = 'link')[1:2]),
                                   c('fit_link','se_link')))

## Create the 95% confidence interval (2*standard error fit) and backtransform 
# to response variable using the inverse link function
gam_CHLOROA_newdata <- mutate(gam_CHLOROA_newdata,
                              fit_resp  = ilink(fit_link),
                              right_upr = ilink(fit_link + (2 * se_link)),
                              right_lwr = ilink(fit_link - (2 * se_link)))
# Check the confidence interval. It should not extend below 0 (negative chla
# values are impossible). But it does (slightly), and i don't know what to do about it.
min(gam_CHLOROA_newdata$right_lwr) # Negative :(
min(gam_CHLOROA_newdata$fit_resp) # Negative too, crap
max(gam_CHLOROA_newdata$right_upr)
max(gam_CHLOROA_newdata$fit_resp)

# First plot

ggplot(gam_CHLOROA_newdata, aes(x = Day, y = fit_resp))+
  geom_ribbon(aes(x=Day, ymin=right_lwr, ymax=right_upr), fill = 'grey70', alpha=0.7) +
  geom_line(linewidth = 1) +
  geom_point(data = Table_CHLOROA, aes(x = Day, y = CHLOROA, color = Year), shape = 21) +
  scale_color_viridis_c('Year') +
  facet_wrap(facets = c('Code_point_Libelle')) +
  theme_classic()

### Checking the model
ModelOutputs<-data.frame(Fitted=fitted(gam_CHLOROA),
                         Residuals=resid(gam_CHLOROA))

# Diagnostic plots

# A color palette for sites
pheno_palette4 <- c('red3', 'orangered', 
                    '#2156A1', '#5995E3')

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- gam_CHLOROA$model
# then we add the values of fitted and residuals 
# (but are they in the same order as the model? -> need to check that)
qq_data <- bind_cols(qq_data, ModelOutputs)

# Plot : verify that true data matches (more or less) model fit
ggplot(qq_data)+
  geom_point(aes(x = Day, y = CHLOROA), color = 'red') +
  geom_point(aes(x = Day, y = Fitted), color = 'blue') +
  facet_wrap(facets = c('Code_point_Libelle')) +
  theme_classic() +
  labs(y = "Chl a concentration (microgram/L)", x = "Calendar day")

# The highest values aren't matched well by the model. Maybe it will be better
# when the Year is taken into account as random effect

# And (qq-)plot
qqplot_custom <- ggplot(qq_data) +
  stat_qq(aes(sample=Residuals, color = Code_point_Libelle), alpha = .7) +
  stat_qq_line(aes(sample=Residuals, color = Code_point_Libelle)) +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  facet_wrap(facets = c('Code_point_Libelle')) +
  labs(title = 'qq-plot - Chlorophyll a GAM (2016-2022)', y="Sample Quantiles", 
       x="Theoretical Quantiles")

qqplot_custom

# Save the plot
# ggsave('Plots/GAMs/Pigments/qqplot_custom_chla_4sites.tiff', dpi = 300,
#        height = 175, width = 250, units = 'mm', compression = 'lzw')

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data)+
  geom_point(aes(x = Fitted, y = Residuals, 
                 color = Code_point_Libelle), alpha = .7) +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  facet_wrap(facets = c('Code_point_Libelle')) +
  labs(title = 'RvF-plot - Chlorophyll a GAM (2016-2022)', y="Residuals",
       x="Fitted Values")

RvFplot_custom
# Trumpet shaped plots here... Not great.

# Save the plot
# ggsave('Plots/GAMs/Pigments/RvFplot_custom_chla_4sites.tiff', dpi = 300,
#        height = 175, width = 250, units = 'mm', compression = 'lzw')

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data, aes(x = Residuals, fill = Code_point_Libelle))+
  geom_histogram(binwidth = 1)+
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  facet_wrap(facets = c('Code_point_Libelle')) +
  labs(title = 'Histogram of residuals - Chlorophyll a GAM 
(2016-2022)', x='Residuals', 
       y = 'Count')

HistRes_custom

# Save the plot
# ggsave('Plots/GAMs/Pigments/HistRes_custom_16sites.tiff', dpi = 300,
#        height = 175, width = 250, units = 'mm', compression = 'lzw')

# Plotting the whole model

ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_CHLOROA_newdata, aes(x = Day, ymin = right_lwr, 
                                              ymax = right_upr,
                                              color = Code_point_Libelle, 
                                              fill = Code_point_Libelle),
              linewidth = .75, alpha = .2) +
  geom_line(data = gam_CHLOROA_newdata, aes(x = Day, y = fit_resp,
                                            color = Code_point_Libelle),
            linewidth = 1) +
  # plot the data
  geom_point(data = Table_CHLOROA, aes(x = Day, y = CHLOROA,
                                       color = Code_point_Libelle, 
                                       fill = Code_point_Libelle), size = 2,
             alpha = .3) +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  facet_wrap(facets = c('Code_point_Libelle')) +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(title = 'Chlorophyll a concentration', x = 'Calendar day', 
       y = 'Pigment concentration (microgram/L)') +
  theme_classic()

# ggsave('GAM_chloroa.tiff', height = 164, width = 164,
#                dpi = 300, unit = 'mm', compression = 'lzw')


### Calculating confidence intervals ####

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
Vb <- vcov(gam_CHLOROA)
# New data
newd <- with(Table_CHLOROA_factor, data.frame(Day = rep(seq(1, 365, length = 365),
                                                     length(unique(Year))*
                                                       length(unique(Code_point_Libelle))
)))

newd <- group_by(newd, Day) %>%
  arrange(Day, by_group = TRUE)

# Add the Year vector
newd$Year <- rep(unique(Table_CHLOROA_factor$Year), 
                 365*length(unique(Table_CHLOROA_factor$Code_point_Libelle)))

newd <- group_by(newd, Year, Day) %>%
  arrange(Year, by_group = TRUE)

# And the site vector
newd$Code_point_Libelle <- rep(unique(Table_CHLOROA_factor$Code_point_Libelle), 
                               365*length(unique(Table_CHLOROA_factor$Year)))

newd <- ungroup(newd) %>%
  group_by(Code_point_Libelle) %>%
  arrange(Code_point_Libelle)

# Prediction by the model on the ***response*** scale, on the new data
pred <- predict(gam_CHLOROA, newd, se.fit = TRUE, type = 'response')
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
Cg <- predict(gam_CHLOROA, newd, type = "lpmatrix")
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


# Constructing a CI based on maximum and minimum simultaneous interval for
# plotting
# This is the response already because it was predicted on the response scale
response_pred_plot <- pred %>%
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

# Saving response_pred_plot so we don't have to re-run the model every time
# write.csv2(response_pred_plot, 'response_pred_plot_Allox_20241122_4sites.csv', row.names = FALSE,
#            fileEncoding = "ISO-8859-1")

# Transform values so they are expressed in the response variable units
# (by applying the inverse link function)
# response_pred_plot <- as.data.frame(pred_plot) %>%
#   mutate(across(c('median.fit', 'lwrS', 'uprS', 'lwrP', 'uprP'), ~ ilink(.)))

# And save that shall we
# write.csv2(response_pred_plot, 'response_pred_plot_Allox_20241122_4sites.csv',
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
# We can see that the simultaneous interval goes bonkers here

ggplot() +
  # plot the GAM
  geom_ribbon(data = response_pred_plot, aes(x = Day, ymin = lwrS, 
                                              ymax = uprS,
                                              color = Code_point_Libelle, 
                                              fill = Code_point_Libelle),
              linewidth = .75, alpha = .2) +
  geom_line(data = response_pred_plot, aes(x = Day, y = median.fit,
                                            color = Code_point_Libelle),
            linewidth = 1) +
  # plot the data
  geom_point(data = Table_CHLOROA, aes(x = Day, y = CHLOROA,
                                       color = Code_point_Libelle, 
                                       fill = Code_point_Libelle), size = 2,
             alpha = .3) +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  scale_y_continuous(limits = c(-1.5,15)) +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(title = 'GAM of Chlorophyll a concentration (2016-2022)', x = 'Calendar day', 
       y = 'Pigment concentration (microgram/L)') +
  theme_classic()

# Save the plot (even if it kinda looks like shit)
# ggsave('Plots/GAMs/Pigments/GAM_chla_4sites.tiff', height = 150, width = 150,
#        units = 'mm', dpi = 300, compression = 'lzw')

### GAM of Alloxanthin ####

# Transform Year as a factor for the random effect
Table_Allo_factor <- Table_Allo %>%
  mutate(Year = as_factor(Year)) %>%
  # And evacuate 4 dates for which the alloxanthin concentration is 0 (unlikely)
  filter(Allo > 0)

gam_Allo <- gam(data = Table_Allo_factor, 
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
gam.check(gam_Allo)

# gam.check indicates things look ok

# Create a new 'empty' dataset for storing model prediction (fit)
gam_Allo_newdata <- expand_grid(Day=seq(1, 365),
                                   # We add a Year vector as it has become a factor
                                   # of the model
                                   Year=seq(min(Table_Allo$Year), 
                                            max(Table_Allo$Year)),
                                   # We also add the site data
                                   Code_point_Libelle = unique(Table_Allo$Code_point_Libelle))

## Get the inverse link function of the model
# With this function, we can transform the prediction we make on the link
#scale to the response scale
ilink <- gam_Allo$family$linkinv

## Predict : add fit and se.fit on the **link** scale
gam_Allo_newdata <- bind_cols(gam_Allo_newdata, 
                                 setNames(as_tibble(
                                   predict(gam_Allo, gam_Allo_newdata, 
                                           se.fit = TRUE, type = 'link')[1:2]),
                                   c('fit_link','se_link')))

## Create the 95% confidence interval (2*standard error fit) and backtransform 
# to response variable using the inverse link function
gam_Allo_newdata <- mutate(gam_Allo_newdata,
                              fit_resp  = ilink(fit_link),
                              right_upr = ilink(fit_link + (2 * se_link)),
                              right_lwr = ilink(fit_link - (2 * se_link)))
# Check the confidence interval. It should not extend below 0 (negative pigment
# concentrations are impossible)
min(gam_Allo_newdata$right_lwr) # Nice :)
min(gam_Allo_newdata$fit_resp)
max(gam_Allo_newdata$right_upr) # Nice too
max(gam_Allo_newdata$fit_resp)
# But we can see that the confidence interval is ridiculously small around the
# min and max of the model fit.

# Plot

ggplot(gam_Allo_newdata, aes(x = Day, y = fit_resp))+
  geom_ribbon(aes(x=Day, ymin=right_lwr, ymax=right_upr), fill = 'grey70', alpha=0.7) +
  geom_line(linewidth = 1) +
  geom_point(data = Table_Allo, aes(x = Day, y = Allo, color = Year), shape = 21) +
  scale_color_viridis_c('Year') +
  facet_wrap(facets = c('Code_point_Libelle')) +
  theme_classic()

# Saving plot
# ggsave('gam_Allo_all.tiff', dpi = 300, height = 120, width = 160, 
#        units = 'mm', compression = 'lzw')

### Checking the model
ModelOutputs<-data.frame(Fitted=fitted(gam_Allo),
                         Residuals=resid(gam_Allo))

# Diagnostic plots

# A color palette for sites
pheno_palette4 <- c('red3', 'orangered', 
                    '#2156A1', '#5995E3')

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- gam_Allo$model
# then we add the values of fitted and residuals 
# (but are they in the same order as the model? -> need to check that)
qq_data <- bind_cols(qq_data, ModelOutputs)

# Plot : verify that true data matches (more or less) model fit
ggplot(qq_data)+
  geom_point(aes(x = Day, y = Allo), color = 'red') +
  geom_point(aes(x = Day, y = Fitted), color = 'blue') +
  facet_wrap(facets = c('Code_point_Libelle')) +
  theme_classic() +
  labs(y = "Chl a concentration (microgram/L)", x = "Calendar day")

# It matches! great!

# And (qq-)plot
qqplot_custom <- ggplot(qq_data) +
  stat_qq(aes(sample=Residuals, color = Code_point_Libelle), alpha = .7) +
  stat_qq_line(aes(sample=Residuals, color = Code_point_Libelle)) +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  facet_wrap(facets = c('Code_point_Libelle')) +
  labs(title = 'qq-plot for alloxanthin GAM', y="Sample Quantiles",
       x="Theoretical Quantiles")

qqplot_custom

# Save the plot
# ggsave('qqplot_allox_4sites.tiff', dpi = 300, height = 175, width = 250,
#                units = 'mm', compression = 'lzw')

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data)+
  geom_point(aes(x = Fitted, y = Residuals, 
                 color = Code_point_Libelle), alpha = .7) +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free') +
  labs(title = 'Residuals vs fitted (alloxanthin GAM)',
       y="Residuals", x="Fitted Values")

RvFplot_custom
# Not much structure : quite good

# Save the plot
# ggsave('RvFplot_allox_4sites.tiff', dpi = 300, height = 175, width = 250,
#                units = 'mm', compression = 'lzw')

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data, aes(x = Residuals, fill = Code_point_Libelle))+
  geom_histogram(binwidth = 1)+
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  facet_wrap(facets = c('Code_point_Libelle')) +
  labs(title = 'Histogram of residuals (alloxanthin GAM)',
    x='Residuals', y = 'Count')

HistRes_custom
# Mostly centered on 0 -> nice

# Save the plot
# ggsave('HistRes_allox_4sites.tiff', dpi = 300, height = 175, width = 250,
#                units = 'mm', compression = 'lzw')

### Calculating confidence intervals ####

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
newd <- with(Table_Allo_factor, data.frame(Day = rep(seq(1, 365, length = 365),
                                               length(unique(Year))*
                                                 length(unique(Code_point_Libelle))
)))

newd <- group_by(newd, Day) %>%
  arrange(Day, by_group = TRUE)

# Add the Year vector
newd$Year <- rep(unique(Table_Allo_factor$Year), 
                 365*length(unique(Table_Allo_factor$Code_point_Libelle)))

newd <- group_by(newd, Year, Day) %>%
  arrange(Year, by_group = TRUE)

# And the site vector
newd$Code_point_Libelle <- rep(unique(Table_Allo_factor$Code_point_Libelle), 
                               365*length(unique(Table_Allo_factor$Year)))

newd <- ungroup(newd) %>%
  group_by(Code_point_Libelle) %>%
  arrange(Code_point_Libelle)

# Prediction by the model on the ***response*** scale, on the new data
pred <- predict(gam_Allo, newd, se.fit = TRUE, type = 'response')
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


# Constructing a CI based on maximum and minimum simultaneous interval for
# plotting
# This is the response already because it was predicted on the response scale
response_pred_plot <- pred %>%
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

# Saving response_pred_plot so we don't have to re-run the model every time
# write.csv2(response_pred_plot, 'response_pred_plot_Allox_20241122_4sites.csv', row.names = FALSE,
#            fileEncoding = "ISO-8859-1")

# Transform values so they are expressed in the response variable units
# (by applying the inverse link function)
# response_pred_plot <- as.data.frame(pred_plot) %>%
#   mutate(across(c('median.fit', 'lwrS', 'uprS', 'lwrP', 'uprP'), ~ ilink(.)))

# And save that shall we
# write.csv2(response_pred_plot, 'response_pred_plot_Allox_20241122_4sites.csv',
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
# We can see that the simultaneous interval goes bonkers here

### Plotting the whole model ####

ggplot() +
  # plot the GAM
  geom_ribbon(data = response_pred_plot, aes(x = Day, ymin = lwrP, 
                                              ymax = uprP,
                                              color = Code_point_Libelle, 
                                              fill = Code_point_Libelle),
              linewidth = .75, alpha = .2) +
  geom_line(data = response_pred_plot, aes(x = Day, y = median.fit,
                                            color = Code_point_Libelle),
            linewidth = 1) +
  # plot the data
  geom_point(data = Table_Allo_factor, aes(x = Day, y = Allo,
                                       color = Code_point_Libelle, 
                                       fill = Code_point_Libelle), size = 2,
             alpha = .3) +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(title = 'GAM of alloxanthin concentration', x = 'Calendar day', 
       y = 'Pigment concentration (microgram/L)') +
  theme_classic()

# ggsave('GAM_Allo.tiff', height = 164, width = 164,
#                dpi = 300, unit = 'mm', compression = 'lzw')

### GAM of Alloxanthin/CHLOROA ratio ####

gam_Allo_ratio <- gam(data = Table_Allo_ratio, 
                   # Only a spline for the day of the year
                   # The use of a cyclic basis spline helps to make ends meet at the
                   # first and last days of the year
                   # 'k = -1' allows the model to fix the 'best' number of basic
                   # functions (= knots)
                   formula = Allo_ratio~s(Day, bs = 'cc', k = -1,
                                       # separate each site
                                       by = Code_point_Libelle) 
                   # We add the Year as a random effect. This will help to assess and
                   # smooth the effects of interannual variability in the phenology
                   + s(Year, bs = 're', k = -1),
                   # Using a Gaussian distribution
                   family = gaussian(),
                   # Restricted maximum likelihood estimation (recommended method)
                   method = 'REML')

summary(gam_Allo_ratio)
gam.check(gam_Allo_ratio)

# gam.check indicates that there might be an issue with the random effect
# smooth for 'Year' (low p-value)

# Create a new 'empty' dataset for storing model prediction (fit)
gam_Allo_ratio_newdata <- expand_grid(Day=seq(1, 365),
                                   # We add a Year vector as it has become a factor
                                   # of the model
                                   Year=seq(min(Table_Allo_ratio$Year), 
                                            max(Table_Allo_ratio$Year)),
                                   # We also add the site data
                                   Code_point_Libelle = unique(Table_Allo_ratio$Code_point_Libelle))

## Get the inverse link function of the model
# With this function, we can transform the prediction we make on the link
#scale to the response scale
ilink <- gam_Allo_ratio$family$linkinv

## Predict : add fit and se.fit on the **link** scale
gam_Allo_ratio_newdata <- bind_cols(gam_Allo_ratio_newdata, 
                                 setNames(as_tibble(
                                   predict(gam_Allo_ratio, gam_Allo_ratio_newdata, 
                                           se.fit = TRUE, type = 'link')[1:2]),
                                   c('fit_link','se_link')))

## Create the 95% confidence interval (2*standard error fit) and backtransform 
# to response variable using the inverse link function
gam_Allo_ratio_newdata <- mutate(gam_Allo_ratio_newdata,
                              fit_resp  = ilink(fit_link),
                              right_upr = ilink(fit_link + (2 * se_link)),
                              right_lwr = ilink(fit_link - (2 * se_link)))
# Check the confidence interval. It should not extend below 0 (negative counts
# are impossible)
min(gam_Allo_ratio_newdata$right_lwr) # Nice :)
min(gam_Allo_ratio_newdata$fit_resp)
max(gam_Allo_ratio_newdata$right_upr) # Nice too
max(gam_Allo_ratio_newdata$fit_resp)
# But we can see that the confidence interval is ridiculously small around the
# min and max of the model fit.

# Plot

ggplot(gam_Allo_ratio_newdata, aes(x = Day, y = fit_resp))+
  geom_ribbon(aes(x=Day, ymin=right_lwr, ymax=right_upr), fill = 'grey70', alpha=0.7) +
  geom_line(linewidth = 1) +
  geom_point(data = Table_Allo_ratio, aes(x = Day, y = Allo_ratio, color = Year), shape = 21) +
  scale_color_viridis_c('Year') +
  facet_wrap(facets = c('Code_point_Libelle')) +
  theme_classic()

# Saving plot
# ggsave('gam_Allo_ratio_all.tiff', dpi = 300, height = 120, width = 160, 
#        units = 'mm', compression = 'lzw')

### Checking the model
ModelOutputs<-data.frame(Fitted=fitted(gam_Allo_ratio),
                         Residuals=resid(gam_Allo_ratio))

# Diagnostic plots

# A color palette for sites
pheno_palette4 <- c('red3', 'orangered', 
                    '#2156A1', '#5995E3')

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- gam_Allo_ratio$model
# then we add the values of fitted and residuals 
# (but are they in the same order as the model? -> need to check that)
qq_data <- bind_cols(qq_data, ModelOutputs)

# Plot : verify that true data matches (more or less) model fit
ggplot(qq_data)+
  geom_point(aes(x = Day, y = Allo_ratio), color = 'red') +
  geom_point(aes(x = Day, y = Fitted), color = 'blue') +
  facet_wrap(facets = c('Code_point_Libelle')) +
  theme_classic() +
  labs(y = "Chl a concentration (microgram/L)", x = "Calendar day")

# It matches! great!

# And (qq-)plot
qqplot_custom <- ggplot(qq_data) +
  stat_qq(aes(sample=Residuals, color = Code_point_Libelle), alpha = .7) +
  stat_qq_line(aes(sample=Residuals, color = Code_point_Libelle)) +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  facet_wrap(facets = c('Code_point_Libelle')) +
  labs(y="Sample Quantiles",x="Theoretical Quantiles")

qqplot_custom

# Save the plot
# ggsave('qqplot_custom_16sites.tiff', dpi = 300, height = 175, width = 250,
#                units = 'mm', compression = 'lzw')

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data)+
  geom_point(aes(x = Fitted, y = Residuals, 
                 color = Code_point_Libelle), alpha = .7) +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  facet_wrap(facets = c('Code_point_Libelle')) +
  labs(y="Residuals",x="Fitted Values")

RvFplot_custom
# Trumpet shaped plots here...

# Save the plot
# ggsave('RvFplot_custom_16sites.tiff', dpi = 300, height = 175, width = 250,
#                units = 'mm', compression = 'lzw')

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data, aes(x = Residuals, fill = Code_point_Libelle))+
  geom_histogram(binwidth = 1)+
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  facet_wrap(facets = c('Code_point_Libelle')) +
  labs(x='Residuals', y = 'Count')

HistRes_custom

# Save the plot
# ggsave('HistRes_custom_16sites.tiff', dpi = 300, height = 175, width = 250,
#                units = 'mm', compression = 'lzw')

# Plotting the whole model

ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Allo_ratio_newdata, aes(x = Day, ymin = right_lwr, 
                                              ymax = right_upr,
                                              color = Code_point_Libelle, 
                                              fill = Code_point_Libelle),
              linewidth = .75, alpha = .2) +
  geom_line(data = gam_Allo_ratio_newdata, aes(x = Day, y = fit_resp,
                                            color = Code_point_Libelle),
            linewidth = 1) +
  # plot the data
  geom_point(data = Table_Allo_ratio, aes(x = Day, y = Allo_ratio,
                                       color = Code_point_Libelle, 
                                       fill = Code_point_Libelle), size = 2,
             alpha = .3) +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(title = 'Alloxanthin/chlorophyll a ratio', x = 'Calendar day', 
       y = 'Alloxanthin/chlorophyll a (-)') +
  theme_classic()

# ggsave('GAM_Allo_ratio.tiff', height = 164, width = 164,
#                dpi = 300, unit = 'mm', compression = 'lzw')

#### Comparison with GAMs of Dinophysis and Mesodinium ####

## Import data ####
# GAM data ####
## Dinophysis GAM data
gam_Dino <- read.csv2('response_pred_plot_20241120_16sites.csv',
                                header = TRUE, fileEncoding = 'ISO-8859-1')
# Separate sites
gam_Dino_Antifer <- filter(gam_Dino, 
                           Code_point_Libelle == 'Antifer ponton pétrolier')
gam_Dino_Cabourg <- filter(gam_Dino, Code_point_Libelle == 'Cabourg')
gam_Dino_MeR <- filter(gam_Dino, Code_point_Libelle == 'Men er Roue')
gam_Dino_OL <- filter(gam_Dino, Code_point_Libelle == 'Ouest Loscolo')

## Mesodinium GAM data
gam_Meso <- read.csv2('response_pred_plot_MESO_20241121.csv', header = TRUE,
                       fileEncoding = 'ISO-8859-1')

# Separate sites
gam_Meso_Antifer <- filter(gam_Meso, 
                           Code_point_Libelle == 'Antifer ponton pétrolier')
gam_Meso_Cabourg <- filter(gam_Meso, Code_point_Libelle == 'Cabourg')
gam_Meso_MeR <- filter(gam_Meso, Code_point_Libelle == 'Men er Roue')
gam_Meso_OL <- filter(gam_Meso, Code_point_Libelle == 'Ouest Loscolo')

# Alloxanthin data as well!
# Import data
gam_Allo_plot <- read.csv2('response_pred_plot_Allox_20241122_4sites.csv',
                           header = TRUE, fileEncoding = 'ISO-8859-1')

gam_Allo_Antifer <- filter(gam_Allo_plot, 
                           Code_point_Libelle == 'Antifer ponton pétrolier')
gam_Allo_Cabourg <- filter(gam_Allo_plot, Code_point_Libelle == 'Cabourg')
gam_Allo_MeR <- filter(gam_Allo_plot, Code_point_Libelle == 'Men er Roue')
gam_Allo_OL <- filter(gam_Allo_plot, Code_point_Libelle == 'Ouest Loscolo')

# Count data ####
# The sites must be entered as factors
Season_Dino <- read.csv2('Season_Dino.csv', header = TRUE, 
                         fileEncoding = 'ISO-8859-1') %>%
  mutate(Code_point_Libelle = as.factor(Code_point_Libelle))

# Separate sites
Season_Dino_Antifer <- filter(Season_Dino, 
                           Code_point_Libelle == 'Antifer ponton pétrolier')
Season_Dino_Cabourg <- filter(Season_Dino, Code_point_Libelle == 'Cabourg')
Season_Dino_MeR <- filter(Season_Dino, Code_point_Libelle == 'Men er Roue')
Season_Dino_OL <- filter(Season_Dino, Code_point_Libelle == 'Ouest Loscolo')

# Alright!

### Mesodinium
# The table of Mesodinium counts
Season_Meso <- read.csv2('Season_Meso.csv', header = TRUE, 
                              fileEncoding = 'ISO-8859-1')

# Separate sites
Season_Meso_Antifer <- filter(Season_Meso, 
                              Code_point_Libelle == 'Antifer ponton pétrolier')
Season_Meso_Cabourg <- filter(Season_Meso, Code_point_Libelle == 'Cabourg')
Season_Meso_MeR <- filter(Season_Meso, Code_point_Libelle == 'Men er Roue')
Season_Meso_OL <- filter(Season_Meso, Code_point_Libelle == 'Ouest Loscolo')

# For Alloxanthin data also (not count data but eh)
Season_Allo_Antifer <- filter(Table_Allo_factor, 
                              Code_point_Libelle == 'Antifer ponton pétrolier')
Season_Allo_Cabourg <- filter(Table_Allo_factor, Code_point_Libelle == 'Cabourg')
Season_Allo_MeR <- filter(Table_Allo_factor, Code_point_Libelle == 'Men er Roue')
Season_Allo_OL <- filter(Table_Allo_factor, Code_point_Libelle == 'Ouest Loscolo')


### Plots ####

pheno_palette4 <- c('red3', 'orangered', 
                    '#2156A1', '#5995E3')
## Antifer ####

# Alloxanthin
plot_Allo_Antifer <- ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Allo_Antifer, aes(x = Day, ymin = lwrP, 
                                              ymax = uprP),
              linewidth = .75, alpha = .2,
              color = 'red3', fill = 'red3') +
  geom_line(data = gam_Allo_Antifer, aes(x = Day, y = median.fit),
            linewidth = 1, color = 'red3') +
  # plot the data
  geom_point(data = Season_Allo_Antifer, aes(x = Day, y = Allo), size = 2,
             alpha = .3, color = 'red3') +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(title = 'Antifer ponton pétrolier',
       subtitle = 'Alloxanthin', x = 'Calendar day', 
       y = 'Pigment concentration (microgram/L)') +
  theme_classic()

plot_Allo_Antifer

# Mesodinium
plot_Meso_Antifer <- ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Meso_Antifer, aes(x = Day, ymin = lwrP, 
                                           ymax = uprP),
              linewidth = .75, alpha = .2,
              color = 'red3', fill = 'red3') +
  geom_line(data = gam_Meso_Antifer, aes(x = Day, y = median.fit),
            linewidth = 1, color = 'red3') +
  # plot the data
  geom_point(data = Season_Meso_Antifer, aes(x = Day, y = true_count), size = 2,
             alpha = .3, color = 'red3') +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(subtitle = 'Mesodinium', x = 'Calendar day', 
       y = 'Mesodinium cells observed in 10 mL') +
  theme_classic()

plot_Meso_Antifer

# Dinophysis
plot_Dino_Antifer <- ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Dino_Antifer, aes(x = Day, ymin = lwrP, 
                                           ymax = uprP),
              linewidth = .75, alpha = .2,
              color = 'red3', fill = 'red3') +
  geom_line(data = gam_Dino_Antifer, aes(x = Day, y = median.fit),
            linewidth = 1, color = 'red3') +
  # plot the data
  geom_point(data = Season_Dino_Antifer, aes(x = Day, y = true_count), size = 2,
             alpha = .3, color = 'red3') +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(subtitle = 'Dinophysis', x = 'Calendar day', 
       y = 'Dinophysis cells observed in 10 mL') +
  theme_classic()

plot_Dino_Antifer

# Arrange the plot
ggarrange(plot_Allo_Antifer, plot_Meso_Antifer, plot_Dino_Antifer, nrow = 3,
          align = 'v')

# Save the plot
# ggsave('Succession_plot_Antifer.tiff', height = 270, width = 140,
#                dpi = 300, unit = 'mm', compression = 'lzw')

## Cabourg ####

# Alloxanthin
plot_Allo_Cabourg <- ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Allo_Cabourg, aes(x = Day, ymin = lwrP, 
                                           ymax = uprP),
              linewidth = .9, alpha = .2,
              color = 'orangered', fill = 'orangered') +
  geom_line(data = gam_Allo_Cabourg, aes(x = Day, y = median.fit),
            linewidth = 1.25, color = 'orangered') +
  # plot the data
  geom_point(data = Season_Allo_Cabourg, aes(x = Day, y = Allo), size = 3,
             alpha = .3, color = 'orangered') +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(title = 'Cabourg',
       subtitle = 'Alloxanthin', x = 'Calendar day', 
       y = 'Alloxanthin') +
  theme_classic()

plot_Allo_Cabourg

# Mesodinium
# We crop out some points for visualisation
Season_Meso_Cabourg_crop <- Season_Meso_Cabourg %>%
  filter(true_count < 100)

plot_Meso_Cabourg <- ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Meso_Cabourg, aes(x = Day, ymin = lwrP, 
                                           ymax = uprP),
              linewidth = .9, alpha = .2,
              color = 'orangered', fill = 'orangered') +
  geom_line(data = gam_Meso_Cabourg, aes(x = Day, y = median.fit),
            linewidth = 1.25, color = 'orangered') +
  # plot the data
  geom_point(data = Season_Meso_Cabourg_crop, aes(x = Day, y = true_count), size = 3,
             alpha = .3, color = 'orangered') +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(subtitle = 'Mesodinium', x = 'Calendar day', 
       y = 'Mesodinium cells') +
  theme_classic()

plot_Meso_Cabourg

# Dinophysis
# We crop out some points for visualisation
Season_Dino_Cabourg_crop <- Season_Dino_Cabourg %>%
  filter(true_count < 50)

# Import the derivative plot
dino_deriv <- read.csv2('gam_Dino_multiyear_deriv_plot_select_20241204.csv',
           header = TRUE, fileEncoding = 'ISO-8859-1')
dino_deriv_Cabourg <- filter(dino_deriv, Code_point_Libelle == 'Cabourg')

plot_Dino_Cabourg <- ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Dino_Cabourg, aes(x = Day, ymin = lwrP, 
                                           ymax = uprP),
              linewidth = .75, alpha = .2,
              color = 'orangered', fill = 'orangered') +
  geom_line(data = gam_Dino_Cabourg, aes(x = Day, y = median.fit),
            linewidth = 1, color = 'orangered') +
  # plot the data
  geom_point(data = Season_Dino_Cabourg_crop, aes(x = Day, y = true_count), size = 3,
             alpha = .3, color = 'orangered') +
  # Chenge the color palette
  # new_scale_color() +
  # scale_color_cmocean(name = 'curl', direction = -1, guide = "none") +
  # # Rugplot with the derivative of the Dinophysis GAM
  # geom_rug(data = dino_deriv_Cabourg, aes(x = Day, color = median.deriv),
  #          linewidth = 1,
  #          length = unit(0.2, 'cm'),
  #) +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(subtitle = 'Dinophysis', x = 'Calendar day', 
       y = 'Dinophysis cells') +
  theme_classic()

plot_Dino_Cabourg

# Arrange the plot
succession_plot_cabourg <- ggarrange(plot_Allo_Cabourg + rremove("xlab"),
                                     plot_Meso_Cabourg + rremove("xlab"),
                                     plot_Dino_Cabourg + rremove("xlab"),
                                     nrow = 3,
                                     align = 'v')

# Common y-axis
annotate_figure(succession_plot_cabourg,
                bottom = textGrob('Julian day'))

# Save the plot
# ggsave('Succession_plot_Cabourg2.tiff', dpi = 300, height = 164, width = 164,
#        units = 'mm', compression = 'lzw')

## Men er Roue ####

# Alloxanthin
plot_Allo_MeR <- ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Allo_MeR, aes(x = Day, ymin = lwrP, 
                                           ymax = uprP),
              linewidth = .75, alpha = .2,
              color = '#2156A1', fill = '#2156A1') +
  geom_line(data = gam_Allo_MeR, aes(x = Day, y = median.fit),
            linewidth = 1, color = '#2156A1') +
  # plot the data
  geom_point(data = Season_Allo_MeR, aes(x = Day, y = Allo), size = 2,
             alpha = .3, color = '#2156A1') +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(title = 'Men er Roue',
       subtitle = 'Alloxanthin', x = 'Calendar day', 
       y = 'Pigment concentration (microgram/L)') +
  theme_classic()

plot_Allo_MeR

# Mesodinium
plot_Meso_MeR <- ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Meso_MeR, aes(x = Day, ymin = lwrP, 
                                           ymax = uprP),
              linewidth = .75, alpha = .2,
              color = '#2156A1', fill = '#2156A1') +
  geom_line(data = gam_Meso_MeR, aes(x = Day, y = median.fit),
            linewidth = 1, color = '#2156A1') +
  # plot the data
  geom_point(data = Season_Meso_MeR, aes(x = Day, y = true_count), size = 2,
             alpha = .3, color = '#2156A1') +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(subtitle = 'Mesodinium', x = 'Calendar day', 
       y = 'Mesodinium cells observed in 10 mL') +
  theme_classic()

plot_Meso_MeR

# Dinophysis
plot_Dino_MeR <- ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Dino_MeR, aes(x = Day, ymin = lwrP, 
                                           ymax = uprP),
              linewidth = .75, alpha = .2,
              color = '#2156A1', fill = '#2156A1') +
  geom_line(data = gam_Dino_MeR, aes(x = Day, y = median.fit),
            linewidth = 1, color = '#2156A1') +
  # plot the data
  geom_point(data = Season_Dino_MeR, aes(x = Day, y = true_count), size = 2,
             alpha = .3, color = '#2156A1') +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(subtitle = 'Dinophysis', x = 'Calendar day', 
       y = 'Dinophysis cells observed in 10 mL') +
  theme_classic()

plot_Dino_MeR

# Arrange the plot
ggarrange(plot_Allo_MeR, plot_Meso_MeR, plot_Dino_MeR, nrow = 3,
          align = 'v')

# Save the plot
# ggsave('Succession_plot_MeR.tiff', height = 270, width = 140,
#                dpi = 300, unit = 'mm', compression = 'lzw')
## Ouest Loscolo ####

# Alloxanthin
plot_Allo_OL <- ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Allo_OL, aes(x = Day, ymin = lwrP, 
                                       ymax = uprP),
              linewidth = .9, alpha = .2,
              color = '#5995E3', fill = '#5995E3') +
  geom_line(data = gam_Allo_OL, aes(x = Day, y = median.fit),
            linewidth = 1.25, color = '#5995E3') +
  # plot the data
  geom_point(data = Season_Allo_OL, aes(x = Day, y = Allo), size = 3,
             alpha = .3, color = '#5995E3') +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(title = 'Ouest Loscolo',
       subtitle = 'Alloxanthin', x = NULL, 
       y = 'Alloxanthin') +
  theme_classic()

plot_Allo_OL

# Mesodinium
plot_Meso_OL <- ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Meso_OL, aes(x = Day, ymin = lwrP, 
                                       ymax = uprP),
              linewidth = .9, alpha = .2,
              color = '#5995E3', fill = '#5995E3') +
  geom_line(data = gam_Meso_OL, aes(x = Day, y = median.fit),
            linewidth = 1.25, color = '#5995E3') +
  # plot the data
  geom_point(data = Season_Meso_OL, aes(x = Day, y = true_count), size = 3,
             alpha = .3, color = '#5995E3') +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(subtitle = 'Mesodinium', x = NULL, 
       y = 'Mesodinium cells') +
  theme_classic()

plot_Meso_OL

# Dinophysis
dino_deriv_OL <- filter(dino_deriv, Code_point_Libelle == 'Ouest Loscolo')

plot_Dino_OL <- ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Dino_OL, aes(x = Day, ymin = lwrP, 
                                       ymax = uprP),
              linewidth = .9, alpha = .2,
              color = '#5995E3', fill = '#5995E3') +
  geom_line(data = gam_Dino_OL, aes(x = Day, y = median.fit),
            linewidth = 1.25, color = '#5995E3') +
  # plot the data
  geom_point(data = Season_Dino_OL, aes(x = Day, y = true_count), size = 3,
             alpha = .3, color = '#5995E3') +
  # # Chenge the color palette
  # new_scale_color() +
  # scale_color_cmocean(name = 'curl', direction = -1, guide = "none") +
  # # Rugplot with the derivative of the Dinophysis GAM
  # geom_rug(data = dino_deriv_OL, aes(x = Day, color = median.deriv),
  #          linewidth = 1,
  #          length = unit(0.2, 'cm'),
  # ) +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(subtitle = 'Dinophysis', x = 'Calendar day', 
       y = 'Dinophysis cells') +
  theme_classic()

plot_Dino_OL

# Arrange the plot
succession_plot_OL <- ggarrange(plot_Allo_OL + rremove("xlab"),
                                     plot_Meso_OL + rremove("xlab"),
                                     plot_Dino_OL + rremove("xlab"),
                                     nrow = 3,
                                     align = 'v')

# Common y-axis
annotate_figure(succession_plot_OL,
                bottom = textGrob('Julian day'))

# Save the plot
# ggsave('Succession_plot_OL2.tiff', dpi = 300, height = 164, width = 164,
#        units = 'mm', compression = 'lzw')