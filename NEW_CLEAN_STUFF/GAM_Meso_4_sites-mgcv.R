###### GAM for Mesodinium phenology with mgcv ###
## V. POCHIC
# 2025-06-04

# This is more or less the same script than for the Dinophysis GAM, but for 
# its ciliate prey Mesodinium. We restrict the analysis to 4 sites (in Southern
# Brittany and Seine Bay), and to years from 2016 to 2022.

#### Packages and functions ####
library(tidyverse)
library(ggplot2)
library(mgcv)
library(gratia)
library(viridis)
library(cmocean)

#### Import data ####

# Table with phytoplankton counts
Table_phyto_taxon <- read.csv2('Data/REPHY_outputs/Table1_phyto_taxon.csv', fileEncoding = "ISO-8859-1")

# Select only 4 sites of interest for Mesodinium from 2016 --> 2022
Table_phyto_sites <- filter(Table_phyto_taxon, Code_point_Libelle %in% 
                           c(# Baie de Seine
                             'Antifer ponton pétrolier', 'Cabourg',
                             # Bretagne Sud
                             'Men er Roue', 'Ouest Loscolo')
                         ) %>%
  filter(Year >= 2016 & Year <= 2022)

########## Organizing the data #######

#### More zeros ####

### The goal here is to create a table where all absences of detection are 
# listed as zeros
Table_Meso_zeros <- Table_phyto_sites %>%
  pivot_wider(names_from = Taxon, values_from = Comptage, values_fill = 0) %>%
  # Create Latitude and Longitude variables for selection (because lon and lat
  # introduce a bunch of unwanted taxa )
  mutate(Longitude = as.numeric(lon)) %>%
  mutate(Latitude = as.numeric(lat)) %>%
  select(starts_with('Mesodinium') |
           contains(c('Code.Region', 'Code_point_Libelle', 'Year', 'Month',
                      'Date', 'Heure', 'Code.parametre', 'SALI', 'TEMP', 
                      'Longitude', 'Latitude', 'ID.interne.passage')))

Table_Meso_zeros <- Table_Meso_zeros %>%
  # Getting rid of some unwelcome guests (who got onboard because the taxon
  # name contains 'sali')
  select(-(contains('Diplo')))

# Save Table_Meso_zeros so we don't have to re-run the model to plot again
# write.csv2(Table_Meso_zeros, 'Data/REPHY_outputs/Table_Meso_zeros.csv',
#            row.names = FALSE, fileEncoding = "ISO-8859-1")


### Clean the environment
rm(Table_phyto_taxon)
rm(Table_phyto_sites)

#### Seasonality ####

### Create a seasonality dataset (with the calendar day/julian day variable)
Season_Meso <- Table_Meso_zeros %>%
  # Only FLORTOT
  filter(Code.parametre == 'FLORTOT') %>%
  # create a calendar day variable
  mutate(Day = as.numeric(yday(Date))) %>%
  mutate(Date = ymd(Date)) %>%
  # Transform counts that do not correspond to the REPHY protocol (not dividable
  # by 100) into 0s
  mutate(across(.cols = c('Mesodinium', 'Mesodinium rubrum'), 
                .fns = ~ ifelse(. %% 100 != 0, 0, .))) %>%
  # And creating some count variables for Mesodinium as a genus
  mutate(Mesodinium_genus = rowSums(across(contains('Mesodinium')))) %>%
  # create the log of abundance + 1
  mutate(log_c = log10(Mesodinium_genus+1)) %>%
  # create a 'true count' variable
  # this variable corresponds to the number of cells that were actually
  # counted by the operator (in 10 mL). This variable will follow a Poisson
  # distribution, contrary to Mesodinium_genus because the conversion from
  # 10 mL to 1L (*100) prevents some intermediate values (e.g., 150 cells.L-1)
  mutate(true_count = Mesodinium_genus/100) %>%
  # Filter out exceptionnaly high counts (>500 cells observed in 10 mL)
  # this represents 1 event
  filter(true_count < 500) %>%
  # converting the site to factor for the model
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo')) %>%
  # Expanding on date variable (creating a fortnight variable)
  mutate(Week = week(Date)) %>%
  # 'ceiling' takes the upper integer of a decimal number
  mutate(Fortnight = ceiling(Week/2))

# Checking the homogeneity of the sampling between years
# We make a histogram of samples for the 16 years (2007-2022)

# Color palette
pheno_palette4 <- c('red3', 'orangered', '#2156A1', '#5995E3')

ggplot(data = Season_Meso) +
  geom_histogram(aes(x = as_factor(Year), fill = Code_point_Libelle),
                 stat = 'count') +
  facet_wrap(facets = 'Code_point_Libelle') +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  scale_x_discrete(breaks = c('2016', '2020', '2022')) +
  theme_classic()
# Mainly homogenous. Atlantic sites 
# experienced a drop in 2020 due to covid19 pandemic.

# Homogeneity across the year, for the 26 fortnights in a year 
# (fortnightly sampling)
ggplot(data = Season_Meso) +
  geom_histogram(aes(x = as_factor(Fortnight), fill = Code_point_Libelle),
                 stat = 'count') +
  facet_wrap(facets = 'Code_point_Libelle') +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  scale_x_discrete(breaks = c('1', '5', '10', '15', '20', '26')) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  theme_classic()
# Relatively OK

# This seems quite fine!!!

# Let's save that
# write.csv2(Season_Meso, 'Data/REPHY_outputs/Season_Meso_20250604.csv', row.names = FALSE,
#            fileEncoding = "ISO-8859-1")


### Plotting the dataset for one site as an example

plotMeso_OL_nozeros <- filter(Season_Meso, (Code_point_Libelle == 'Ouest Loscolo'
                                            & true_count != 0))
plotMeso_OL <- filter(Season_Meso, (Code_point_Libelle == 'Ouest Loscolo'))

# Plot the Dinophysis observations as function of date (timeseries view)
ggplot(plotMeso_OL_nozeros)+
  geom_point(aes(x=Date, y=true_count), color = '#2156A1',
             size = 3, alpha = .8) +
  theme_classic() +
  labs(title = 'Ouest Loscolo',
         y="Mesodinium cells observed in 10 mL", x="Date")


### Plot the observations per year

ggplot(plotMeso_OL_nozeros)+
  geom_point(aes(x=Day, y=true_count), color = '#2156A1', size = 2,
             alpha = .8) +
  theme_classic() +
  facet_wrap(facets = 'Year', scales = 'free_y') +
  labs(title = 'Ouest Loscolo',
       y="Mesodinium cells observed in 10 mL", x="Julian day")

### In 1 composite year

ggplot(plotMeso_OL_nozeros)+
  geom_point(aes(x=Day, y=true_count), color = '#2156A1', size = 3,
             alpha = .8) +
  theme_classic() +
  labs(title = 'Ouest Loscolo - composite year',
       y="Mesodinium cells observed in 10 mL", x="Calendar day")

### Composite year with zeros

ggplot(plotMeso_OL)+
  geom_point(aes(x=Day, y=true_count), color = '#2156A1', size = 3,
             alpha = .8) +
  theme_classic() +
  labs(title = 'Ouest Loscolo - composite year',
       y="Mesodinium cells observed in 10 mL", x="Calendar day")

########## Running the GAM ###########
#### Model and basic model checks ####

# Formulate the GAM

# We will use the Year as a random effect and for this we need it as a factor
Season_Meso_factor <- Season_Meso %>%
  mutate(Year = as_factor(Year))


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
gam.check(gam_Meso)

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
# But we can see that the confidence interval is quite small around the
# min and max of the model fit.

# Saving data to make another plot in another script
# write.csv2(gam_Meso_newdata, 
#           'Data/GAM_outputs/Mesodinium/gam_Meso_multiyear_data.csv', 
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')

### Checking the model
ModelOutputs<-data.frame(Fitted=fitted(gam_Meso),
                         Residuals=resid(gam_Meso))

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- gam_Meso$model
# then we add the values of fitted and residuals 
# (but are they in the same order as the model? -> need to check that)
qq_data <- bind_cols(qq_data, ModelOutputs)

# Plot : verify that true data matches (more or less) model fit
ggplot(qq_data)+
  geom_point(aes(x = Day, y = true_count), color = 'red') +
  geom_point(aes(x = Day, y = Fitted), color = 'blue') +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  theme_classic() +
  labs(y = "Mesodinium count", x = "Calendar day")

# It matches! great!

# Now for the "official" qq plot

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
  labs(title = 'qq-plot for Mesodinium GAM (2016-2022)', 
       y="Sample Quantiles",x="Theoretical Quantiles")

qqplot_custom

# Save the plot
# ggsave('Plots/GAMs/Mesodinium/qqplot_Meso_custom_4sites.tiff', dpi = 300, height = 175, width = 250,
#                units = 'mm', compression = 'lzw')

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data_reordered)+
  geom_point(aes(x=Fitted,y=Residuals, color  =Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free') +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  labs(title = 'RvF-plot for Mesodinium GAM (2016-2022)',
       y="Residuals",x="Fitted Values")

RvFplot_custom

# Save the plot
# ggsave('Plots/GAMs/Mesodinium/RvFplot_custom_4sites.tiff', dpi = 300, 
# height = 175, width = 250, units = 'mm', compression = 'lzw')

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data_reordered, aes(x = Residuals, 
                                                fill = Code_point_Libelle))+
  geom_histogram(binwidth = 1)+
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free') +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic() +
  labs(title = 'Histogram of residuals for 
Mesodinium GAM (2016-2022)',
       x='Residuals', y = 'Count')

HistRes_custom

# Save the plot
# ggsave('Plots/GAMs/Mesodinium/HistRes_custom_4sites.tiff', dpi = 300,
# height = 175, width = 250,  units = 'mm', compression = 'lzw')

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

# Save response_pred
# write.csv2(response_pred, 'Data/GAM_outputs/Mesodinium/response_pred_GAMMeso_20250604.csv',
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

# Saving pred_plot so we don't have to re-run the model every time
# write.csv2(pred_plot, 'Data/GAM_outputs/Mesodinium/pred_plot_Meso_20250604_4sites.csv', 
# row.names = FALSE, fileEncoding = "ISO-8859-1")

# Transform values so they are expressed in the response variable units
# (by applying the inverse link function)
response_pred_plot <- as.data.frame(pred_plot) %>%
  mutate(across(c('median.fit', 'lwrS', 'uprS', 'lwrP', 'uprP'), ~ ilink(.)))

# And save that shall we
# write.csv2(response_pred_plot, 'Data/GAM_outputs/Mesodinium/response_pred_plot_Meso_20250604_4sites.csv',
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

# plotting
ggplot(response_pred_plot, aes(x = Day)) +
  geom_line(aes(y = median.fit), lwd = 1, color = 'black') +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2, fill = "red") +
  geom_ribbon(aes(ymin = lwrP, ymax = uprP), alpha = 0.2, fill = "red") +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  labs(y = "Mesodinium (cells.L-1)",
       x = "Julian day")


###### Plot aesthetics #########

# Import 'response_pred_plot' if necessary
# response_pred_plot <- read.csv2('Data/GAM_outputs/Mesodinium/response_pred_plot_Meso_20250604_4sites.csv',
#                                 header = TRUE, fileEncoding = 'ISO-8859-1')

# And Season_Meso
# Season_Meso <- read.csv2('Data/REPHY_outputs/Season_Meso_20250604.csv',
#                          header = TRUE, fileEncoding = 'ISO-8859-1')

# Reordering the factor 'Code_point_Libelle' so that the sites appear in the
# plot in the desired order
response_pred_plot <- response_pred_plot %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# And do the exact same thing in the Season_Dino dataset
Season_Meso <- Season_Meso %>%
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


##### Plot this #####
ggplot(response_pred_plot, aes(x = Day, y = median.fit, 
                               color = Code_point_Libelle,
                               fill = Code_point_Libelle)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
  geom_path(lwd = 1) +
  # scale_y_continuous(limits = c(0,40)) +
  geom_point(data = Season_Meso_crop, aes(x = Day, y = true_count), 
             size = .8, alpha = .5) +
  # Add a "ghost point" to force the minimum y-axis range to 5
  geom_point(aes(x = 1, y = 5), color = 'transparent', fill = 'transparent',
             size = .8, alpha = .5) +
  labs(y = "Mesodinium cells observed",
       x = "Julian day",
       title = "Poisson GAM of Mesodinium phenology"
  ) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic()

# Saving plot
# ggsave('Plots/GAMs/gam_Meso_4sites_crop.tiff', dpi = 300, height = 175, width = 250,
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
  geom_point(data = Season_Meso_crop, aes(x = Day, y = Mesodinium_genus), 
             size = .8, alpha = .5) +
  # Add a "ghost point" to force the minimum y-axis range to 5
  geom_point(aes(x = 1, y = 500), color = 'transparent', fill = 'transparent',
             size = .8, alpha = .5) +
  # Labels
  labs(y = c(expression(paste("Mesodinium cells.L"^'-1'))),
       x = "Julian day",
       title = "Poisson GAM of Mesodinium phenology") +
  # facets
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  theme_classic()

# Saving plot
# ggsave('Plots/GAMs/Mesodinium/gam_Meso_4sites_crop_cells_per_L.tiff', 
# dpi = 300, height = 175, width = 200,
#        units = 'mm', compression = 'lzw')