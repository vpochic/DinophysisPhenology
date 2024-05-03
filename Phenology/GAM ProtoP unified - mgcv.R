###### GAM for Protoperidinium phenology with mgcv - unified 2 ###
## V. POCHIC
# 2024-05-02


#### Packages and functions ####
library(tidyverse)
library(ggplot2)
library(mgcv)
library(viridis)
library(nlme)

### Additional stuff (deriv function notably)

source('derivFun.R')
source('tsDiagGamm.R')

#### Import data ####

# Table with phytoplankton counts
Table_phyto_taxon <- read.csv2('Table1_phyto_taxon.csv', fileEncoding = "ISO-8859-1")

# Select only stations of interest from 2006 --> 2022
Table_phyto_OL <- filter(Table_phyto_taxon, Code_point_Libelle %in% 
                           c(# Baie de Seine
                             'Antifer ponton pétrolier', 'Cabourg',
                             # Bretagne Sud
                             'Men er Roue', 'Ouest Loscolo',
                             # Pertuis charentais
                             'Auger', 'Le Cornard',
                             # Arcachon
                             'Arcachon - Bouée 7', 'Teychan bis',
                             # Mediterranée
                             'Parc Leucate 2', 'Bouzigues (a)', 'Sète mer',
                             'Diana centre')
                         ) %>%
  filter(Year >= 2007 & Year <= 2022)

#### More zeros ####

### The goal here is to create a table where all absences of detection are 
# listed as zeros
Table_ProtoP_zeros <- Table_phyto_OL %>%
  pivot_wider(names_from = Taxon, values_from = Comptage, values_fill = 0) %>%
  select(starts_with('Protoperidinium') | 
           contains(c('Code.Region', 'Code_point_Libelle', 'Year', 'Month',
                      'Date', 'Code.parametre', 'SALI', 'TEMP')))

Table_ProtoP_zeros <- Table_ProtoP_zeros %>%
  # Getting rid of some unwelcome guests (who got onboard because the taxon
  # name contains 'sali')
  select(-(contains('Diplo')))

# Save Table_ProtoP_zeros so we don't have to re-run the model to plot again
# write.csv2(Table_ProtoP_zeros, 'Table_ProtoP_zeros.csv', row.names = FALSE,
# fileEncoding = "ISO-8859-1")


### Clean the environment
rm(Table_phyto_taxon)
rm(Table_phyto_OL)

########## Organizing the data #######
#### Seasonality ####

### Create a seasonality dataset (with the calendar day/julian day variable)
Season_ProtoP <- Table_ProtoP_zeros %>%
  # Only FLORTOT
  filter(Code.parametre == 'FLORTOT') %>%
  # create a calendar day variable
  mutate(Day = as.numeric(yday(Date))) %>%
  mutate(Date = ymd(Date)) %>%
  #### Note that for the 2 sites in Arcachon, counts are sometimes done in
  # 100 mL, so we need to filter out those counts
  # For this, we go and look for any count values IN ANY DINOPHYSIS TAXON that
  # cannot correspond to 10 mL counts (i.e., not multiples of 100)
  ### We replace these values with 0. Note that for a given date, we can have
  # one Dinophysis species counted in 10 mL and another in 100 mL
  ## That's why we have to proceed taxon by taxon
  mutate(across(.cols = contains('Protoperidinium'), .fns = ~ ifelse(. %% 100 != 0, 0, .))) %>%
  # And creating some count variables for Dinophysis as a genus
  mutate(Protoperidinium_genus = rowSums(across(contains('Protoperidinium')))) %>%
  # create the log of abundance + 1
  mutate(log_c = log10(Protoperidinium_genus+1)) %>%
  # create a 'true count' variable
  # this variable corresponds to the number of cells that were actually
  # counted by the operator (in 10 mL). This variable will follow a Poisson
  # distribution, contrary to Protoperidinium_genus because the conversion from
  # 10 mL to 1L (*100) prevents some intermediate values (e.g., 150 cells.L-1)
  mutate(true_count = Protoperidinium_genus/100) %>%
  # Filter out exceptionnaly high counts (>500 cells observed in 10 mL)
  # this represents 3 events (2 in Antifer and 1 in Cabourg)
  filter(true_count < 500) %>%
  # converting the site to factor for the model
  mutate(Code_point_Libelle = as.factor(Code_point_Libelle))

# Checking the homogeneity of the sampling between years
# We make a histogram of samples for the 17 years (2006-2022)

hist(Season_ProtoP$Year, breaks = 16)
# Approximately homogeneous around 25 per year, a bit more sampling in 2006
# (almost 35), a bit less in 2020 (15-20) due to Covid.

# Homogeneity across the year, for the 26 fortnights in a year 
# (fortnightly sampling)
hist(Season_ProtoP$Day, breaks = 26)
# Approximately OK, except for the last fortnight of the year 
# (Christmas + NY's eve)

# This seems quite fine!!!

########## Running the GAM ###########
#### Model and basic model checks ####

# Formulate the GAM

gam_ProtoP <- gam(data = Season_ProtoP, 
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

summary(gam_ProtoP)
gam.check(gam_ProtoP)

# gam.check indicates that there might be an issue with the random effect
# smooth for 'Year' (low p-value)

# Create a new 'empty' dataset for storing model prediction (fit)
gam_ProtoP_newdata <- expand_grid(Day=seq(1, 365),
                                # We add a Year vector as it has become a factor
                                # of the model
                                Year=seq(min(Season_ProtoP$Year), 
                                        max(Season_ProtoP$Year)),
                                # We also add the site data
                                Code_point_Libelle = unique(Season_ProtoP$Code_point_Libelle))

## Get the inverse link function of the model
# With this function, we can transform the prediction we make on the link
#scale to the response scale
ilink <- gam_ProtoP$family$linkinv

## Predict : add fit and se.fit on the **link** scale
gam_ProtoP_newdata <- bind_cols(gam_ProtoP_newdata, 
                              setNames(as_tibble(
                                predict(gam_ProtoP, gam_ProtoP_newdata, 
                                        se.fit = TRUE, type = 'link')[1:2]),
                                c('fit_link','se_link')))

## Create the 95% confidence interval (2*standard error fit) and backtransform 
# to response variable using the inverse link function
gam_ProtoP_newdata <- mutate(gam_ProtoP_newdata,
                           fit_resp  = ilink(fit_link),
                           right_upr = ilink(fit_link + (2 * se_link)),
                           right_lwr = ilink(fit_link - (2 * se_link)))
# Check the confidence interval. It should not extend below 0 (negative counts
# are impossible)
min(gam_ProtoP_newdata$right_lwr) # Nice :)
min(gam_ProtoP_newdata$fit_resp)
max(gam_ProtoP_newdata$right_upr) # Nice too
max(gam_ProtoP_newdata$fit_resp)
# But we can see that the confidence interval is ridiculously small around the
# min and max of the model fit.

# Plot

ggplot(gam_ProtoP_newdata, aes(x = Day, y = fit_resp))+
  geom_ribbon(aes(x=Day, ymin=right_lwr, ymax=right_upr), fill = 'grey70', alpha=0.7) +
  geom_line(linewidth = 1) +
  geom_point(data = Season_ProtoP, aes(x = Day, y = true_count, color = Year), shape = 21) +
  scale_color_viridis_c('Year') +
  facet_wrap(facets = c('Code_point_Libelle')) +
  theme_classic()

# Saving plot
# ggsave('gam_Dino_OL_all.tiff', dpi = 300, height = 120, width = 160, 
#        units = 'mm', compression = 'lzw')

### Checking the model
ModelOutputs<-data.frame(Fitted=fitted(gam_ProtoP),
                         Residuals=resid(gam_ProtoP))

# P3 is residuals vs fitted
p3<-ggplot(ModelOutputs)+
  geom_point(aes(x=Fitted,y=Residuals))+
  theme_classic()+
  labs(y="Residuals",x="Fitted Values")

# P4 is the qq-plot
p4<-ggplot(ModelOutputs) +
  stat_qq(aes(sample=Residuals))+
  stat_qq_line(aes(sample=Residuals))+
  theme_classic()+
  labs(y="Sample Quartiles",x="Theoretical Quartiles")

p3
# Resid. vs fitted plot is approximately ok, even if there are some structures
# due to the weights applied
p4
# QQ-plot is quite shit.

rm(p3)
rm(p4)

#### Confidence intervals ####

# We can try to obtain better confidence intervals.
# This section is almost entirely taken from a post by Gavin Simpson on his
# blog 'From the bottom of the heap'

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
Vb <- vcov(gam_ProtoP)
# New data
newd <- with(Season_ProtoP, data.frame(Day = rep(seq(1, 365, length = 365),
                                               length(unique(Year))*
                                                 length(unique(Code_point_Libelle))
                                               )))

newd <- group_by(newd, Day) %>%
  arrange(Day, by_group = TRUE)

# Add the Year vector
newd$Year <- rep(unique(Season_ProtoP$Year), 
                 365*length(unique(Season_ProtoP$Code_point_Libelle)))

newd <- group_by(newd, Year, Day) %>%
  arrange(Year, by_group = TRUE)

# And the site vector
newd$Code_point_Libelle <- rep(unique(Season_ProtoP$Code_point_Libelle), 
                               365*length(unique(Season_ProtoP$Year)))

newd <- ungroup(newd) %>%
  group_by(Code_point_Libelle) %>%
  arrange(Code_point_Libelle)

# Prediction by the model on the ***link*** scale, on the new data
pred <- predict(gam_ProtoP, newd, se.fit = TRUE, type = 'link')
# Isolate standard error of the predicted fit
se.fit <- pred$se.fit

# Set the pseudo-random seed to make results reproducible (?)
set.seed(42)
# specify the number of simulations to generate
N <- 10000

# N draws from a multivariate normal distributed matrix with mean 0 (?)
BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)

# Remove
rm(Vb)

# Calculating a function now (?)
Cg <- predict(gam_ProtoP, newd, type = "lpmatrix")
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
# write.csv2(pred_plot, 'pred_plot_20240418_2.csv', row.names = FALSE,
#            fileEncoding = "ISO-8859-1")

# Transform values so they are expressed in the response variable units
# (by applying the inverse link function)
response_pred_plot <- as.data.frame(pred_plot) %>%
  mutate(across(c('median.fit', 'lwrS', 'uprS', 'lwrP', 'uprP'), ~ ilink(.)))

# And save that shall we
# write.csv2(response_pred_plot, 'response_pred_plot_ProtoP_20240502.csv',
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
  facet_wrap(facets = c('Code_point_Libelle')) +
  labs(y = "Protoperidinium (cells.L-1)",
       x = "Day of the year")

# Draws from the bayesian posterior distribution to be plotted on the graph :
# not optimized yet

# # Drawing a sample from the posterior distribution of the model
# sims <- rmvn(N, mu = coef(gam_ProtoP), sig = Vb)
# fits <- Cg %*% t(sims)
# 
# # Draw 150 simulations randomly from the sample to figure them on the graph
# nrnd <- 150
# rnd <- sample(N, nrnd)
# stackFits <- stack(as.data.frame(fits[, rnd]))
# stackFits <- transform(stackFits, Day = rep(newd$Day, length(rnd)))
# stackFits_response <- stackFits %>%
#   mutate(values_response = ilink(values))


###### Plot aesthetics #########

# Defining a beautiful color palette
pheno_palette11 <- c('red3', 'orangered1', 'chocolate4', 'chocolate',
                     'dodgerblue4', 'dodgerblue1', 'chartreuse4', 'chartreuse2',
                     'darkorchid4', 'darkorchid1', 'deeppink2'
                     )

pheno_palette12 <- c('red3', 'orangered1', 'dodgerblue4', 'dodgerblue1', 
                     'chartreuse4', 'chartreuse2','goldenrod3', 'darkgoldenrod1',
                     'darkorchid4', 'darkorchid1', 'firebrick1', 'deeppink2'
                     )

pheno_palette13 <- c('red3', 'orangered1', 'chocolate4', 'chocolate',
                     'dodgerblue4', 'dodgerblue1', 'chartreuse4', 'chartreuse2',
                     'goldenrod3', 'darkgoldenrod1',
                     'darkorchid4', 'darkorchid1', 'deeppink2'
                     )

# Import 'response_pred_plot' if necessary
# response_pred_plot <- read.csv2('response_pred_plot_ProtoP_20240502.csv',
#                                 header = TRUE, fileEncoding = 'ISO-8859-1')

# Reordering the factor 'Code_point_Libelle' so that the sites appear in the
# plot in the desired order
response_pred_plot <- response_pred_plot %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

# And do the exact same thing in the Season_ProtoP dataset
Season_ProtoP <- Season_ProtoP %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

#### Cropping out data points for plotting ####
# Crop out highest data points of each site to have better view of the model
Season_ProtoP_crop <- Season_ProtoP %>%
  # Antifer
  filter(ifelse(true_count > 40 & Code_point_Libelle == 'Antifer ponton pétrolier',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Cabourg
  filter(ifelse(true_count > 30 & Code_point_Libelle == 'Cabourg',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Men er Roue
  filter(ifelse(true_count > 40 & Code_point_Libelle == 'Men er Roue',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Ouest Loscolo
  filter(ifelse(true_count > 70 & Code_point_Libelle == 'Ouest Loscolo',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Le Cornard
  filter(ifelse(true_count > 70 & Code_point_Libelle == 'Le Cornard',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Auger
  filter(ifelse(true_count > 70 & Code_point_Libelle == 'Auger',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Arcachon - Bouée 7
  filter(ifelse(true_count > 30 & Code_point_Libelle == 'Arcachon - Bouée 7',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Teychan
  filter(ifelse(true_count > 30 & Code_point_Libelle == 'Teychan bis',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Parc Leucate
  filter(ifelse(true_count > 20 & Code_point_Libelle == 'Parc Leucate 2',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Bouzigues
  filter(ifelse(true_count > 10 & Code_point_Libelle == 'Bouzigues (a)',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Sète mer
  filter(ifelse(true_count > 15 & Code_point_Libelle == 'Sète mer',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Diana
  filter(ifelse(true_count > 40 & Code_point_Libelle == 'Diana centre',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE))


##### Plot this #####
ggplot(response_pred_plot, aes(x = Day, y = median.fit, 
                               color = Code_point_Libelle,
                               fill = Code_point_Libelle)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
  # geom_path(data = stackFits_response, mapping = aes(y = values_response, x = Day, group= ind),
  #           alpha = 0.1, colour = "grey20") +
  geom_path(lwd = 1) +
  # scale_y_continuous(limits = c(0,40)) +
  geom_point(data = Season_ProtoP, aes(x = Day, y = true_count), 
             size = .8, alpha = .5) +
  labs(y = "Protoperidinium cells observed",
       x = "Day of the year",
       title = "Poisson GAM of Protoperidinium phenology" #,
       #subtitle = sprintf("Each line is one of %i draws from the Bayesian posterior distribution of the model", nrnd)
  ) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  scale_fill_discrete(type = pheno_palette12, guide = 'none') +
  theme_classic()

# Saving plot
# ggsave('gam_ProtoP_allsites_crop2.tiff', dpi = 300, height = 175, width = 250,
#        units = 'mm', compression = 'lzw')

# This is quite nice, but Arcachon and Teychan samples are not the same volume
# as others (10 times more), so the scale is misleading.
# We need to plot the data as cells.L to rectify this.

pred_plot_cells_per_L <- response_pred_plot %>%
  mutate(across(-c('Code_point_Libelle', 'Day'), ~ ifelse(
    # Si
    Code_point_Libelle == 'Arcachon - Bouée 7' |
      # Ou
      Code_point_Libelle == 'Teychan bis',
    # Alors
    .*10,
    # Sinon
    .*100
  )))

# New plot in cells per L

ggplot(pred_plot_cells_per_L, aes(x = Day, y = median.fit, 
                               color = Code_point_Libelle,
                               fill = Code_point_Libelle)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
  # geom_path(data = stackFits_response, mapping = aes(y = values_response, x = Day, group= ind),
  #           alpha = 0.1, colour = "grey20") +
  geom_path(lwd = 1) +
  # scale_y_continuous(limits = c(0,40)) +
  geom_point(data = Season_ProtoP_crop, aes(x = Day, y = Protoperidinium_genus), 
             size = .8, alpha = .5) +
  labs(y = c(expression(paste("Protoperidinium cells.L"^'-1'))),
       x = "Day of the year",
       title = "Poisson GAM of Protoperidinium phenology" #,
       #subtitle = sprintf("Each line is one of %i draws from the Bayesian posterior distribution of the model", nrnd)
  ) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette13, guide = 'none') +
  scale_fill_discrete(type = pheno_palette13, guide = 'none') +
  theme_classic()

# Saving plot
# ggsave('gam_ProtoP_allsites.tiff', dpi = 300, height = 270, width = 300,
#        units = 'mm', compression = 'lzw')

# Checking how many of the N simulations fit entirely in the
# 1) point-wise confidence interval
# 2) simultaneous confidence interval

# function for determining if function fits entirely in the CI
inCI <- function(x, upr, lwr) {
  all(x >= lwr & x <= upr)
}

# Does it fit in the point-wise CI ?
fitsInPCI <- apply(fits, 2L, inCI, upr = pred$uprP, lwr = pred$lwrP)
# Does it fit in the simultaneous CI ?
fitsInSCI <- apply(fits, 2L, inCI, upr = pred$uprS, lwr = pred$lwrS)

sum(fitsInPCI) / length(fitsInPCI)      # Point-wise
sum(fitsInSCI) / length(fitsInSCI)      # Simultaneous

#### DayF : plotting phenology plots so they start at Day 50 ####

# We create a variable named DayF that is a modified calendar Day :
# DayF starts at Day 50 where DayF = 1 and ends at Day 49 with DayF = 365
Season_ProtoP_DayF <- Season_ProtoP %>%
  mutate(DayF = ifelse(# If Day = 50, DayF is set to 1
    Day == 50, 1,
    # Else, if Day > 50, DayF is set to Day - 49 (so Day 51 becomes Day 2 etc.)
    ifelse(Day > 50, Day - 49,
           # Else (between Day 1 and 49), Day is set to 365-49+Day, so that Day 1 
           # becomes Day 343, and Day 49 becomes Day 365, and full cycle again
           365 - 49 + Day)))

# Same with the cropped version of the data
Season_ProtoP_DayF_crop <- Season_ProtoP_crop %>%
  mutate(DayF = ifelse(# If Day = 50, DayF is set to 1
    Day == 50, 1,
    # Else, if Day > 50, DayF is set to Day - 49 (so Day 51 becomes Day 2 etc.)
    ifelse(Day > 50, Day - 49,
           # Else (between Day 1 and 49), Day is set to 365-49+Day, so that Day 1 
           # becomes Day 343, and Day 49 becomes Day 365, and full cycle again
           365 - 49 + Day)))

# Same for repsonse_pred_plot
response_pred_plot_DayF <- response_pred_plot %>%
  mutate(DayF = ifelse(# If Day = 50, DayF is set to 1
    Day == 50, 1,
    # Else, if Day > 50, DayF is set to Day - 49 (so Day 51 becomes Day 2 etc.)
    ifelse(Day > 50, Day - 49,
           # Else (between Day 1 and 49), Day is set to 365-49+Day, so that Day 1 
           # becomes Day 343, and Day 49 becomes Day 365, and full cycle again
           365 - 49 + Day)))

# We need to specify the labels in advance.
labels_DayF <- c(50, 100, 200, 300, 365)

### Now plot this
ggplot(response_pred_plot_DayF, aes(x = DayF, y = median.fit, 
                                    color = Code_point_Libelle,
                                    fill = Code_point_Libelle)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
  # geom_path(data = stackFits_response, mapping = aes(y = values_response, x = Day, group= ind),
  #           alpha = 0.1, colour = "grey20") +
  geom_line(lwd = 1) +
  # scale_y_continuous(limits = c(0,40)) +
  geom_point(data = Season_ProtoP_DayF_crop, aes(x = DayF, y = true_count), 
             size = .8, alpha = .5) +
  labs(y = "Protoperidinium cells observed",
       x = "Day of the year",
       title = "Poisson GAM of Protoperidinium phenology" #,
       #subtitle = sprintf("Each line is one of %i draws from the Bayesian posterior distribution of the model", nrnd)
  ) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  scale_fill_discrete(type = pheno_palette12, guide = 'none') +
  # Specifiy for the plot that it needs to plot the right Day (not DayF) as x-axis
  # labels
  scale_x_continuous(breaks = c(1, 51, 151, 251, 316), labels = labels_DayF) +
  theme_classic()

# Saving plot
# ggsave('gam_ProtoP_allsites_crop_DayF.tiff', dpi = 300, height = 270, width = 300,
#        units = 'mm', compression = 'lzw')
