###### GAM for Mesodinium phenology with mgcv - unified ###
## V. POCHIC
# 2024-04-18

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
  filter(Year >= 2007)

#### More zeros ####

### The goal here is to create a table where all absences of detection are 
# listed as zeros
Table_Meso_zeros <- Table_phyto_OL %>%
  pivot_wider(names_from = Taxon, values_from = Comptage, values_fill = 0) %>%
  select(starts_with('Mesodinium') | 
           contains(c('Code.Region', 'Code_point_Libelle', 'Year', 'Month',
                      'Date', 'Code.parametre', 'SALI', 'TEMP')))

Table_Meso_zeros <- Table_Meso_zeros %>%
  # Getting rid of some unwelcome guests (who got onboard because the taxon
  # name contains 'sali')
  select(-(contains('Diplo')))

# Save the table
# write.csv2(Table_Meso_zeros, 'Table_Meso_zeros.csv', row.names = FALSE,
#            fileEncoding = 'ISO-8859-1')

### Clean the environment
rm(Table_phyto_taxon)
rm(Table_phyto_OL)

########## Organizing the data #######
#### Seasonality ####

### Create a seasonality dataset (with the calendar day/julian day variable)
Season_Meso <- Table_Meso_zeros %>%
  # Only FLORTOT
  filter(Code.parametre == 'FLORTOT') %>%
  # create a calendar day variable
  mutate(Day = as.numeric(yday(Date))) %>%
  mutate(Date = ymd(Date)) %>%#### Note that for the 2 sites in Arcachon, counts are sometimes done in
  # 100 mL, so we need to filter out those counts
  # For this, we go and look for any count values IN ANY DINOPHYSIS TAXON that
  # cannot correspond to 10 mL counts (i.e., not multiples of 100)
  ### We replace these values with 0. Note that for a given date, we can have
  # one Dinophysis species counted in 10 mL and another in 100 mL
  ## That's why we have to proceed taxon by taxon
  mutate(across(.cols = c('Mesodinium', 'Mesodinium rubrum'), 
                .fns = ~ ifelse(. %% 100 != 0, 0, .))) %>%
  # And creating some count variables for Mesodinium as a genus
  mutate(Mesodinium_genus = rowSums(across(contains('Mesodinium')))) %>%
  # create the log of abundance + 1
  mutate(log_c = log10(Mesodinium_genus+1)) %>%
  # create a 'true count' variable
  # this variable corresponds to the number of cells that were actually
  # counted by the operator (in 10 mL). This variable will follow a Poisson
  # distribution, contrary to Dinophysis_genus because the conversion from
  # 10 mL to 1L (*100) prevents some intermediate values (e.g., 150 cells.L-1)
  mutate(true_count = Mesodinium_genus/100) %>%
  # Filter out exceptionnaly high counts (>500 cells observed in 10 mL)
  # this represents 3 events (2 in Antifer and 1 in Cabourg)
  # filter(true_count < 500) %>%
  # converting the site to factor for the model
  mutate(Code_point_Libelle = as.factor(Code_point_Libelle))


# Checking the homogeneity of the sampling between years
# We make a histogram of samples for the 17 years (2006-2022)

hist(Season_Meso$Year, breaks = 16)
# Approximately homogeneous around 25 per year, a bit more sampling in 2006
# (almost 35), a bit less in 2020 (15-20) due to Covid.

# Homogeneity across the year, for the 26 fortnights in a year 
# (fortnightly sampling)
hist(Season_Meso$Day, breaks = 26)
# Approximately OK, except for the last fortnight of the year 
# (Christmas + NY's eve)

# This seems quite fine!!!

#### Where the zeros lie ####

# Distribution of data for the cell count of Mesodinium
hist(Season_Meso$Mesodinium_genus, breaks = 15)


# According to this distribution, the data seem to follow a Poisson
# distribution that is zero-inflated

# We want a distribution of the zeros to estimate when (in the 'example year')
# they have the most chance of being true or false zeros

Season_zeros <- Season_Meso %>%
  filter(Mesodinium_genus == 0)

hist(Season_zeros$Day, breaks = 26)
# There we have it : a histogram of distribution of 0 values for Mesodinium
# for 26 periods of 15 days in an "example year".

# We can note that it is approximately the opposite of the Mesodinium
# seasonality

# To deal with the zero-inflation, we want to reduce the weight of the zeros
# when we suspect they might be false zeros.

# Because the sampling protocol has a rather poor resolution, sometimes 
# Mesodinium might be counted as absent when it is in fact present.
# Additionnaly, the quantification limit (100 cells.L-1) is quite high, this
# might favor false-zero observations.
# Finally, Mesodinium is a 'rare' species (almost always low cell counts), so
# it is particularly suject to these sampling biases.

# What we can do is attribute each zero-value a weight depending on the time
# bin they belong to.
# But the right way to do it is to calculate the chance of having a zero given
# the frequency of sampling in a given bin (=15-day period)
# (e.g. few zeros in the Christmas period, because few samplings)

# Create a 'Fortnight' variable so we divide the data in 26 bins
Where_zeros <- Season_Meso %>%
  mutate(Week = week(Date)) %>%
  # 'ceiling' takes the upper integer of a decimal number
  mutate(Fortnight = ceiling(Week/2)) %>%
  arrange(Fortnight)

# calculate, for each fortnight, the ratio of zeros/samplings
Where_zeros <- Where_zeros %>%
  group_by(Fortnight) %>%
  # Create a vector that 'brands' the zeros (absence of Mesodinium)
  # with 1 (0 for other values)
  mutate(is.zero = ifelse(Mesodinium_genus == 0, 1, 0)) %>%
  # Count the zeros and all the rows (= all samples) in each fortnight
  summarise(zeros = sum(is.zero), samplings = n(), 
            .groups = 'keep') %>%
  mutate(ratio_zeros = zeros/samplings)

# Plot of the ratio of zeros
plot(x = Where_zeros$Fortnight, y = Where_zeros$ratio_zeros)

# We can consider the likelihood of a zero being a true zero based on this
# distribution of zeros in a (composite) year.
# The more zeros are observed in a fortnight, the more likely they are to be
# 'true' zeros.

# Let's now associate the weights to the data. For now we only add differential
# weights to the zeros, and put a weight of 1 on positive abundances
Season_Meso <- Season_Meso %>%
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2))

# Join the dataset of zeros
Season_Meso <- left_join(Season_Meso, Where_zeros, by = 'Fortnight', 
                         suffix = c('','')) %>%
  # And add a weights vector that is equal to ratio_zeros if 0, and 1 else
  mutate(weight = ifelse(Mesodinium_genus == 0, ratio_zeros, 1)) %>%
  # And try the same with stronger weights (ratio_zeros squared)
  mutate(weight_squared = ifelse(Mesodinium_genus == 0, ratio_zeros^2, 1))

# Clean up
rm(Where_zeros)
rm(Season_zeros)

### Now in Season_Meso, we have different weights for zero counts, depending
# on our 'confidence' in them being true zeros


###### Detecting outliers in the time series ######

# Zero-inflation is not the only problem, we also have some extreme values
# (outliers)

# We do not want to simply remove them, as these are TRUE observations, that
# describe an extreme event that did occur (they are not measurement errors)
# However, we do not want them to skew the phenology, as these extreme values
# can apply a strong leverage effect to the model.

# We want to attribute a weight to outliers that is a function of their distance
# to the median cell count (that is, function of their 'unlikeliness' or their
# exceptional nature)

# We apply the MAD (Median Absolute Deviation) method to detect outliers
outliers <- function(x) {
  median_dev <- abs(x - median(x, na.rm = T))
  indexes <- which(median_dev > (2.323561*1.486) * mad(x,na.rm = T))
  return(list(indexes, median_dev))
}

# Applying the function
outliers(Season_Meso$Mesodinium_genus)
# Because of the prevalence of zeros in the dataset, the function returns as an
# outlier anything that is not zero.
# Let's try it again on the dataset without zeros

Season_Meso_nozero <- filter(Season_Meso, Mesodinium_genus != 0)
outliers <- outliers(Season_Meso_nozero$Mesodinium_genus) 
# Nice :)
# The function returns 19 potential outliers out of 162 positive Mesodinium
# cell counts (around 12%)
# Let's check how they look like
# Isolate outliers in the dataset
Season_Meso_outliers <- filter(Season_Meso_nozero, row_number() %in% outliers[[1]])

# Minimum cell count for an outlier
min(Season_Meso_outliers$Mesodinium_genus)
# Maximum cell count for an outlier
max(Season_Meso_outliers$Mesodinium_genus)

# Let's apply weights to outliers
Season_Meso_nozero <- Season_Meso_nozero %>%
  # Caclulate the median deviation (how far a point is from the median of the 
  # dataset?)
  mutate(median_dev = abs(Season_Meso_nozero$Mesodinium_genus 
                          - median(Season_Meso_nozero$Mesodinium_genus, na.rm = T))) %>%
  # Flag outliers with a dummy variable (1 = outlier, 0 = not outlier)
  mutate(outlier = ifelse(median_dev > (2.323561*1.486) 
                          * mad(Season_Meso_nozero$Mesodinium_genus,na.rm = T),
                          1, 0))

# Distribution of median deviations (from 0 to max)
hist(Season_Meso_nozero$median_dev)
# Max median deviation in the dataset
max(Season_Meso_nozero$median_dev)
# We got a distribution, now we would like to scale values between 0 and 1 to
# apply weights

# Then, scale between 0 and 1
Season_Meso_weights <- Season_Meso_nozero %>%
  # We apply a weight to each outlier that corresponds to :
  # the median of the series divided by the outlier count's median deviation.
  # This way, the furthest an outlier is from the median, the lowest weight it 
  # will have.
  mutate(outlier_weight = ifelse(
    # Condition
    outlier == 1,
    # If yes
    median(Season_Meso_nozero$Mesodinium_genus)/median_dev,
    # Else (non-outliers have a weight of 1)
    1))

# Distribution of outlier weight
hist(Season_Meso_weights$outlier_weight)
#The vast majority have a weight of 1 or close to 1


# Let's check the weights for detected outliers
test <- filter(Season_Meso_weights, outlier == 1)
min(test$outlier_weight)
# The minimum is not 0, which is nice (we don't want to exclude any point)
max(test$outlier_weight)
# Between 5 and 25% of the normal weight (1). Seems reasonable.

# Now incorporating the outlier weights in the dataset

# Select only relevant columns in Season_Meso_weights
Season_Meso_weights <-  select(Season_Meso_weights, c('Date', 'outlier', 'outlier_weight'))

# Then join datasets (with weights for zeros and for outliers)
Season_Meso_joined <- left_join(Season_Meso, Season_Meso_weights, by = 'Date') %>%
  # NAs were introduced for non-matching values (zero counts), get rid of them
  mutate(outlier = ifelse(is.na(outlier), 0, outlier)) %>%
  # And apply a "unified weight" that encompasses both weights
  mutate(unif_weight = ifelse(outlier != 1, weight_squared, outlier_weight))

# Clean up
rm(Season_Meso_weights)
rm(Season_Meso_nozero)
rm(Season_Meso_outliers)
rm(outliers)
rm(test)

########## Running the GAM ###########
#### Model and basic model checks ####
# Formulate the GAM

Season_Meso <- filter(Season_Meso, Code_point_Libelle != 'Parc Leucate 2'
                       & Code_point_Libelle != 'Sète mer')
                      # Only zeros in these sites (not suitable for a Poisson
                      # distribution)

gam_Meso <- gam(data = Season_Meso, 
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

summary(gam_Meso)
gam.check(gam_Meso)

# gam.check indicates that there might be an issue with the random effect
# smooth for 'Year' (low p-value)

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
                                        se.fit = TRUE, type = 'link')[1:2]),
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
# min and max of the model fit.

# Plot

ggplot(gam_Meso_newdata, aes(x = Day, y = fit_resp))+
  geom_ribbon(aes(x=Day, ymin=right_lwr, ymax=right_upr), fill = 'grey70', alpha=0.7) +
  geom_line(linewidth = 1) +
  geom_point(data = Season_Meso, aes(x = Day, y = true_count, color = Year), shape = 21) +
  scale_color_viridis_c('Year') +
  facet_wrap(facets = c('Code_point_Libelle')) +
  theme_classic()

# Saving plot
# ggsave('gam_Meso_OL_all.tiff', dpi = 300, height = 120, width = 160, 
#        units = 'mm', compression = 'lzw')

### Checking the model
ModelOutputs<-data.frame(Fitted=fitted(gam_Meso),
                         Residuals=resid(gam_Meso))

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
N <- 10000

# N draws from a multivariate normal distributed matrix with mean 0 (?)
BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)

# Calculating a function now (?)
Cg <- predict(gam_Meso, newd, type = "lpmatrix")
simDev <- Cg %*% t(BUdiff)

# Remove useless things to free memory space
rm(Cg)
rm(BUdiff)
rm(Vb)

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

# Now that crit is calculated, remove the big tables that take all the space
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
# write.csv2(pred_plot, 'pred_plot_MESO_20240418.csv', row.names = FALSE,
#            fileEncoding = "ISO-8859-1")

# Transform values so they are expressed in the response variable units
# (by applying the inverse link function)
response_pred_plot <- as.data.frame(pred_plot) %>%
  mutate(across(c('median.fit', 'lwrS', 'uprS', 'lwrP', 'uprP'), ~ ilink(.)))

# Let's save the response_pred_plot (might be useful later)
# write.csv2(response_pred_plot, 'response_pred_plot_MESO_20240424.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')

# You can check that the backtransformation by the inverse link function
# worked by comparing the pred_plot and response_pred_plot tables.
# The response_pred_plot values should be higher.

# Checking the confidence intervals
# minimum values
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
  labs(y = "Mesodinium (cells.L-1)",
       x = "Day of the year")

# Draws from the bayesian posterior distribution to be plotted on the graph :
# not optimized yet

# # Drawing a sample from the posterior distribution of the model
# sims <- rmvn(N, mu = coef(gam_Meso), sig = Vb)
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
pheno_palette10 <- c('red3', 'orangered1', 'dodgerblue4', 'dodgerblue1', 
                     'chartreuse4', 'chartreuse2','goldenrod3', 'darkgoldenrod1',
                     'darkorchid1', 'deeppink2'
                     )

# Reordering the factor 'Code_point_Libelle' so that the sites appear in the
# plot in the desired order
response_pred_plot <- response_pred_plot %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Bouzigues (a)', 'Diana centre'))

# And do the exact same thing in the Season_Meso dataset
Season_Meso <- Season_Meso %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Bouzigues (a)', 'Diana centre'))

#### Cropping out data points for plotting ####
# Crop out highest data points of each site to have better view of the model
Season_Meso_crop <- Season_Meso %>%
  # Antifer
  filter(ifelse(true_count > 10 & Code_point_Libelle == 'Antifer ponton pétrolier',
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
  filter(ifelse(true_count > 90 & Code_point_Libelle == 'Men er Roue',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Ouest Loscolo
  filter(ifelse(true_count > 90 & Code_point_Libelle == 'Ouest Loscolo',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Le Cornard
  filter(ifelse(true_count > 90 & Code_point_Libelle == 'Le Cornard',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Auger
  filter(ifelse(true_count > 90 & Code_point_Libelle == 'Auger',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Arcachon
  filter(ifelse(true_count > 25 & Code_point_Libelle == 'Arcachon - Bouée 7',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  
  # Teychan
  filter(ifelse(true_count > 35 & Code_point_Libelle == 'Teychan bis',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Bouzigues
  filter(ifelse(true_count > 90 & Code_point_Libelle == 'Bouzigues (a)',
                # if condition met, drop the line
                FALSE,
                # else, keep the line
                TRUE)) %>%
  # Diana
  filter(ifelse(true_count > 90 & Code_point_Libelle == 'Diana centre',
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
  geom_point(data = Season_Meso_crop, aes(x = Day, y = true_count), 
             size = .8, alpha = .5) +
  labs(y = "Mesodinium cells observed in 10 mL",
       x = "Day of the year",
       title = "Poisson GAM of Mesodinium phenology" #,
       #subtitle = sprintf("Each line is one of %i draws from the Bayesian posterior distribution of the model", nrnd)
  ) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette10, guide = 'none') +
  scale_fill_discrete(type = pheno_palette10, guide = 'none') +
  theme_classic()

# Saving plot
# ggsave('gam_Meso_allsites_crop.tiff', dpi = 300, height = 175, width = 250,
#        units = 'mm', compression = 'lzw')

# Saving plot
# ggsave('full_weighted_gam_Meso_OL.tiff', height = 150, width = 225, units = 'mm', compression = 'lzw')


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
