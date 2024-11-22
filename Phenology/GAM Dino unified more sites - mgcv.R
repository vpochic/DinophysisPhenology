###### GAM for Dinophysis phenology with mgcv - unified 2 ###
## V. POCHIC
# 2024-11-20

# In this version of the script we add the sites in Pas de Calais 
# ('Point 1 Boulogne' and 'At so') and in Northern Britanny ('les Hébihens' and 
# 'Loguivy') to include sites with almost no occurrences of Dinophysis

#### Packages and functions ####
library(tidyverse)
library(ggplot2)
library(mgcv)
library(gratia)
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
                             'Arcachon - Bouée 7', 'Teychan bis',
                             # Mediterranée
                             'Parc Leucate 2', 'Bouzigues (a)', 'Sète mer',
                             'Diana centre')
                         ) %>%
  filter(Year >= 2007 & Year <= 2022)

#### More zeros ####

### The goal here is to create a table where all absences of detection are 
# listed as zeros
Table_Dino_zeros <- Table_phyto_OL %>%
  pivot_wider(names_from = Taxon, values_from = Comptage, values_fill = 0) %>%
  # Create Latitude and Longitude variables for selection (because lon and lat
  # introduce a bunch of unwanted taxa )
  mutate(Longitude = as.numeric(lon)) %>%
  mutate(Latitude = as.numeric(lat)) %>%
  select(starts_with('Dinophysis') |
           contains(c('Code.Region', 'Code_point_Libelle', 'Year', 'Month',
                      'Date', 'Heure', 'Code.parametre', 'SALI', 'TEMP', 
                      'Longitude', 'Latitude', 'ID.interne.passage')))

Table_Dino_zeros <- Table_Dino_zeros %>%
  # Getting rid of some unwelcome guests (who got onboard because the taxon
  # name contains 'sali')
  select(-(contains('Diplo')))

# Save Table_Dino_zeros so we don't have to re-run the model to plot again
# write.csv2(Table_Dino_zeros, 'Table_Dino_zeros3.csv', row.names = FALSE,
# fileEncoding = "ISO-8859-1")


### Clean the environment
rm(Table_phyto_taxon)
rm(Table_phyto_OL)

########## Organizing the data #######
#### Seasonality ####

### Create a seasonality dataset (with the calendar day/julian day variable)
Season_Dino <- Table_Dino_zeros %>%
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
  mutate(across(.cols = c('Dinophysis', 'Dinophysis + phalacroma', 'Dinophysis acuta',
                'Dinophysis acuminata', 'Dinophysis caudata', 'Dinophysis tripos',
                'Dinophysis sacculus', 'Dinophysis fortii',
                'Dinophysis hastata + odiosa'), .fns = ~ ifelse(. %% 100 != 0, 0, .))) %>%
  # And creating some count variables for Dinophysis as a genus
  mutate(Dinophysis_genus = rowSums(across(contains('Dinophysis')))) %>%
  # create the log of abundance + 1
  mutate(log_c = log10(Dinophysis_genus+1)) %>%
  # create a 'true count' variable
  # this variable corresponds to the number of cells that were actually
  # counted by the operator (in 10 mL). This variable will follow a Poisson
  # distribution, contrary to Dinophysis_genus because the conversion from
  # 10 mL to 1L (*100) prevents some intermediate values (e.g., 150 cells.L-1)
  mutate(true_count = Dinophysis_genus/100) %>%
  # Filter out exceptionnaly high counts (>500 cells observed in 10 mL)
  # this represents 3 events (2 in Antifer and 1 in Cabourg)
  filter(true_count < 500) %>%
  # converting the site to factor for the model
  mutate(Code_point_Libelle = as.factor(Code_point_Libelle))

# Checking the homogeneity of the sampling between years
# We make a histogram of samples for the 17 years (2006-2022)

hist(Season_Dino$Year, breaks = 17)
# Approximately homogeneous around 25 per year, a bit more sampling in 2006
# (almost 35), a bit less in 2020 (15-20) due to Covid.

# Homogeneity across the year, for the 26 fortnights in a year 
# (fortnightly sampling)
hist(Season_Dino$Day, breaks = 26)
# Approximately OK, except for the last fortnight of the year 
# (Christmas + NY's eve)

# This seems quite fine!!!

# Let's save that
# write.csv2(Season_Dino, 'Season_Dino.csv', row.names = FALSE,
#            fileEncoding = "ISO-8859-1")


### Plotting the dataset for one site as an example

plotDino_OL_nozeros <- filter(Season_Dino, (Code_point_Libelle == 'Ouest Loscolo'
                                            & true_count != 0))
plotDino_OL <- filter(Season_Dino, (Code_point_Libelle == 'Ouest Loscolo'))

# Plot the Dinophysis observations as function of date (timeseries view)
ggplot(plotDino_OL_nozeros)+
  geom_point(aes(x=Date, y=true_count), color = '#2156A1',
             size = 3, alpha = .8) +
  theme_classic() +
  labs(title = 'Ouest Loscolo',
         y="Dinophysis cells observed in 10 mL", x="Date")

# Save the plot
# ggsave('plotDino_timeseries_OL.tiff', height = 150, width = 300,
# dpi = 300, unit = 'mm', compression = 'lzw')


### Plot the observations per year

ggplot(plotDino_OL_nozeros)+
  geom_point(aes(x=Day, y=true_count), color = '#2156A1', size = 2,
             alpha = .8) +
  theme_classic() +
  facet_wrap(facets = 'Year', scales = 'free_y') +
  labs(title = 'Ouest Loscolo',
       y="Dinophysis cells observed in 10 mL", x="Calendar day")

# Save the plot
# ggsave('plotDino_year_OL.tiff', height = 150, width = 300,
# dpi = 300, unit = 'mm', compression = 'lzw')

### In 1 composite year

ggplot(plotDino_OL_nozeros)+
  geom_point(aes(x=Day, y=true_count), color = '#2156A1', size = 3,
             alpha = .8) +
  theme_classic() +
  labs(title = 'Ouest Loscolo - composite year',
       y="Dinophysis cells observed in 10 mL", x="Calendar day")

# Save the plot
# ggsave('plotDino_composite-year_OL.tiff', height = 150, width = 300,
# dpi = 300, unit = 'mm', compression = 'lzw')

### Composite year with zeros

ggplot(plotDino_OL)+
  geom_point(aes(x=Day, y=true_count), color = '#2156A1', size = 3,
             alpha = .8) +
  theme_classic() +
  labs(title = 'Ouest Loscolo - composite year',
       y="Dinophysis cells observed in 10 mL", x="Calendar day")

# Save the plot
# ggsave('plotDino_composite-year_OL_zeros.tiff', height = 150, width = 300,
# dpi = 300, unit = 'mm', compression = 'lzw')

#### Where the zeros lie ####

# Distribution of data for the cell count of Dinophysis
hist(Season_Dino$Dinophysis_genus, breaks = 15)

# According to this distribution, the data seem to follow a Poisson
# distribution that is zero-inflated

# We want a distribution of the zeros to estimate when (in the 'example year')
# they have the most chance of being true or false zeros

Season_zeros <- Season_Dino %>%
  filter(Dinophysis_genus == 0)

hist(Season_zeros$Day, breaks = 26)
# There we have it : a histogram of distribution of 0 values for Dinophysis
# for 26 periods of 15 days in an "example year".

# We can note that it is approximately the opposite of the Dinophysis
# seasonality

# To deal with the zero-inflation, we want to reduce the weight of the zeros
# when we suspect they might be false zeros.

# Because the sampling protocol has a rather poor resolution, sometimes 
# Dinophysis might be counted as absent when it is in fact present.
# Additionnaly, the quantification limit (100 cells.L-1) is quite high, this
# might favor false-zero observations.
# Finally, Dinophysis is a 'rare' species (almost always low cell counts), so
# it is particularly suject to these sampling biases.

# What we can do is attribute each zero-value a weight depending on the time
# bin they belong to.
# But the right way to do it is to calculate the chance of having a zero given
# the frequency of sampling in a given bin (=15-day period)
# (e.g. few zeros in the Christmas period, because few samplings)

# Create a 'Fortnight' variable so we divide the data in 26 bins
Where_zeros <- Season_Dino %>%
  mutate(Week = week(Date)) %>%
  # 'ceiling' takes the upper integer of a decimal number
  mutate(Fortnight = ceiling(Week/2)) %>%
  arrange(Fortnight)

# calculate, for each fortnight, the ratio of zeros/samplings
Where_zeros <- Where_zeros %>%
  group_by(Fortnight) %>%
  # Create a vector that 'brands' the zeros (absence of Dinophysis)
  # with 1 (0 for other values)
  mutate(is.zero = ifelse(Dinophysis_genus == 0, 1, 0)) %>%
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
Season_Dino <- Season_Dino %>%
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2))

# Join the dataset of zeros
Season_Dino <- left_join(Season_Dino, Where_zeros, by = 'Fortnight', 
                         suffix = c('','')) %>%
  # And add a weights vector that is equal to ratio_zeros if 0, and 1 else
  mutate(weight = ifelse(Dinophysis_genus == 0, ratio_zeros, 1)) %>%
  # And try the same with stronger weights (ratio_zeros squared)
  mutate(weight_squared = ifelse(Dinophysis_genus == 0, ratio_zeros^2, 1))

# Clean up
rm(Where_zeros)
rm(Season_zeros)

### Now in Season_Dino, we have different weights for zero counts, depending
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
outliers(Season_Dino$Dinophysis_genus)
# Because of the prevalence of zeros in the dataset, the function returns as an
# outlier anything that is not zero.
# Let's try it again on the dataset without zeros

Season_Dino_nozero <- filter(Season_Dino, Dinophysis_genus != 0)
outliers <- outliers(Season_Dino_nozero$Dinophysis_genus) 
# Nice :)
# The function returns 19 potential outliers out of 162 positive Dinophysis
# cell counts (around 12%)
# Let's check how they look like
# Isolate outliers in the dataset
Season_Dino_outliers <- filter(Season_Dino_nozero, row_number() %in% outliers[[1]])

# Minimum cell count for an outlier
min(Season_Dino_outliers$Dinophysis_genus)
# Maximum cell count for an outlier
max(Season_Dino_outliers$Dinophysis_genus)

# Let's apply weights to outliers
Season_Dino_nozero <- Season_Dino_nozero %>%
  # Caclulate the median deviation (how far a point is from the median of the 
  # dataset?)
  mutate(median_dev = abs(Season_Dino_nozero$Dinophysis_genus 
                          - median(Season_Dino_nozero$Dinophysis_genus, na.rm = T))) %>%
  # Flag outliers with a dummy variable (1 = outlier, 0 = not outlier)
  mutate(outlier = ifelse(median_dev > (2.323561*1.486) 
                          * mad(Season_Dino_nozero$Dinophysis_genus,na.rm = T),
                          1, 0))

# Distribution of median deviations (from 0 to max)
hist(Season_Dino_nozero$median_dev)
# Max median deviation in the dataset
max(Season_Dino_nozero$median_dev)
# We got a distribution, now we would like to scale values between 0 and 1 to
# apply weights

# Then, scale between 0 and 1
Season_Dino_weights <- Season_Dino_nozero %>%
  # We apply a weight to each outlier that corresponds to :
  # the median of the series divided by the outlier count's median deviation.
  # This way, the furthest an outlier is from the median, the lowest weight it 
  # will have.
  mutate(outlier_weight = ifelse(
    # Condition
    outlier == 1,
    # If yes
    median(Season_Dino_nozero$Dinophysis_genus)/median_dev,
    # Else (non-outliers have a weight of 1)
    1))

# Distribution of outlier weight
hist(Season_Dino_weights$outlier_weight)
#The vast majority have a weight of 1 or close to 1


# Let's check the weights for detected outliers
test <- filter(Season_Dino_weights, outlier == 1)
min(test$outlier_weight)
# The minimum is not 0, which is nice (we don't want to exclude any point)
max(test$outlier_weight)
# Between 5 and 25% of the normal weight (1). Seems reasonable.

# Now incorporating the outlier weights in the dataset

# Select only relevant columns in Season_Dino_weights
Season_Dino_weights <-  select(Season_Dino_weights, c('Date', 'outlier', 'outlier_weight'))

# Then join datasets (with weights for zeros and for outliers)
Season_Dino_joined <- left_join(Season_Dino, Season_Dino_weights, by = 'Date') %>%
  # NAs were introduced for non-matching values (zero counts), get rid of them
  mutate(outlier = ifelse(is.na(outlier), 0, outlier)) %>%
  # And apply a "unified weight" that encompasses both weights
  mutate(unif_weight = ifelse(outlier != 1, weight_squared, outlier_weight))

# Clean up
rm(Season_Dino_weights)
rm(Season_Dino_nozero)
rm(Season_Dino_outliers)
rm(outliers)
rm(test)

########## Running the GAM ###########
#### Model and basic model checks ####

# Exclude some events in Arcachon where the counts do not correspond to
# 100 mL counts (31 events in total in Arcachon and Teychan)
# Detail :
# 22/263 of positive events in Arcachon - Bouée 7 (8.3%) and
# 9/167 of positive events in Teychan bis (5.3%)
Season_Dino <- filter(Season_Dino, Dinophysis_genus %% 10 == 0)
# This command ensures that all values of Dinophysis_genus are multiples of 10
# (therefore, that all true_count values are integers)
### With the modification of the code (2024/04/17), it doesn't remove anything
# (this step is done earlier)

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

# gam.check indicates that there might be an issue with the random effect
# smooth for 'Year' (low p-value)

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
# min and max of the model fit.

# Plot

ggplot(gam_Dino_newdata, aes(x = Day, y = fit_resp))+
  geom_ribbon(aes(x=Day, ymin=right_lwr, ymax=right_upr), fill = 'grey70', alpha=0.7) +
  geom_line(linewidth = 1) +
  geom_point(data = Season_Dino, aes(x = Day, y = true_count, color = Year), shape = 21) +
  scale_color_viridis_c('Year') +
  facet_wrap(facets = c('Code_point_Libelle')) +
  theme_classic()

# Saving plot
# ggsave('gam_Dino_OL_all.tiff', dpi = 300, height = 120, width = 160, 
#        units = 'mm', compression = 'lzw')

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
  facet_wrap(facets = c('Code_point_Libelle')) +
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

# Save the plot
# ggsave('qqplot_custom_16sites.tiff', dpi = 300, height = 175, width = 250,
#                units = 'mm', compression = 'lzw')

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data_reordered)+
  geom_point(aes(x=Fitted,y=Residuals, color  =Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free') +
  scale_color_discrete(type = pheno_palette16, guide = 'none') +
  theme_classic() +
  labs(y="Residuals",x="Fitted Values")

RvFplot_custom

# Save the plot
# ggsave('RvFplot_custom_16sites.tiff', dpi = 300, height = 175, width = 250,
#                units = 'mm', compression = 'lzw')

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data_reordered, aes(x = Residuals, 
                                                fill = Code_point_Libelle))+
  geom_histogram(binwidth = 1)+
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free') +
  scale_fill_discrete(type = pheno_palette16, guide = 'none') +
  theme_classic() +
  labs(x='Residuals', y = 'Count')

HistRes_custom

# Save the plot
# ggsave('HistRes_custom_16sites.tiff', dpi = 300, height = 175, width = 250,
#                units = 'mm', compression = 'lzw')

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
# write.csv2(pred_plot, 'pred_plot_20241120_16sites.csv', row.names = FALSE,
#            fileEncoding = "ISO-8859-1")

# Transform values so they are expressed in the response variable units
# (by applying the inverse link function)
response_pred_plot <- as.data.frame(pred_plot) %>%
  mutate(across(c('median.fit', 'lwrS', 'uprS', 'lwrP', 'uprP'), ~ ilink(.)))

# And save that shall we
# write.csv2(response_pred_plot, 'response_pred_plot_20241120_16sites.csv',
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
  labs(y = "Dinophysis (cells.L-1)",
       x = "Day of the year")

# Draws from the bayesian posterior distribution to be plotted on the graph :
# not optimized yet

# # Drawing a sample from the posterior distribution of the model
# sims <- rmvn(N, mu = coef(gam_Dino), sig = Vb)
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
pheno_palette16 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646',
                     '#642C3A', '#DEB1CC', '#FC4D6B', '#791D40')

# Import 'response_pred_plot' if necessary
response_pred_plot <- read.csv2('response_pred_plot_20241120_16sites.csv',
                                header = TRUE, fileEncoding = 'ISO-8859-1')

# And Season_Dino
Season_Dino <- read.csv2('Season_Dino.csv',
                         header = TRUE, fileEncoding = 'ISO-8859-1')

# Reordering the factor 'Code_point_Libelle' so that the sites appear in the
# plot in the desired order
response_pred_plot <- response_pred_plot %>%
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


##### Plot this #####
ggplot(response_pred_plot, aes(x = Day, y = median.fit, 
                               color = Code_point_Libelle,
                               fill = Code_point_Libelle)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
  # geom_path(data = stackFits_response, mapping = aes(y = values_response, x = Day, group= ind),
  #           alpha = 0.1, colour = "grey20") +
  geom_path(lwd = 1) +
  # scale_y_continuous(limits = c(0,40)) +
  geom_point(data = Season_Dino_crop, aes(x = Day, y = true_count), 
             size = .8, alpha = .5) +
  # Add a "ghost point" to force the minimum y-axis range to 5
  geom_point(aes(x = 1, y = 5), color = 'transparent', fill = 'transparent',
             size = .8, alpha = .5) +
  labs(y = "Dinophysis cells observed",
       x = "Day of the year",
       title = "Poisson GAM of Dinophysis phenology" #,
       #subtitle = sprintf("Each line is one of %i draws from the Bayesian posterior distribution of the model", nrnd)
  ) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette16, guide = 'none') +
  scale_fill_discrete(type = pheno_palette16, guide = 'none') +
  theme_classic()

# Saving plot
# ggsave('gam_Dino_16sites_crop.tiff', dpi = 300, height = 175, width = 250,
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
  geom_point(data = Season_Dino_crop, aes(x = Day, y = Dinophysis_genus), 
             size = .8, alpha = .5) +
  labs(y = c(expression(paste("Dinophysis cells.L"^'-1'))),
       x = "Day of the year",
       title = "Poisson GAM of Dinophysis phenology" #,
       #subtitle = sprintf("Each line is one of %i draws from the Bayesian posterior distribution of the model", nrnd)
  ) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette13, guide = 'none') +
  scale_fill_discrete(type = pheno_palette13, guide = 'none') +
  theme_classic()

# Saving plot
ggsave('gam_Dino_allsites_crop_cells_per_L.tiff', dpi = 300, height = 270, width = 300,
       units = 'mm', compression = 'lzw')

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
Season_Dino_DayF <- Season_Dino %>%
  mutate(DayF = ifelse(# If Day = 50, DayF is set to 1
    Day == 50, 1,
    # Else, if Day > 50, DayF is set to Day - 49 (so Day 51 becomes Day 2 etc.)
    ifelse(Day > 50, Day - 49,
           # Else (between Day 1 and 49), Day is set to 365-49+Day, so that Day 1 
           # becomes Day 343, and Day 49 becomes Day 365, and full cycle again
           365 - 49 + Day)))

# Same with the cropped version of the data
Season_Dino_DayF_crop <- Season_Dino_crop %>%
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
  geom_point(data = Season_Dino_DayF_crop, aes(x = DayF, y = true_count), 
             size = .8, alpha = .5) +
  labs(y = "Dinophysis cells observed",
       x = "Day of the year",
       title = "Poisson GAM of Dinophysis phenology" #,
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
# ggsave('gam_Dino_allsites_crop_DayF.tiff', dpi = 300, height = 270, width = 300,
#        units = 'mm', compression = 'lzw')

#### Alternative GAM for teaching purposes ####

gam_Dino_OL <- gam(data = plotDino_OL, 
                   # Only a spline for the day of the year
                   # The use of a cyclic basis spline helps to make ends meet at the
                   # first and last days of the year
                   # 'k = -1' allows the model to fix the 'best' number of basic
                   # functions (= knots)
                   formula = true_count~s(Day, bs = 'cc', k = -1) 
                   # We add the Year as a random effect. This will help to assess and
                   # smooth the effects of interannual variability in the phenology
                   + s(Year, bs = 're', k = -1),
                   # Introducing the weights
                   # weights = unif_weight,
                   # Using a Poisson distribution for count data
                   family = poisson(),
                   # Restricted maximum likelihood estimation (recommended method)
                   method = 'REML')

summary(gam_Dino_OL)
gam.check(gam_Dino_OL)

# gam.check indicates that there might be an issue with the random effect
# smooth for 'Year' (low p-value)

# Create a new 'empty' dataset for storing model prediction (fit)
gam_Dino_OL_newdata <- expand_grid(Day=seq(1, 365),
                                   # We add a Year vector as it has become a factor
                                   # of the model
                                   Year=seq(min(Season_Dino$Year), 
                                            max(Season_Dino$Year)))

## Get the inverse link function of the model
# With this function, we can transform the prediction we make on the link
#scale to the response scale
ilink <- gam_Dino_OL$family$linkinv

## Predict : add fit and se.fit on the **link** scale
gam_Dino_OL_newdata <- bind_cols(gam_Dino_OL_newdata, 
                                 setNames(as_tibble(
                                   predict(gam_Dino_OL, gam_Dino_OL_newdata, 
                                           se.fit = TRUE, type = 'link')[1:2]),
                                   c('fit_link','se_link')))

## Create the 95% confidence interval (2*standard error fit) and backtransform 
# to response variable using the inverse link function
gam_Dino_OL_newdata <- mutate(gam_Dino_OL_newdata,
                              fit_resp  = ilink(fit_link),
                              right_upr = ilink(fit_link + (2 * se_link)),
                              right_lwr = ilink(fit_link - (2 * se_link)))
# Check the confidence interval. It should not extend below 0 (negative counts
# are impossible)
min(gam_Dino_OL_newdata$right_lwr) # Nice :)
min(gam_Dino_OL_newdata$fit_resp)
max(gam_Dino_OL_newdata$right_upr) # Nice too
max(gam_Dino_OL_newdata$fit_resp)
# But we can see that the confidence interval is ridiculously small around the
# min and max of the model fit.

# Plot

ggplot(gam_Dino_OL_newdata, aes(x = Day, y = fit_resp))+
  geom_ribbon(aes(x=Day, ymin=right_lwr, ymax=right_upr), fill = 'grey70', alpha=0.7) +
  geom_line(linewidth = 1) +
  geom_point(data = Season_Dino, aes(x = Day, y = true_count, color = Year), shape = 21) +
  scale_color_viridis_c('Year') +
  theme_classic()

# Saving plot
# ggsave('gam_Dino_OL_all.tiff', dpi = 300, height = 120, width = 160, 
#        units = 'mm', compression = 'lzw')

### Checking the model
ModelOutputs<-data.frame(Fitted=fitted(gam_Dino_OL),
                         Residuals=resid(gam_Dino_OL))

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- gam_Dino_OL$model
# then we add the values of fitted and residuals 
# (but are they in the same order as the model? -> need to check that)
qq_data <- bind_cols(qq_data, ModelOutputs)

# Plot : verify that true data matches (more or less) model fit
ggplot(qq_data)+
  geom_point(aes(x = Day, y = true_count), color = 'red') +
  geom_point(aes(x = Day, y = Fitted), color = 'blue') +
  theme_classic() +
  labs(y = "Dinophysis count", x = "Calendar day")

# It matches! great!

# And (qq-)plot
qqplot_custom <- ggplot(qq_data) +
  stat_qq(aes(sample=Residuals), color = '#2156A1', alpha = .7) +
  stat_qq_line(aes(sample=Residuals), color = '#2156A1') +
  theme_classic() +
  labs(y="Sample Quantiles",x="Theoretical Quantiles")

qqplot_custom

# Save the plot
# ggsave('qqplot_custom_16sites.tiff', dpi = 300, height = 175, width = 250,
#                units = 'mm', compression = 'lzw')

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data)+
  geom_point(aes(x=Fitted,y=Residuals), 
             color = '#2156A1', alpha = .7) +
  theme_classic() +
  labs(y="Residuals",x="Fitted Values")

RvFplot_custom

# Save the plot
# ggsave('RvFplot_custom_16sites.tiff', dpi = 300, height = 175, width = 250,
#                units = 'mm', compression = 'lzw')

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data, aes(x = Residuals))+
  geom_histogram(binwidth = 1, fill = '#2156A1')+
  theme_classic() +
  labs(x='Residuals', y = 'Count')

HistRes_custom

# Save the plot
# ggsave('HistRes_custom_16sites.tiff', dpi = 300, height = 175, width = 250,
#                units = 'mm', compression = 'lzw')

# Plotting the whole model

# Adding the date variable
# To do this, we will exploit the fact that adding a number to an object of type
# 'Date' just adds the appropriate number of days.
gam_Dino_OL_newdata <- gam_Dino_OL_newdata %>%
  # Create a month and day character variable, that is always '01-01'
  mutate(MonthDay = '01-01') %>%
  # Then unite it with Year to create something that will resemble a ymd date
  mutate(CharYear = as.character(Year)) %>%
  unite(Date, CharYear, MonthDay, sep = '-') %>%
  # And make it a 'Date' object
  mutate(Date = ymd(Date)) %>%
  # Now, add Day-1 to each
  mutate(Date = Date + (as.numeric(Day)-1))

# Ouest Loscolo
ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Dino_OL_newdata, aes(x = Date, ymin = right_lwr, 
                                     ymax = right_upr),
              color = '#7F96B6',
              fill = '#BBD4F2',
              linewidth = .75, alpha = .8) +
  geom_line(data = gam_Dino_OL_newdata, aes(x = Date, y = fit_resp),
            linewidth = 1,
            color = '#435E7B') +
  # plot the data
  geom_point(data = trueData_OL, aes(x = Date, y = true_count),
             color = '#11203E', size = 4,
             alpha = .3) +
  # x axis as date and adding limits
  scale_x_date() +
  # cut the y scale at 22
  scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(x = 'Date', y = 'Dinophysis cells observed in 10 mL') +
  theme_classic()

#### Plotting additional components of the GAM ####

### 1) Plotting smooths of linear predictors ####

# Define function "EvaluateSmooths"
# (Function from this post by unique2 on StackOverflow. Many thanks to them!)
# https://stackoverflow.com/questions/19735149/is-it-possible-to-plot-the-smooth-components-of-a-gam-fit-with-ggplot2)

EvaluateSmooths = function(model, select=NULL, x=NULL, n=100) {
  if (is.null(select)) {
    select = 1:length(model$smooth)
  }
  do.call(rbind, lapply(select, function(i) {
    smooth = model$smooth[[i]]
    data = model$model
    
    if (is.null(x)) {
      min = min(data[smooth$term])
      max = max(data[smooth$term])
      x = seq(min, max, length=n)
    }
    if (smooth$by == "NA") {
      by.level = "NA"
    } else {
      by.level = smooth$by.level
    }
    range = data.frame(x=x, by=by.level)
    names(range) = c(smooth$term, smooth$by)
    
    mat = PredictMat(smooth, range)
    par = smooth$first.para:smooth$last.para
    
    y = mat %*% model$coefficients[par]
    
    se = sqrt(rowSums(
      (mat %*% model$Vp[par, par, drop = FALSE]) * mat
    ))
    
    return(data.frame(
      label=smooth$label
      , x.var=smooth$term
      , x.val=x
      , by.var=smooth$by
      , by.val=by.level
      , value = y
      , se = se
    ))
  }))
}

# Call the function on the desired GAM
gamOL_smooths = EvaluateSmooths(gam_Dino_OL)
# In case we want only the smooth for the predictor 'Day'
gamOL_smooths_Day = filter(gamOL_smooths, x.var == 'Day')

# Plot the smooths for variables Day and Year
ggplot(gamOL_smooths, aes(x.val, value)) + 
  geom_line(color = '#2156A1') + 
  geom_line(aes(y=value + 2*se), linetype="dashed", color = '#2156A1') + 
  geom_line(aes(y=value - 2*se), linetype="dashed", color = '#2156A1') + 
  facet_wrap(facets = 'x.var', scales = 'free') +
  labs(x = 'Value of linear predictor', y = 'Smooth value', title = 'GAM smooths') +
  theme_classic()

### 2) Plotting basis functions of the GAM ####

# We use the 'basis' function from package gratia
basisOL <- basis(gam_Dino_OL)
# Let's say we only want basis functions for the 'Day' variable
basisOL_day <- filter(basisOL, .smooth == 's(Day)')

# Define a nice color palette (from the Bretagne color palettes!)
palette_splines8 <- c('#0A1635', '#76A7E2', '#FBA823', '#8D6456',
                      '#8B064B', '#FF6448', '#649003', '#1F3700')

# Plot the basis functions for variable Day
ggplot(basisOL_day, aes(x = Day, y = ilink(.value), color = .bf)) + 
  geom_line(linewidth = .75, alpha = .8) +
  scale_color_discrete(type = palette_splines8, guide = 'none') +
  labs(x = 'Calendar day', y = 'Basis spline value', 
       title = 'The 8 basis functions of the GAM') +
  theme_classic()

# Pretty nice
# Let's do a plot in 3 steps

# First step, only data points (cut at the top for clarity)
ggplot(basisOL_day, aes(x = Day, y = ilink(value), color = bf)) +
  # plot the data
  geom_point(data = plotDino_OL, aes(x=Day, y=true_count), 
             color = '#2156A1', size = 3,
             alpha = .3) +
  # cut the y scale at 20
  scale_y_continuous(limits = c(0,20)) +
  # Text
  labs(x = 'Calendar day', y = 'Dinophysis cells observed in 10 mL') +
  theme_classic()

# Save the plot (step 1)
# ggsave('educGAM_step1.tiff', height = 150, width = 300,
# dpi = 300, unit = 'mm', compression = 'lzw')

# Second step, plotting the basis functions of the GAM
ggplot(basisOL_day, aes(x = Day, y = ilink(.value), color = .bf)) +
  # plot the data
  geom_point(data = plotDino_OL, aes(x=Day, y=true_count), 
             color = '#2156A1', size = 3,
             alpha = .3) +
  # plot the basis functions
  geom_line(linewidth = 1, alpha = .8) +
  # color scale
  scale_color_discrete(type = palette_splines8, guide = 'none') +
  # cut the y scale at 20
  scale_y_continuous(limits = c(0,20)) +
  # Text
  labs(x = 'Calendar day', y = 'Dinophysis cells observed in 10 mL') +
  theme_classic()

# Save the plot (step 2)
# ggsave('educGAM_step2.tiff', height = 150, width = 300,
# dpi = 300, unit = 'mm', compression = 'lzw')

# Third step, with the smooth
# /!\ This step uses objects defined in the previous section of the script /!\
ggplot(basisOL_day, aes(x = Day, y = ilink(.value), color = .bf)) +
  # plot the data
  geom_point(data = plotDino_OL, aes(x=Day, y=true_count), 
             color = '#2156A1', size = 3,
             alpha = .3) +
  # plot the basis functions
  geom_line(linewidth = 1, alpha = .8) +
  # plot the smooth
  geom_line(data = gamOL_smooths_Day, aes(x.val, ilink(value)),
            color = '#2156A1', linewidth = 1.5) +
  # And the std error
  geom_line(data = gamOL_smooths_Day, aes(x = x.val, y=ilink(value + 2*se)), 
            linetype="dashed", color = '#2156A1', alpha = .8, linewidth = 1) + 
  geom_line(data = gamOL_smooths_Day, aes(x = x.val, y=ilink(value - 2*se)), 
            linetype="dashed", color = '#2156A1', alpha = .8, linewidth = 1) + 
  # color scale
  scale_color_discrete(type = palette_splines8, guide = 'none') +
  # cut the y scale at 20
  scale_y_continuous(limits = c(0,20)) +
  # Text
  labs(x = 'Calendar day', y = 'Dinophysis cells observed in 10 mL') +
  theme_classic()

# Save the plot (step 3)
# ggsave('educGAM_step3_frontslide.tiff', height = 300, width = 900,
# dpi = 300, unit = 'mm', compression = 'lzw')

# Fourth step, only the final GAM smooth
# /!\ This step uses objects defined in the previous section of the script /!\
ggplot() +
  # plot the data
  geom_point(data = plotDino_OL, aes(x=Day, y=true_count), 
             color = '#2156A1', size = 3,
             alpha = .3) +
  # plot the smooth
  geom_line(data = gamOL_smooths_Day, aes(x.val, ilink(value)),
            color = '#2156A1', linewidth = 1.5) +
  # And the std error
  geom_line(data = gamOL_smooths_Day, aes(x = x.val, y=ilink(value + 2*se)), 
            linetype="dashed", color = '#2156A1', alpha = .8, linewidth = 1) + 
  geom_line(data = gamOL_smooths_Day, aes(x = x.val, y=ilink(value - 2*se)), 
            linetype="dashed", color = '#2156A1', alpha = .8, linewidth = 1) + 
  # cut the y scale at 20
  scale_y_continuous(limits = c(0,20)) +
  # Text
  labs(x = 'Calendar day', y = 'Dinophysis cells observed in 10 mL') +
  theme_classic()

ggplot() +
  # plot the data
  geom_point(data = plotDino_OL, aes(x=Day, y=true_count), 
             color = '#2156A1', size = 3,
             alpha = .3) +
  # plot the smooth
  geom_line(data = gamOL_smooths_Day, aes(x.val, ilink(value)),
            color = '#2156A1', linewidth = 1.5) +
  # And the std error
  geom_line(data = gamOL_smooths_Day, aes(x = x.val, y=ilink(value + 2*se)), 
            linetype="dashed", color = '#2156A1', alpha = .8, linewidth = 1) + 
  geom_line(data = gamOL_smooths_Day, aes(x = x.val, y=ilink(value - 2*se)), 
            linetype="dashed", color = '#2156A1', alpha = .8, linewidth = 1) +
  # cut the y scale at 20
  scale_y_continuous(limits = c(0,20)) +
  # Text
  labs(x = 'Calendar day', y = 'Dinophysis cells observed in 10 mL') +
  theme_classic()

# Save the plot (step 4)
# ggsave('educGAM_step4.tiff', height = 150, width = 300,
# dpi = 300, unit = 'mm', compression = 'lzw')

#### Plotting the GAM on the different years ####
## The logic here is to see how well the GAM fits the data year by year
# Can we see some years where the phenology significantly deviates from
# what the GAM predicts?
## Ouest Loscolo ####
response_pred_plot <- read.csv2('response_pred_plot_20240813_14sites.csv', 
                                header = TRUE, fileEncoding = 'ISO-8859-1')

rpp_OL <- filter(response_pred_plot, Code_point_Libelle == 'Ouest Loscolo')
Season_Dino_OL <- filter(Season_Dino, Code_point_Libelle == 'Ouest Loscolo')

ggplot(Season_Dino_OL, aes(x = Day, y = true_count)) +
  geom_point(color = 'dodgerblue3', alpha = .8) +
  geom_line(data = rpp_OL, aes(x = Day, y = median.fit), color = 'dodgerblue3') +
  geom_ribbon(data = rpp_OL, aes(x = Day, y = median.fit,
                                 ymin = lwrS, ymax = uprS), 
              color = 'dodgerblue3', alpha = 0.2) +
  facet_wrap(facets = c('Year'), scales = 'free_y') +
  # scale_y_continuous(limits = c(0,30)) +
  theme_classic()
## Cabourg ####

rpp_Cab <- filter(response_pred_plot, Code_point_Libelle == 'Cabourg')
Season_Dino_Cab <- filter(Season_Dino, Code_point_Libelle == 'Cabourg')

ggplot(Season_Dino_Cab, aes(x = Day, y = true_count)) +
  geom_point(color = 'orangered', alpha = .8) +
  geom_line(data = rpp_Cab, aes(x = Day, y = median.fit), color = 'orangered') +
  geom_ribbon(data = rpp_Cab, aes(x = Day, y = median.fit,
                                 ymin = lwrS, ymax = uprS), 
              color = 'orangered', alpha = 0.2) +
  facet_wrap(facets = c('Year'), scales = 'free_y') +
  # scale_y_continuous(limits = c(0,30)) +
  theme_classic()
