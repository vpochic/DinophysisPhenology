###### GAM for Dinophysis phenology with mgcv ###
## V. POCHIC
# 2025-06-04

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

# Select only stations of interest from 2006 --> 2022
Table_phyto_sites <- filter(Table_phyto_taxon, Code_point_Libelle %in% 
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

########## Organizing the data #######

#### More zeros ####

### The goal here is to create a table where all absences of detection are 
# listed as zeros
Table_Dino_zeros <- Table_phyto_sites %>%
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
rm(Table_phyto_sites)

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
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre')) %>%
  # Expanding on date variable (creating a fortnight variable)
  mutate(Week = week(Date)) %>%
  # 'ceiling' takes the upper integer of a decimal number
  mutate(Fortnight = ceiling(Week/2))

# Checking the homogeneity of the sampling between years
# We make a histogram of samples for the 16 years (2007-2022)

# Color palette
pheno_palette16 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646',
                     '#642C3A', '#DEB1CC', '#FC4D6B', '#791D40')

ggplot(data = Season_Dino) +
  geom_histogram(aes(x = as_factor(Year), fill = Code_point_Libelle),
                 stat = 'count') +
  facet_wrap(facets = 'Code_point_Libelle') +
  scale_fill_discrete(type = pheno_palette16, guide = 'none') +
  scale_x_discrete(breaks = c('2007', '2010', '2015', '2020', '2022')) +
  theme_classic()
# Some heterogeneity in some sites (Antifer, notably). Atlantic sites 
# experienced a drop in 2020 due to covid19 pandemic.

# Homogeneity across the year, for the 26 fortnights in a year 
# (fortnightly sampling)
ggplot(data = Season_Dino) +
  geom_histogram(aes(x = as_factor(Fortnight), fill = Code_point_Libelle),
                 stat = 'count') +
  facet_wrap(facets = 'Code_point_Libelle') +
  scale_fill_discrete(type = pheno_palette16, guide = 'none') +
  scale_x_discrete(breaks = c('1', '5', '10', '15', '20', '26')) +
  theme_classic()
# Approximately OK, homogeneity is the rule

# This seems quite fine!!!

# Let's save that
# write.csv2(Season_Dino, 'Data/REPHY_outputs/Season_Dino_20250604.csv', row.names = FALSE,
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

#### Environmental variables ####

# Here we will call a table of environmental data to plot a heatmap
# underneath our Dinophysis data
# This table was created by the script "Dino_phenology_heatmaps.R"

Table_hydro_daily <- read.csv2('Data/REPHY_outputs/Table_hydro_daily_20240821.csv', 
                                     header = TRUE,
                                     fileEncoding = 'ISO-8859-1')

# Now, we plot! First, what we want to do is to have an idea of the distribution
# of Dinophysis maxima for each site. We will focus on Channel/Atlantic sites,
# and only those where we observe Dinophysis
Table_hydro_daily_select <- filter(Table_hydro_daily,
                                   Code_point_Libelle %in%
                                     c('Antifer ponton pétrolier', 'Cabourg',
                                       'Men er Roue', 'Ouest Loscolo',
                                       'Le Cornard', 'Auger',
                                       'Arcachon - Bouée 7', 'Teychan bis'))
Maxima_Dino_stats_select <- filter(Maxima_Dino_stats,
                                   Code_point_Libelle %in%
                                     c('Antifer ponton pétrolier', 'Cabourg',
                                       'Men er Roue', 'Ouest Loscolo',
                                       'Le Cornard', 'Auger',
                                       'Arcachon - Bouée 7', 'Teychan bis'))
Maxima_Dino_select <- filter(Maxima_Dino,
                                   Code_point_Libelle %in%
                                     c('Antifer ponton pétrolier', 'Cabourg',
                                       'Men er Roue', 'Ouest Loscolo',
                                       'Le Cornard', 'Auger',
                                       'Arcachon - Bouée 7', 'Teychan bis'))

# New color scale
pheno_palette8 <- c('red3', 'orangered', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646')

# First, let's plot the distribution of the maxima (without the heatmap)
ggplot(Maxima_Dino_stats_select, aes(color = Code_point_Libelle)) +
  # All observations as small translucid points
  geom_point(data = Maxima_Dino_select, aes(x = Day, y = Code_point_Libelle, 
                                     color = Code_point_Libelle), 
             size = 3, alpha = .5) +
  scale_color_discrete(type = pheno_palette8, guide = 'none') +
  scale_x_continuous(limits = c(1,365), breaks = c(1, 100, 200, 300, 365)) +
  # reverse y axis to get sites in the desired order
  scale_y_discrete(limits = rev) +
  labs(x = 'Day of the year', y = NULL) +
  theme_classic()

# Now we try to underlie a heatmap of chl a seasonality
ggplot(Maxima_Dino_stats_select) +
  geom_tile(data = Table_hydro_daily_select, 
            aes(x = Day, y = Code_point_Libelle, fill = CHLOROA.med), 
            alpha = .8) +
  # color palette for fill
  # All observations as small translucid points
  geom_point(data = Maxima_Dino_select, aes(x = Day, y = Code_point_Libelle, 
                                            color = Code_point_Libelle), 
             size = 3, alpha = .5) +
  scale_color_discrete(type = pheno_palette8, guide = 'none') +
  scale_x_continuous(limits = c(1,365), breaks = c(1, 100, 200, 300, 365)) +
  # reverse y axis to get sites in the desired order
  scale_y_discrete(limits = rev) +
  scale_fill_cmocean(name = 'algae') +
  labs(x = 'Day of the year', y = NULL, fill = 'Median [chl a]') +
  theme(plot.title = element_text(size = 11), 
        # Axis
        axis.title.x = element_text(size=10), 
        axis.title.y =element_text(size=10), 
        axis.text = element_text(size=8, color = 'black'),
        axis.line.x = element_line(linewidth = .2, color = 'black'),
        axis.line.y = element_line(linewidth = .2, color = 'black'),
        # Legend
        legend.background = element_rect(linewidth = .5, color = 'grey10'),
        legend.title = element_text(size = 10, color = 'grey5'),
        legend.frame = element_rect(linewidth = .5, color = 'grey10'),
        legend.ticks = element_line(linewidth = .2, color = 'grey25'),
        legend.position = 'bottom',
        # Panel
        panel.grid = element_blank(),
        panel.background = element_blank(),
        # Facet labels
        strip.background = element_rect(fill = 'transparent',
                                        linewidth = 1,
                                        color = 'grey10'),
        strip.text = element_text(color = 'grey5', size = 7.5))

# Nice one. Let's save that
# ggsave('Dino_max_chloroa.tiff', dpi = 300, height = 180, width = 150, 
#        units = 'mm', compression = 'lzw')

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
# write.csv2(response_pred, 'Data/GAM_outputs/response_pred_GAMDino_20250604.csv', row.names = FALSE,
#            fileEncoding = 'ISO-8859-1')

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
# write.csv2(pred_plot, 'Data/GAM_outputs/pred_plot_20250604_16sites.csv', row.names = FALSE,
#            fileEncoding = "ISO-8859-1")

# Transform values so they are expressed in the response variable units
# (by applying the inverse link function)
response_pred_plot <- as.data.frame(pred_plot) %>%
  mutate(across(c('median.fit', 'lwrS', 'uprS', 'lwrP', 'uprP'), ~ ilink(.)))

# And save that shall we
# write.csv2(response_pred_plot, 'Data/GAM_outputs/response_pred_plot_20250604_16sites.csv',
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


###### Plot aesthetics #########

# Defining a beautiful color palette
pheno_palette16 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646',
                     '#642C3A', '#DEB1CC', '#FC4D6B', '#791D40')

# Import 'response_pred_plot' if necessary
# response_pred_plot <- read.csv2('Data/GAM_outputs/response_pred_plot_20250604_16sites.csv',
#                                 header = TRUE, fileEncoding = 'ISO-8859-1')

# And Season_Dino
# Season_Dino <- read.csv2('Data/REPHY_outputs/Season_Dino_20250604.csv',
#                          header = TRUE, fileEncoding = 'ISO-8859-1')

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
  scale_color_discrete(type = pheno_palette16, guide = 'none') +
  scale_fill_discrete(type = pheno_palette16, guide = 'none') +
  theme_classic()

# Saving plot
# ggsave('Plots/GAMs/gam_Dino_16sites_crop.tiff', dpi = 300, height = 175, width = 250,
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
  scale_color_discrete(type = pheno_palette16, guide = 'none') +
  scale_fill_discrete(type = pheno_palette16, guide = 'none') +
  theme_classic()

# Saving plot
# ggsave('Plots/GAMs/gam_Dino_allsites_crop_cells_per_L.tiff', dpi = 300, height = 175, width = 250,
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
