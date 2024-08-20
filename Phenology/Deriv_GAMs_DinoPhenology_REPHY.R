#### Derivatives of fitted GAMs for Dinophysis phenology ##
### V. POCHIC
# 2024-08-20

# /!\ This script requires data tables generated with the 'GAM Dino unified more sites' 
# and 'GAM Meso unified' scripts /!\

# The idea behind this is to identify periods of Dinophysis accumulation (deriv-
# ative > 0) and loss (derivative < 0 ) in the different Dinophysis phenologies.

#### Packages and functions ####
library(tidyverse)
library(mgcv)
library(gratia)
library(ggnewscale)
library(cmocean)
library(RColorBrewer)

### deriv function
source('derivFun.R')
#### We'll need the GAM models (obviously) ####

### Import the 'seasonalized' data ####

Season_Dino <- read.csv2('Season_Dino.csv', header = TRUE, 
                         fileEncoding = 'ISO-8859-1')

# The sites must be entered as factors
Season_Dino <- Season_Dino %>%
  mutate(Code_point_Libelle = as.factor(Code_point_Libelle))

# Alright!

### Mesodinium
# The table of Mesodinium counts
Table_Meso_zeros <- read.csv2('Table_Meso_zeros.csv', header = TRUE, 
                              fileEncoding = 'ISO-8859-1',
                              # ensures that special characters (like +) are not
                              # changed during import
                              check.names = FALSE)

# Seasonality
Season_Meso <- Table_Meso_zeros %>%
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
  mutate(across(.cols = c('Mesodinium', 'Mesodinium rubrum'), .fns = ~ ifelse(. %% 100 != 0, 0, .))) %>%
  # And creating some count variables for Dinophysis as a genus
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
  filter(true_count < 500) %>%
  # converting the site to factor for the model
  mutate(Code_point_Libelle = as.factor(Code_point_Libelle))

# Great

#### Fitting GAMs ####

### Dinophysis
gam_Dino <- gam(data = Season_Dino, 
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

# We can have a look at the diagnostic plots to be sure it worked as intended
# Model outputs
ModelOutputs<-data.frame(Fitted=fitted(gam_Dino),
                         Residuals=resid(gam_Dino))

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- gam_Dino$model
# then we add the values of fitted and residuals 
# (they are in the same order as in the model)
qq_data <- bind_cols(qq_data, ModelOutputs)

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

# color palette for 16 sites
pheno_palette16 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646',
                     '#642C3A', '#DEB1CC', '#FC4D6B', '#791D40')

# And (qq-)plot
qqplot_custom <- ggplot(qq_data_reordered) +
  stat_qq(aes(sample=Residuals, color = Code_point_Libelle), alpha = .7) +
  stat_qq_line(aes(sample=Residuals, color = Code_point_Libelle)) +
  facet_wrap(facets = c('Code_point_Libelle')) +
  scale_color_discrete(type = pheno_palette16, guide = 'none') +
  theme_classic() +
  labs(y="Sample Quantiles",x="Theoretical Quantiles")

qqplot_custom

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data_reordered)+
  geom_point(aes(x=Fitted,y=Residuals, color  =Code_point_Libelle), 
             alpha = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free') +
  scale_color_discrete(type = pheno_palette16, guide = 'none') +
  theme_classic() +
  labs(y="Residuals",x="Fitted Values")

RvFplot_custom

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data_reordered, aes(x = Residuals, 
                                                fill = Code_point_Libelle))+
  geom_histogram(binwidth = 1)+
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free') +
  scale_fill_discrete(type = pheno_palette16, guide = 'none') +
  theme_classic() +
  labs(x='Residuals', y = 'Count')

HistRes_custom

# Remove all the diagnostic plots
rm(HistRes_custom)
rm(RvFplot_custom)
rm(qqplot_custom)

### Mesodinium
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

### QQ plots of both GAMs are quite shitty.


#### Computing first derivative of GAMs ####

## We'll be using the 'derivatives()' function from the 'gratia' package

# First, we need new data for the model to compute the derivative on
gam_Dino_newdata <- expand_grid(Day=seq(1, 365),
                                # We add a Year vector as it has become a factor
                                # of the model
                                Year=seq(min(Season_Dino$Year), 
                                         max(Season_Dino$Year)),
                                # We also add the site data
                                Code_point_Libelle = unique(Season_Dino$Code_point_Libelle))


### Dinophysis
gam_Dino.d <- derivatives(# The object of which we want to calculate derivatives
  gam_Dino,
  # On which smooth term do we want to calculate the derivative?
  # term = 'Day',
  # # partial match for the term
  # partial_match = TRUE,
  # New data to be called upon
  data = gam_Dino_newdata,
  # Order of the derivative (here, we want first order)
  order = 1,
  # The type of finite difference used (?)
  type = 'central',
  # n is the number of points to evaluate the derivative at
  n = 365,
  # The type of confidence interval to compute. We want simultaneous
  interval = 'simultaneous',
  # the number of simulations used in computing the simultaneous intervals
  n_sim = 10000, # let's not start too big...
  # level of confidence interval (0.95 for 95% CI is default)
  level = 0.95,
  # Use the bayesian covariance matrix?
  unconditional = TRUE, # let's try...
  # Number of cores (to try and do this in less than 56 hours)
  ncores = 4
  )

# It works!
# Let's try to plot that shall we

### Plotting

# Aesthetics
# Nice color palette with 16 colors
pheno_palette16 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646',
                     '#642C3A', '#DEB1CC', '#FC4D6B', '#791D40')

# Relevel factors so we have sampling sites in the desired order
gam_Dino.d <- gam_Dino.d %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

# A nice plot
ggplot(gam_Dino.d, aes(x = data, y = derivative, 
                           color = Code_point_Libelle,
                           fill = Code_point_Libelle)) +
  # Confidence interval
         geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  # Derivative fit
         geom_path(lwd = 1) +
  # Draw a line at 0 to separate accumulation from loss
         geom_line(aes(x = data, y = 0), color = 'grey10', linewidth = .7) +
  labs(y = "1st derivative of Dinophysis GAM",
       x = "Day of the year",
       title = "1st derivative of Dinophysis GAM"
  ) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette16, guide = 'none') +
  scale_fill_discrete(type = pheno_palette16, guide = 'none') +
  theme_classic()

# That seems to work nicely

# Saving plot
# ggsave('DinoDeriv_16sites_phenology.tiff', dpi = 300, height = 225, width = 300,
#        units = 'mm', compression = 'lzw')

### We can also do that with only the sites we will analyse further
gam_Dino.d_select <- gam_Dino.d %>%
  # filter
  filter(Code_point_Libelle %in% c('Point 1 Boulogne', 'At so',
                                   'Antifer ponton pétrolier', 'Cabourg',
                                   'les Hébihens', 'Loguivy',
                                   'Men er Roue', 'Ouest Loscolo',
                                   'Le Cornard', 'Auger',
                                   'Arcachon - Bouée 7', 'Teychan bis')) %>%
  # and relevel
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis'))

# New color palette needed
pheno_palette12 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646')

# Select plot
ggplot(gam_Dino.d_select, aes(x = data, y = derivative, 
                       color = Code_point_Libelle,
                       fill = Code_point_Libelle)) +
  # Confidence interval
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  # Derivative fit
  geom_path(lwd = 1) +
  # Draw a line at 0 to separate accumulation from loss
  geom_line(aes(x = data, y = 0), color = 'grey10', linewidth = .7) +
  labs(y = "1st derivative of Dinophysis GAM",
       x = "Day of the year",
       title = "1st derivative of Dinophysis GAM"
  ) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y',
             # 4 rows to highlight the latitudinal change in phenology
             nrow = 4) +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  scale_fill_discrete(type = pheno_palette12, guide = 'none') +
  theme_classic()

#### Annotating the derivatives ####
### We can assign ecological meaning to the derivatives
# And the confidence interval allows us to add a degree of certainty.

Deriv_Dino <- gam_Dino.d %>%
  mutate(deriv_sign = ifelse(
    # If deriv > 0 and lower > 0, this is a 'true positive' (coded as 2)
    (derivative > 0 & lower > 0), 2,
    # If deriv > 0 and lower <= 0, this is a 'positive' (coded as 1)
      ifelse(derivative > 0 & lower <= 0, 1,
    # If derivative = 0, it's a 'null' (coded as 0)
        ifelse(derivative == 0, 0,
    # If derivative < 0 and upper >= 0, this is a 'negative' (coded as -1)
          ifelse(derivative < 0 & upper >= 0, -1,
    # If derivative < 0 and upper < 0, this is a 'true negative' (coded as -2)
            -2)
            )
          )
        )
  ) %>%
  mutate(Day = data) %>%
  select(c('Code_point_Libelle', 'Day', 'derivative', 'upper', 'lower',
           'deriv_sign')) %>%
  group_by(Code_point_Libelle, Day) %>%
  summarise(across(c('derivative', 'upper', 'lower',
                     'deriv_sign'),
                   mean),
                   .groups = 'keep')

# Great! We can work with that 'Deriv_Dino' file for plots
# Let's save it
# write.csv2(Deriv_Dino, 'Deriv_Dino_20240424.csv', row.names = FALSE,
#            fileEncoding = 'ISO-8859-1')

#### Plotting the derivative of Dinophysis along the Mesodinium model ####

### Import Mesodinium GAM data for plotting ####
pred_plot <- read.csv2('pred_plot_MESO_20240418.csv', header = TRUE,
                            fileEncoding = 'ISO-8859-1')

# Compute the response pred plot
# We need the inverse link function of gam_Meso
ilink <- gam_Meso$family$linkinv

# Transform values so they are expressed in the response variable units
# (by applying the inverse link function)
response_pred_plot <- as.data.frame(pred_plot) %>%
  mutate(across(c('median.fit', 'lwrS', 'uprS', 'lwrP', 'uprP'), ~ ilink(.)))

# Let's save the response_pred_plot for the Meso GAM
# write.csv2(response_pred_plot, 'response_pred_plot_MESO_20240424.csv',
#            row.names = FALSE, fileEncoding = 'ISO-8859-1')

# We have our response pred plot that we can plot to check if everything is ok:

#### Plotting the vanilla Meso GAM ####
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

## Cropping out data points for plotting ####
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


## Plot this #####
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

#### Applying a rugplot of Dinophysis derivative under the Mesodinium GAM ####

### Joining datasets
combined_plot <- right_join(response_pred_plot,
                           Deriv_Dino,
                           by = c('Code_point_Libelle', 'Day'),
                           suffix = c('','')) %>%
  # relevel factors
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))
  

### A different color palette for the derivative rugplot (discrete option)
deriv_palette <- c(# -2, -1
  'indianred3', 'indianred1',
  # +1, +2
  'lightblue2', 'lightslateblue')

# Plot
ggplot(combined_plot) +
  # First part of plot (a rug plot)
  # We use a color palette from the cmocean package
  scale_color_cmocean(name = 'curl', direction = -1) +
  geom_rug(aes(x = Day,
               color = derivative),
           linewidth = .3
  ) +
  # Labels
  labs(y = "Mesodinium cells observed in 10 mL",
       x = "Day of the year",
       title = "Mesodinium and Dinophysis dynamics",
       color = c(expression(paste('Dinophysis loss/accumulation (d'^'-1',')')))
       #subtitle = sprintf("Each line is one of %i draws from the Bayesian posterior distribution of the model", nrnd)
  ) +
  # Change the color scale
  new_scale_color() +
  # Second part of plot
  geom_ribbon(aes(x = Day,
                  ymin = lwrS, ymax = uprS,
                  color = Code_point_Libelle,
                  fill = Code_point_Libelle), alpha = 0.1) +
  # geom_path(data = stackFits_response, mapping = aes(y = values_response, x = Day, group= ind),
  #           alpha = 0.1, colour = "grey20") +
  geom_path(aes(x = Day, y = median.fit,
                color = Code_point_Libelle),
            lwd = 1) +
  geom_point(data = Season_Meso_crop, aes(x = Day, y = true_count,
                                          color = Code_point_Libelle), 
             size = .8, alpha = .5) +
  # Facet
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  scale_fill_discrete(type = pheno_palette12, guide = 'none') +
  theme(
    # panel
    panel.background = element_rect(fill = 'transparent', 
                                    linewidth = .3, 
                                    color = 'grey10'),
    # legend
    legend.background = element_rect(linewidth = .5, color = 'grey10'),
    legend.title = element_text(size = 11, color = 'grey5'),
    legend.frame = element_rect(linewidth = .5, color = 'grey10'),
    legend.ticks = element_line(linewidth = .2, color = 'grey25'),
    legend.position = 'bottom',
    # grid
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    # Facet labels
    strip.background = element_rect(fill = 'grey80',
                                    linewidth = .2,
                                    color = 'grey10'),
    strip.text = element_text(color = 'grey5')
  )

# Saving plot
# ggsave('Meso_and_DinoDeriv_12sites.tiff', dpi = 300, height = 250, width = 300,
#        units = 'mm', compression = 'lzw')



#### Only for Seine Bay and South Brittany ####

# Only the 4 sampling sites that we want
select_plot <- combined_plot %>%
  filter(Code_point_Libelle %in% c('Antifer ponton pétrolier', 'Cabourg',
                                   'Men er Roue', 'Ouest Loscolo'))

Season_Meso_select <- Season_Meso_crop %>%
  filter(Code_point_Libelle %in% c('Antifer ponton pétrolier', 'Cabourg',
                                   'Men er Roue', 'Ouest Loscolo'))

## Plot (continuous version) ####
ggplot(select_plot) +
  # First part of plot (a rug plot)
  # We use a color palette from the cmocean package
  scale_color_cmocean(name = 'curl', direction = -1) +
  geom_rug(aes(x = Day,
               color = derivative),
           linewidth = 1
  ) +
  # Labels
  labs(y = "Mesodinium cells observed in 10 mL",
       x = "Day of the year",
       title = "Mesodinium and Dinophysis dynamics",
       color = c(expression(paste('Dinophysis loss/accumulation (d'^'-1',')')))
       #subtitle = sprintf("Each line is one of %i draws from the Bayesian posterior distribution of the model", nrnd)
  ) +
  # Change the color scale
  new_scale_color() +
  # Second part of plot
  geom_ribbon(aes(x = Day,
                  ymin = lwrS, ymax = uprS,
                  color = Code_point_Libelle,
                  fill = Code_point_Libelle), alpha = 0.2) +
  # geom_path(data = stackFits_response, mapping = aes(y = values_response, x = Day, group= ind),
  #           alpha = 0.1, colour = "grey20") +
  geom_path(aes(x = Day, y = median.fit,
                color = Code_point_Libelle),
            lwd = 1) +
  geom_point(data = Season_Meso_select, aes(x = Day, y = true_count,
                                          color = Code_point_Libelle), 
             size = .8, alpha = .5) +
  # Facet
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  scale_fill_discrete(type = pheno_palette12, guide = 'none') +
  theme(
    # panel
    panel.background = element_rect(fill = 'transparent', 
                                    linewidth = .3, 
                                    color = 'grey10'),
    # legend
    legend.background = element_rect(linewidth = .5, color = 'grey10'),
    legend.title = element_text(size = 11, color = 'grey5'),
    legend.frame = element_rect(linewidth = .5, color = 'grey10'),
    legend.ticks = element_line(linewidth = .2, color = 'grey25'),
    legend.position = 'bottom',
    # grid
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    # Facet labels
    strip.background = element_rect(fill = 'grey80',
                                    linewidth = .2,
                                    color = 'grey10'),
    strip.text = element_text(color = 'grey5')
  )

# Saving plot
# ggsave('Meso_and_DinoDeriv_4sites.tiff', dpi = 300, height = 170, width = 160,
#        units = 'mm', compression = 'lzw')

## Plot (discrete version) ####

# Discrete color palette
deriv_palette <- c(# -2, -1
  'indianred3', 'indianred1',
  # +1, +2
  'lightblue2', 'lightslateblue')

# Plot
ggplot(select_plot) +
  # First part of plot (a rug plot)
  # We use a color palette from the cmocean package
  scale_color_discrete(type = deriv_palette) +
  geom_rug(aes(x = Day,
               color = factor(deriv_sign)),
           linewidth = 1
  ) +
  # Labels
  labs(y = "Mesodinium cells observed in 10 mL",
       x = "Day of the year",
       title = "Mesodinium and Dinophysis dynamics",
       color = c(expression(paste('Dinophysis loss/accumulation')))
       #subtitle = sprintf("Each line is one of %i draws from the Bayesian posterior distribution of the model", nrnd)
  ) +
  # Change the color scale
  new_scale_color() +
  # Second part of plot
  geom_ribbon(aes(x = Day,
                  ymin = lwrS, ymax = uprS,
                  color = Code_point_Libelle,
                  fill = Code_point_Libelle), alpha = 0.2) +
  # geom_path(data = stackFits_response, mapping = aes(y = values_response, x = Day, group= ind),
  #           alpha = 0.1, colour = "grey20") +
  geom_path(aes(x = Day, y = median.fit,
                color = Code_point_Libelle),
            lwd = 1) +
  geom_point(data = Season_Meso_select, aes(x = Day, y = true_count,
                                            color = Code_point_Libelle), 
             size = .8, alpha = .5) +
  # Facet
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette12, guide = 'none') +
  scale_fill_discrete(type = pheno_palette12, guide = 'none') +
  theme(
    # panel
    panel.background = element_rect(fill = 'transparent', 
                                    linewidth = .3, 
                                    color = 'grey10'),
    # legend
    legend.background = element_rect(linewidth = .5, color = 'grey10'),
    legend.title = element_text(size = 11, color = 'grey5'),
    legend.frame = element_rect(linewidth = .5, color = 'grey10'),
    legend.ticks = element_line(linewidth = .2, color = 'grey25'),
    legend.position = 'bottom',
    # grid
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    # Facet labels
    strip.background = element_rect(fill = 'grey80',
                                    linewidth = .2,
                                    color = 'grey10'),
    strip.text = element_text(color = 'grey5')
  )

# Saving plot
# ggsave('Meso_and_DinoDeriv_4sites_discrete.tiff', dpi = 300, height = 170, width = 160,
#        units = 'mm', compression = 'lzw')

#### Rugplot of environmental variables under the Dino derivative plot ####
### Import data for environmental variables ####
Table_hydro_fortnightly <- read.csv2('Table_hydro_fortnightly_20240820.csv', 
                                     header = TRUE,
                                     fileEncoding = 'ISO-8859-1')

# Create a daily table for plotting the rugplot
# We create a vector for 'Day of the year'
Daily_basis <- expand_grid(Day = seq(1,365))
# And one for 'Fortnight' that matches
Daily_basis$Fortnight = as.vector(c(rep(1:26, each = 14), 26))
# We add the site data
Daily_basis <- expand_grid(Daily_basis,
                           Code_point_Libelle = unique(Table_hydro_fortnightly$Code_point_Libelle))

# We arrange Daily_basis by Site
Daily_basis <- Daily_basis %>%
  group_by(Code_point_Libelle, Day) %>%
  arrange(Code_point_Libelle)

# Nickel chrome!

Table_hydro_daily <- left_join(Daily_basis, Table_hydro_fortnightly,
                               by = c('Code_point_Libelle', 'Fortnight'),
                               suffix  = c('','')) %>%
  filter(is.na(TEMP.med) == FALSE) %>%
  mutate(Code_point_Libelle = as.factor(Code_point_Libelle)) %>%
  group_by(Day, Code_point_Libelle)

# Let's save this hydrology table for later
# write.csv2(Table_hydro_daily, 'Table_hydro_daily_20240820.csv', row.names = FALSE,
#            fileEncoding = 'ISO-8859-1')

### Plotting ####

# We want to plot the Dino derivative with the temperature as rugplot
# A nice plot
ggplot(gam_Dino.d_select) +
  # First part of the plot: the rug plot
  # We use a color palette from the RColorBrewer package
  scale_color_distiller(palette = 'RdBu', direction = -1) +
  geom_rug(data = Table_hydro_daily, aes(x = Day, color = TEMP.med),
           linewidth = .3,
           length = unit(0.5, 'cm')
  ) +
  # Labels
  labs(y = "1st derivative of Dinophysis GAM",
       x = "Day of the year",
       title = '1st derivative of Dinophysis GAM',
       color = c(expression(paste('Sea surface temperature (°C)')))
  ) +
  # Change the color scale
  new_scale_color() +
  # Second part of plot: GAM derivatives
  # Confidence interval
  geom_ribbon(aes(x = data, ymin = lower, ymax = upper,
                  color = Code_point_Libelle,
                  fill = Code_point_Libelle), alpha = 0.2) +
  # Derivative fit
  geom_path(aes(x = data, y = derivative, 
                  color = Code_point_Libelle), lwd = 1) +
  # Draw a line at 0 to separate accumulation from loss
  geom_line(aes(x = data, y = 0), color = 'grey10', linewidth = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette16, guide = 'none') +
  scale_fill_discrete(type = pheno_palette16, guide = 'none') +
  # Theme
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

# Saving the temperature plot
# ggsave('Dinoderiv_temperature_12sites_large.tiff', dpi = 300, height = 200, width = 250,
#        units = 'mm', compression = 'lzw')

# Now with salinity
ggplot(gam_Dino.d_select) +
  # First part of the plot: the rug plot
  # We use a color palette from the cmocean package
  scale_color_cmocean(name = 'haline') +
  geom_rug(data = Table_hydro_daily, aes(x = Day, color = SALI.med),
           linewidth = .3,
           length = unit(0.5, 'cm')
  ) +
  # Labels
  labs(y = "1st derivative of Dinophysis GAM",
       x = "Day of the year",
       title = '1st derivative of Dinophysis GAM',
       color = c(expression(paste('Sea surface salinity (PSU)')))
  ) +
  # Change the color scale
  new_scale_color() +
  # Second part of plot: GAM derivatives
  # Confidence interval
  geom_ribbon(aes(x = data, ymin = lower, ymax = upper,
                  color = Code_point_Libelle,
                  fill = Code_point_Libelle), alpha = 0.2) +
  # Derivative fit
  geom_path(aes(x = data, y = derivative, 
                color = Code_point_Libelle), lwd = 1) +
  # Draw a line at 0 to separate accumulation from loss
  geom_line(aes(x = data, y = 0), color = 'grey10', linewidth = .7) +
  facet_wrap(facets = c('Code_point_Libelle'), scales = 'free_y') +
  # Set the color palette :
  scale_color_discrete(type = pheno_palette16, guide = 'none') +
  scale_fill_discrete(type = pheno_palette16, guide = 'none') +
  # Theme
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

# Saving the temperature plot
# ggsave('Dinoderiv_salinity_12sites_large.tiff', dpi = 300, height = 200, width = 250,
#        units = 'mm', compression = 'lzw')
