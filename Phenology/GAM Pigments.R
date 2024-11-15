### GAMs of pigment concentrations ###
# Part of the Dinophysis Phenology project
# V. POCHIC 2024-

### Required packages ####

library(tidyverse)
library(mgcv)
library(gratia)
library(ggpubr)
library(grid)

### Import data ####

Table_pigments <- read.csv2('Table_pigments_2007-2022.csv', header = TRUE,
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

gam_CHLOROA <- gam(data = Table_CHLOROA, 
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
                   family = gaussian(),
                   # Restricted maximum likelihood estimation (recommended method)
                   method = 'REML')

summary(gam_CHLOROA)
gam.check(gam_CHLOROA)

# gam.check indicates that there might be an issue with the random effect
# smooth for 'Year' (low p-value)

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
# Check the confidence interval. It should not extend below 0 (negative counts
# are impossible)
min(gam_CHLOROA_newdata$right_lwr) # Nice :)
min(gam_CHLOROA_newdata$fit_resp)
max(gam_CHLOROA_newdata$right_upr) # Nice too
max(gam_CHLOROA_newdata$fit_resp)
# But we can see that the confidence interval is ridiculously small around the
# min and max of the model fit.

# Plot

ggplot(gam_CHLOROA_newdata, aes(x = Day, y = fit_resp))+
  geom_ribbon(aes(x=Day, ymin=right_lwr, ymax=right_upr), fill = 'grey70', alpha=0.7) +
  geom_line(linewidth = 1) +
  geom_point(data = Table_CHLOROA, aes(x = Day, y = CHLOROA, color = Year), shape = 21) +
  scale_color_viridis_c('Year') +
  facet_wrap(facets = c('Code_point_Libelle')) +
  theme_classic()

# Saving plot
# ggsave('gam_CHLOROA_all.tiff', dpi = 300, height = 120, width = 160, 
#        units = 'mm', compression = 'lzw')

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

### GAM of Alloxanthin ####

gam_Allo <- gam(data = Table_Allo, 
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
                   # Using a Gaussian distribution
                   family = gaussian(),
                   # Restricted maximum likelihood estimation (recommended method)
                   method = 'REML')

summary(gam_Allo)
gam.check(gam_Allo)

# gam.check indicates that there might be an issue with the random effect
# smooth for 'Year' (low p-value)

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
# Check the confidence interval. It should not extend below 0 (negative counts
# are impossible)
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
# Some trumpet shaped plots here...

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
# x-axis range not adapted

# Save the plot
# ggsave('HistRes_custom_16sites.tiff', dpi = 300, height = 175, width = 250,
#                units = 'mm', compression = 'lzw')

# Plotting the whole model

ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Allo_newdata, aes(x = Day, ymin = right_lwr, 
                                              ymax = right_upr,
                                              color = Code_point_Libelle, 
                                              fill = Code_point_Libelle),
              linewidth = .75, alpha = .2) +
  geom_line(data = gam_Allo_newdata, aes(x = Day, y = fit_resp,
                                            color = Code_point_Libelle),
            linewidth = 1) +
  # plot the data
  geom_point(data = Table_Allo, aes(x = Day, y = Allo,
                                       color = Code_point_Libelle, 
                                       fill = Code_point_Libelle), size = 2,
             alpha = .3) +
  scale_color_discrete(type = pheno_palette4, guide = 'none') +
  scale_fill_discrete(type = pheno_palette4, guide = 'none') +
  facet_wrap(facets = c('Code_point_Libelle')) +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(title = 'Alloxanthin concentration', x = 'Calendar day', 
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
gam_Dino <- read.csv2('response_pred_plot_20240813_14sites.csv',
                                header = TRUE, fileEncoding = 'ISO-8859-1')
# Separate sites
gam_Dino_Antifer <- filter(gam_Dino, 
                           Code_point_Libelle == 'Antifer ponton pétrolier')
gam_Dino_Cabourg <- filter(gam_Dino, Code_point_Libelle == 'Cabourg')
gam_Dino_MeR <- filter(gam_Dino, Code_point_Libelle == 'Men er Roue')
gam_Dino_OL <- filter(gam_Dino, Code_point_Libelle == 'Ouest Loscolo')

## Mesodinium GAM data
gam_Meso <- read.csv2('pred_plot_MESO_20240418.csv', header = TRUE,
                       fileEncoding = 'ISO-8859-1')

# Separate sites
gam_Meso_Antifer <- filter(gam_Meso, 
                           Code_point_Libelle == 'Antifer ponton pétrolier')
gam_Meso_Cabourg <- filter(gam_Meso, Code_point_Libelle == 'Cabourg')
gam_Meso_MeR <- filter(gam_Meso, Code_point_Libelle == 'Men er Roue')
gam_Meso_OL <- filter(gam_Meso, Code_point_Libelle == 'Ouest Loscolo')

# Count data ####
# The sites must be entered as factors
Season_Dino <- Season_Dino %>%
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

# Separate sites
Season_Meso_Antifer <- filter(Season_Meso, 
                              Code_point_Libelle == 'Antifer ponton pétrolier')
Season_Meso_Cabourg <- filter(Season_Meso, Code_point_Libelle == 'Cabourg')
Season_Meso_MeR <- filter(Season_Meso, Code_point_Libelle == 'Men er Roue')
Season_Meso_OL <- filter(Season_Meso, Code_point_Libelle == 'Ouest Loscolo')



### Plots ####
## Antifer ####

plot_Allo_Antifer <- ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Allo_newdata, aes(x = Day, ymin = right_lwr, 
                                              ymax = right_upr),
              linewidth = .75, alpha = .2,
              color = 'red3', fill = 'red3') +
  geom_line(data = gam_Allo_newdata, aes(x = Day, y = fit_resp),
            linewidth = 1, color = 'red3') +
  # plot the data
  geom_point(data = Table_Allo, aes(x = Day, y = Allo), size = 2,
             alpha = .3, color = 'red3') +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(title = 'Alloxanthin concentration - Antifer', x = 'Calendar day', 
       y = 'Pigment concentration (microgram/L)') +
  theme_classic()

plot_Allo_Antifer
