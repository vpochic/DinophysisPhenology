### A Dinophysis GAM at Basse Michaud ###
# V. POCHIC - 2024-15-11


#### Packages ####
library(tidyverse)
library(ggplot2)
library(mgcv)
library(gratia)
library(viridis)
library(nlme)
library(ggpubr)
library(grid)



#### Data ####

Season_Dino_BM <- read.csv2('Season_Dino_BM.csv', header = TRUE, 
                            fileEncoding = 'ISO-8859-1')
plotDino_BM <- filter(Season_Dino_BM, (Code_point_Libelle == 'Basse Michaud'))
plotDino_OL_select <- filter(Season_Dino_BM, (Code_point_Libelle == 'Ouest Loscolo')
                             & (Year >= 2016))

#### Dinophysis GAM at Basse Michaud ####

gam_Dino_BM <- gam(data = plotDino_BM, 
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

summary(gam_Dino_BM)
gam.check(gam_Dino_BM)

# gam.check indicates that there might be an issue with the random effect
# smooth for 'Year' (low p-value)

# Create a new 'empty' dataset for storing model prediction (fit)
gam_Dino_BM_newdata <- expand_grid(Day=seq(1, 365),
                                   # We add a Year vector as it has become a factor
                                   # of the model
                                   Year=seq(min(Season_Dino_BM$Year), 
                                            max(Season_Dino_BM$Year)))

## Get the inverse link function of the model
# With this function, we can transform the prediction we make on the link
#scale to the response scale
ilink <- gam_Dino_BM$family$linkinv

## Predict : add fit and se.fit on the **link** scale
gam_Dino_BM_newdata <- bind_cols(gam_Dino_BM_newdata, 
                                 setNames(as_tibble(
                                   predict(gam_Dino_BM, gam_Dino_BM_newdata, 
                                           se.fit = TRUE, type = 'link')[1:2]),
                                   c('fit_link','se_link')))

## Create the 95% confidence interval (2*standard error fit) and backtransform 
# to response variable using the inverse link function
gam_Dino_BM_newdata <- mutate(gam_Dino_BM_newdata,
                              fit_resp  = ilink(fit_link),
                              right_upr = ilink(fit_link + (2 * se_link)),
                              right_lwr = ilink(fit_link - (2 * se_link)))
# Check the confidence interval. It should not extend below 0 (negative counts
# are impossible)
min(gam_Dino_BM_newdata$right_lwr) # Nice :)
min(gam_Dino_BM_newdata$fit_resp)
max(gam_Dino_BM_newdata$right_upr) # Nice too
max(gam_Dino_BM_newdata$fit_resp)
# But we can see that the confidence interval is ridiculously small around the
# min and max of the model fit.

# Plot

ggplot(gam_Dino_BM_newdata, aes(x = Day, y = fit_resp))+
  geom_ribbon(aes(x=Day, ymin=right_lwr, ymax=right_upr), fill = 'grey70', alpha=0.7) +
  geom_line(linewidth = 1) +
  geom_point(data = plotDino_BM, aes(x = Day, y = true_count, color = Year), shape = 21) +
  scale_color_viridis_c('Year') +
  theme_classic()

### Checking the model
ModelOutputs<-data.frame(Fitted=fitted(gam_Dino_BM),
                         Residuals=resid(gam_Dino_BM))

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- gam_Dino_BM$model
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
gam_Dino_BM_newdata <- gam_Dino_BM_newdata %>%
  # Create a month and day character variable, that is always '01-01'
  mutate(MonthDay = '01-01') %>%
  # Then unite it with Year to create something that will resemble a ymd date
  mutate(CharYear = as.character(Year)) %>%
  unite(Date, CharYear, MonthDay, sep = '-') %>%
  # And make it a 'Date' object
  mutate(Date = ymd(Date)) %>%
  # Now, add Day-1 to each
  mutate(Date = Date + (as.numeric(Day)-1))

# Basse Michaud
ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Dino_BM_newdata, aes(x = Day, ymin = right_lwr, 
                                              ymax = right_upr),
              color = 'orangered',
              fill = 'orange',
              linewidth = .75, alpha = .8) +
  geom_line(data = gam_Dino_BM_newdata, aes(x = Day, y = fit_resp),
            linewidth = 1,
            color = 'darkred') +
  # plot the data
  geom_point(data = plotDino_BM, aes(x = Day, y = true_count),
             color = 'firebrick1', size = 4,
             alpha = .3) +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(title = 'Basse Michaud', x = 'Date', y = 'Dinophysis cells observed in 10 mL') +
  theme_classic()

## Ouest Loscolo for comparison ####

gam_Dino_OL <- gam(data = plotDino_OL_select, 
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
                                   Year=seq(min(Season_Dino_BM$Year), 
                                            max(Season_Dino_BM$Year)))

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
  geom_point(data = plotDino_OL_select, aes(x = Day, y = true_count, color = Year), shape = 21) +
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

ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Dino_OL_newdata, aes(x = Day, ymin = right_lwr, 
                                              ymax = right_upr),
              color = '#7F96B6',
              fill = '#BBD4F2',
              linewidth = .75, alpha = .8) +
  geom_line(data = gam_Dino_OL_newdata, aes(x = Day, y = fit_resp),
            linewidth = 1,
            color = '#435E7B') +
  # plot the data
  geom_point(data = plotDino_OL_select, aes(x = Day, y = true_count),
             color = '#11203E', size = 4,
             alpha = .3) +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(title = 'Ouest Loscolo (2016-2022)', x = NULL, y = NULL) +
  theme_classic()

### Plots ####

# Basse Michaud
plotgam_BM <- ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Dino_BM_newdata, aes(x = Day, ymin = right_lwr, 
                                              ymax = right_upr),
              color = 'orangered',
              fill = 'orange',
              linewidth = .75, alpha = .8) +
  geom_line(data = gam_Dino_BM_newdata, aes(x = Day, y = fit_resp),
            linewidth = 1,
            color = 'darkred') +
  # plot the data
  geom_point(data = plotDino_BM, aes(x = Day, y = true_count),
             color = 'firebrick1', size = 4,
             alpha = .3) +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(title = 'Basse Michaud (2016-2022)', x = NULL, y = NULL) +
  theme_classic()

plotgam_BM

# Ouest Loscolo
plotgam_OL <- ggplot() +
  # plot the GAM
  geom_ribbon(data = gam_Dino_OL_newdata, aes(x = Day, ymin = right_lwr, 
                                              ymax = right_upr),
              color = '#7F96B6',
              fill = '#BBD4F2',
              linewidth = .75, alpha = .8) +
  geom_line(data = gam_Dino_OL_newdata, aes(x = Day, y = fit_resp),
            linewidth = 1,
            color = '#435E7B') +
  # plot the data
  geom_point(data = plotDino_OL_select, aes(x = Day, y = true_count),
             color = '#11203E', size = 4,
             alpha = .3) +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(title = 'Ouest Loscolo (2016-2022)', x = NULL, y = NULL) +
  theme_classic()

plotgam_OL

# Both sites together
plotgam_both <- ggplot() +
  # plot the GAM
  # OL
  geom_ribbon(data = gam_Dino_OL_newdata, aes(x = Day, ymin = right_lwr, 
                                              ymax = right_upr),
              color = '#7F96B6',
              fill = '#BBD4F2',
              linewidth = .75, alpha = .8) +
  geom_line(data = gam_Dino_OL_newdata, aes(x = Day, y = fit_resp),
            linewidth = 1,
            color = '#435E7B') +
  # BM
  geom_ribbon(data = gam_Dino_BM_newdata, aes(x = Day, ymin = right_lwr, 
                                              ymax = right_upr),
              color = 'orangered',
              fill = 'orange',
              linewidth = .75, alpha = .8) +
  geom_line(data = gam_Dino_BM_newdata, aes(x = Day, y = fit_resp),
            linewidth = 1,
            color = 'darkred') +
  # cut the y scale at 22
  #scale_y_continuous(limits = c(0,22)) +
  # Text
  labs(title = 'Both sites (2016-2022)', x = NULL, y = NULL) +
  theme_classic()

plotgam_both

comparison_plot <- ggarrange(plotgam_OL, plotgam_BM, plotgam_both, nrow = 3, align = 'v')

annotate_figure(comparison_plot, 
                left = textGrob('Dinophysis cells observed in 10 mL', rot = 90, 
                                vjust = 1, gp = gpar(fontsize = 12)),
                bottom = textGrob('Day of the year', gp = gpar(fontsize = 12)))

# ggsave('Comparison_plot_OLBM.tiff', height = 350, width = 164,
#        dpi = 300, unit = 'mm', compression = 'lzw')
