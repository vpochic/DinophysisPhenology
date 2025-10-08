###### Dinophysis phenology: visualisations on multiple years ###
## V. POCHIC
# 2025-10-03

# Additional visualisations using the results of GAM, REPHY data and satellite
# observations, to observe Dinophysis dynamics over several years

#### Packages and functions ####
library(tidyverse)
library(ggplot2)
library(ggnewscale)
library(viridis)
library(cmocean)
library(RColorBrewer)

#### Import data ####

# Table with Dinophysis counts (with added zeros)
Season_Dino <- read.csv2('Data/REPHY_outputs/Season_Dino.csv', 
                         fileEncoding = "ISO-8859-1") %>%
  filter(Year >= 2007 & Year <= 2022)

# Output of Dinophysis GAM
Dino_response_pred <- read.csv2('Data/GAM_outputs/response_pred_GAMDino_20241129.csv', 
                                header = TRUE, fileEncoding = 'ISO-8859-1')

#### Evolution of the timing of maximum count over the study period ####

# One question is: did the phenology of Dinophysis evolve over the study
# period (2007-2022)?
# The GAM we designed is not adequate to answer this question. A **very basic** 
# way to investigate it is to look at the date of the maximum for each site, for
# each year, and check if there is a trend

Season_Dino_nozeros <- Season_Dino %>%
  filter(true_count != 0)

Maxima_Dino <- Season_Dino_nozeros %>%
  # group
  group_by(Code_point_Libelle, Year) %>%
  # extract rows with yearly maxima of 'true_count' (slice() retains groups)
  slice(which.max(true_count)) %>%
  # Code_point_Libelle as factor
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

# Getting some stats
Maxima_Dino_stats <- Maxima_Dino %>%
  group_by(Code_point_Libelle) %>%
  summarise(median_daymax = median(Day), mean_daymax = mean(Day), 
            stdev_daymax = sd(Day), min_daymax = min(Day), max_daymax = max(Day),
            median.lat = median(Latitude),
            .groups = 'keep')



# Color palettes
pheno_palette16 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646',
                     '#642C3A', '#DEB1CC', '#FC4D6B', '#791D40')

# For the sites with regular Dinophysis blooms only
pheno_palette8 <- c('red3', 'orangered', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646')

ggplot(Maxima_Dino, aes(x = Year, y = Day, color = Code_point_Libelle)) +
  # add linear regression fit before the points
  geom_smooth(method = 'lm', alpha = .3) +
  # then the points
  geom_point(size = 4, alpha = .8) +
  # aesthetics
  facet_wrap(facets = 'Code_point_Libelle') +
  scale_color_discrete(type = pheno_palette16, guide = 'none') +
  scale_y_continuous(limits = c(1, 365), 
                     breaks = c(1, 100, 200, 300, 365)) +
  labs(y = 'Day of maximum Dinophysis count') +
  theme_classic()

# That's dope
# ggsave('Maxima_plot_16sites.tiff', dpi = 300, height = 175, width = 250,
#                units = 'mm', compression = 'lzw')

## Let's try something else : for sites who have 2 peaks in their phenology
# (all except Seine Bay sites), let's split the year in 2 to get both peaks
# Our splitting point will be the minimum of the GAM between the 2 peaks

# We're going to identify the turning point in the GAM: the Day for which 
# fit(Day)<fit(Day-1) AND fit(Day)<fit(Day+1)

# Because the phenology of the GAM should be the same for all years, we should
# get the same date (at each site) whatever the year we consider. But we're 
# going to check that.

GAM_pit <- Dino_response_pred %>%
  # Group by site and Year
  group_by(Code_point_Libelle, Year) %>%
  # Filter rows for which the value of fit is inferior to both the previous
  # (lag) and following (lead) value
  filter(lag(fit) > fit & lead(fit) > fit)

# Here we do the opposite: we're looking for the peaks
GAM_peak <- Dino_response_pred %>%
  # Group by site and Year
  group_by(Code_point_Libelle, Year) %>%
  # Filter rows for which the value of fit is inferior to both the previous
  # (lag) and following (lead) value
  filter(lag(fit) < fit & lead(fit) < fit)

# Cleaning up
GAM_pit_select <- filter(GAM_pit, Code_point_Libelle %in% 
                         c('Men er Roue', 'Ouest Loscolo',
                         'Le Cornard', 'Auger',
                         'Arcachon - Bouée 7', 'Teychan bis',
                         'Bouzigues (a)', 'Parc Leucate 2', 'Sète mer')) %>%
# We see there are 2 pits for many sites: 1 corresponds to the pit between the 
# 2 peaks, and 1 to the "renewal" after the winter. We will select only the 
# former (Day > 70 and < 300). This works for all sites except for Parc Leucate.
# We create a special condition for Parc Leucate (pit at day 55)
  filter(ifelse(Code_point_Libelle != 'Parc Leucate 2', Day > 70 & Day < 300,
                Day == 55)) %>%
  # Reorder the site as a factor
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Bouzigues (a)', 'Parc Leucate 2',
                                          'Sète mer'))

GAM_peak_select <- filter(GAM_peak, Code_point_Libelle %in% 
                           c('Antifer ponton pétrolier', 'Cabourg',
                             'Men er Roue', 'Ouest Loscolo',
                             'Le Cornard', 'Auger',
                             'Arcachon - Bouée 7', 'Teychan bis',
                             'Bouzigues (a)', 'Parc Leucate 2',
                             'Sète mer')) %>%
  
  # Reorder the site as a factor
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Bouzigues (a)', 'Parc Leucate 2',
                                          'Sète mer'))

# Next step: split the maxima in 2 periods (before and after the pit)

# True data
Maxima_Dino_2 <- Season_Dino_nozeros %>%
  # group
  group_by(Code_point_Libelle, Year) %>%
  # We'll need to create a period variable (1/2) for each site
  # This piece of code is not pretty but it works
  mutate(period = 
           #Men er Roue
           ifelse(Code_point_Libelle == 'Men er Roue' & Day <= 223,
                  1,
                  ifelse(Code_point_Libelle == 'Men er Roue',
                         2,
           # Ouest Loscolo
           ifelse(Code_point_Libelle == 'Ouest Loscolo' & Day <= 232,
                  1,
                  ifelse(Code_point_Libelle == 'Ouest Loscolo',
                         2,
           # Le Cornard
           ifelse(Code_point_Libelle == 'Le Cornard' & Day <= 210,
                  1,
                  ifelse(Code_point_Libelle == 'Le Cornard',
                         2,
           # Auger
           ifelse(Code_point_Libelle == 'Auger' & Day <= 209,
                  1,
                  ifelse(Code_point_Libelle == 'Auger',
                         2,
           # Arcachon Bouée 7
           ifelse(Code_point_Libelle == 'Arcachon - Bouée 7' & Day <= 195,
                  1,
                  ifelse(Code_point_Libelle == 'Arcachon - Bouée 7',
                         2,
           # Teychan bis
           ifelse(Code_point_Libelle == 'Teychan bis' & Day <= 201,
                  1,
                  ifelse(Code_point_Libelle == 'Teychan bis',
                         2,
           # Parc Leucate 2
           ifelse(Code_point_Libelle == 'Parc Leucate 2' & Day <= 234 
                  # Need to add this condition because of winter peak continuing
                  # into the next year
                  & Day >= 55,
                  1,
                  ifelse(Code_point_Libelle == 'Parc Leucate 2',
                         2,
          # Bouzigues (a)
           ifelse(Code_point_Libelle == 'Bouzigues (a)' & Day <= 232,
                  1,
                  ifelse(Code_point_Libelle == 'Bouzigues (a)',
                         2,
          # Sète mer
           ifelse(Code_point_Libelle == 'Sète mer' & Day <= 209,
                  1,
                  ifelse(Code_point_Libelle == 'Sète mer',
                         2,
          # Diana centre
           ifelse(Code_point_Libelle == 'Diana centre' & Day <= 164,
                  1,
                  ifelse(Code_point_Libelle == 'Diana centre',
                         2,
           # All other sites
           1))))))))))))))))))))) %>%
  # Add period as a grouping variable
  group_by(Code_point_Libelle, Year, period) %>%
  # extract rows with yearly maxima of 'true_count' (slice() retains groups)
  slice(which.max(true_count)) %>%
  # Code_point_Libelle as factor
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
# Keep only Dinophysis cell densities >= X00 cells/L (more than X cells counted)
filter(true_count >= 1)

# GAM fit
GAM_peak_select2 <- GAM_peak_select %>%
  # group
  group_by(Code_point_Libelle) %>%
  # We'll need to create a period variable (1/2) for each site, based on the
  # pits of the GAM
  # This piece of code is not pretty but it works
  mutate(period = 
           #Men er Roue
           ifelse(Code_point_Libelle == 'Men er Roue' & Day <= 223,
                  1,
                  ifelse(Code_point_Libelle == 'Men er Roue',
                         2,
                         # Ouest Loscolo
                         ifelse(Code_point_Libelle == 'Ouest Loscolo' & Day <= 232,
                                1,
                                ifelse(Code_point_Libelle == 'Ouest Loscolo',
                                       2,
                                       # Le Cornard
                                       ifelse(Code_point_Libelle == 'Le Cornard' & Day <= 210,
                                              1,
                                              ifelse(Code_point_Libelle == 'Le Cornard',
                                                     2,
                                                     # Auger
                                                     ifelse(Code_point_Libelle == 'Auger' & Day <= 209,
                                                            1,
                                                            ifelse(Code_point_Libelle == 'Auger',
                                                                   2,
                                                                   # Arcachon Bouée 7
                                                                   ifelse(Code_point_Libelle == 'Arcachon - Bouée 7' & Day <= 195,
                                                                          1,
                                                                          ifelse(Code_point_Libelle == 'Arcachon - Bouée 7',
                                                                                 2,
                                                                                 # Teychan bis
                                                                                 ifelse(Code_point_Libelle == 'Teychan bis' & Day <= 201,
                                                                                        1,
                                                                                        ifelse(Code_point_Libelle == 'Teychan bis',
                                                                                               2,
                                                                                               # All other sites
                                                                                               1))))))))))))) %>%
  # Add period as a grouping variable
  group_by(Code_point_Libelle, Year, period) %>%
  # Code_point_Libelle as factor
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis'))

# Getting some stats for this second thing
Maxima_Dino_2_stats <- Maxima_Dino_2 %>%
  group_by(Code_point_Libelle, period) %>%
  summarise(median_daymax = median(Day), mean_daymax = mean(Day), 
            stdev_daymax = sd(Day), min_daymax = min(Day), max_daymax = max(Day),
            median.lat = median(Latitude),
            .groups = 'keep')

# Let's plot it again

ggplot(data = subset(Maxima_Dino_2,
                     Code_point_Libelle %in% 
                       c('Antifer ponton pétrolier', 'Cabourg',
                         'Men er Roue', 'Ouest Loscolo',
                         'Le Cornard', 'Auger',
                         'Arcachon - Bouée 7', 'Teychan bis')), 
       aes(x = Year, y = Day, color = Code_point_Libelle, 
           shape = as.factor(period))) +
  # add linear regression fit before the points
  geom_smooth(method = 'lm', alpha = .3) +
  # then the points
  geom_point(size = 4, alpha = .8) +
  # aesthetics
  facet_wrap(facets = 'Code_point_Libelle', nrow = 2) +
  scale_color_discrete(type = pheno_palette8, guide = 'none') +
  scale_y_continuous(limits = c(1, 365), 
                     breaks = c(1, 100, 200, 300, 365)) +
  scale_shape_discrete(guide = 'none') +
  # scale_shape_manual(values = c(1, 3), guide = 'none') +
  # labels
  labs(y = 'Day of maximum Dinophysis count') +
  theme_classic()

# That's even doper
# ggsave('Plots/REPHY/Maxima_plot2_8sites_3.tiff', dpi = 300, height = 125, width = 164,
#                units = 'mm', compression = 'lzw')

### Why stop at 8? Let's do the 16 sites!

# New color scale for only the sites with smooths
pheno_palette12 <- c('red3', 'orangered', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646',
                     '#642C3A', '#DEB1CC', '#FC4D6B', '#791D40')


# plot
ggplot(Maxima_Dino_2, 
       aes(x = Year, y = Day, color = Code_point_Libelle, 
           shape = as.factor(period))) +
  # then the points
  geom_point(size = 4, alpha = .8) +
  # add linear regression fit before the points
  geom_smooth(data = subset(Maxima_Dino_2,
                            Code_point_Libelle %in% 
                              c('Antifer ponton pétrolier', 'Cabourg',
                                'Men er Roue', 'Ouest Loscolo',
                                'Le Cornard', 'Auger',
                                'Arcachon - Bouée 7', 'Teychan bis',
                                'Parc Leucate 2', 'Bouzigues (a)', 'Sète mer',
                                'Diana centre')),
                            aes(x = Year, y = Day, color = Code_point_Libelle),
                            method = 'lm', alpha = .18) +
  # new color scale for the points
  scale_color_discrete(type = pheno_palette16, guide = 'none') +
  # aesthetics
  facet_wrap(facets = 'Code_point_Libelle', nrow = 4) +
  scale_y_continuous(limits = c(1, 365), 
                     breaks = c(1, 100, 200, 300, 365)) +
  scale_shape_discrete(guide = 'none') +
  # scale_shape_manual(values = c(1, 3), guide = 'none') +
  # labels
  labs(y = 'Day of maximum Dinophysis count') +
  theme_classic()

# We have some problems ar Parc Leucate 2: the maxima are cut in 2 groups
# because they occur around the 1st of January.
# To fix this, we need to "cut" the seasonality so that maxima in winter don't get 
# separated artificially. Let's do this by setting the first day of our new year
# as day 25.

# UPDATE : with the modification of the code on 2025/08/12, this is no longer a 
# viable solution, because it fcks up the other Med sites.
# Hence, the NewDay variable creation is deactivated (see just below). The
# "solution" is to NOT fit a linear model into the winter peak at Parc Leucate
# (see plot command)
# VP

# Maxima_Dino_2 <- Maxima_Dino_2 %>%
#   # We'll go from 36 to 401 (in our new system)
#   mutate(NewDay = ifelse(Day <= 35, Day + 365, Day))

# Same plot as before but with NewDay as x
ggplot(Maxima_Dino_2, 
       aes(x = Year, y = Day, color = Code_point_Libelle, 
           shape = as.factor(period))) +
  # First the points
  geom_point(size = 2.5, alpha = .8) +
  # then add linear regression fit for some sites only
  geom_smooth(data = subset(Maxima_Dino_2,
                            Code_point_Libelle %in% 
                              c('Antifer ponton pétrolier', 'Cabourg',
                                'Men er Roue', 'Ouest Loscolo',
                                'Le Cornard', 'Auger',
                                'Arcachon - Bouée 7', 'Teychan bis',
                                'Bouzigues (a)', 'Sète mer',
                                'Diana centre')),
              aes(x = Year, y = Day, color = Code_point_Libelle),
              method = 'lm', alpha = .25, linewidth = .5) +
  # Special one for Parc Leucate
  geom_smooth(data = subset(Maxima_Dino_2,
                            Code_point_Libelle == 'Parc Leucate 2' &
                              period == 1),
              aes(x = Year, y = Day),color = '#642C3A',
              method = 'lm', alpha = .25, linewidth = .5) +
  # new color scale for the points
  scale_color_discrete(type = pheno_palette16, guide = 'none') +
  # aesthetics
  facet_wrap(facets = 'Code_point_Libelle', nrow = 4) +
  scale_y_continuous(limits = c(1, 365), 
                     breaks = c(1, 100, 200, 300, 365)) +
  scale_shape_discrete(guide = 'none') +
  # scale_shape_manual(values = c(1, 3), guide = 'none') +
  # labels
  labs(y = 'Day of maximum Dinophysis count (DOY)') +
  theme_classic()

# This seems pretty good!
# Let's save that
# ggsave('Plots/REPHY/Maxima_Dinophysis_16_sites.tiff',
#        dpi = 300, height = 125, width = 164, units = 'mm', compression = 'lzw')

#### Hovmoller diagrams ####

# We'll try to plot Hovmoller diagrams of Dinophysis presence for our sites, 
# to have a point of comparison with the CPR data analysis.

Season_Dino_monthly <- Season_Dino %>%
  group_by(Year, Month, Code_point_Libelle) %>%
  # Create a dummy variable for Dinophysis presence
  mutate(Dinodummy = ifelse(Dinophysis_genus > 0, 1, 0)) %>%
  # Summarise as desired values
  summarise(mean_Dino = mean(Dinophysis_genus),
            # Total sum of Dinophysis counted
            sumDino = sum(Dinophysis_genus),
            # Number of samples with Dinophysis
            nsamples_Dino = sum(Dinodummy),
            # Total number of samples
            nsamples = n(),
            .groups = 'keep') %>%
  # Compute the proportion of Dinophysis
  mutate(prop_pos = nsamples_Dino/nsamples) %>%
  # reorder the sampling site variable
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

## Plots

# Hovmoller diagram of sampling effort
ggplot(data = subset(Season_Dino_monthly, Code_point_Libelle != 0)) +
  geom_tile(aes(x = Year, y = Month, fill = as_factor(nsamples))) +
  scale_fill_viridis(option = 'plasma', discrete = TRUE) +
  facet_wrap(facets = 'Code_point_Libelle') +
  scale_y_continuous(limits = c(0, 13), 
                     breaks = c(1,2,3,4,5,6,7,8,9,10,11,12)) +
  labs(title = 'Sampling effort over the study period,
by sampling site',
       fill = 'Number of samples (log10)') +
  theme_classic() +
  theme(legend.position = 'bottom')

# Save plot
# ggsave('Plots/REPHY/REPHY_Hovmoller_sampling.tiff', height = 225, width = 250, units = 'mm',
#        dpi = 300, compression = 'lzw')

### Then, Hovmoller diagram of Dinophysis presence

ggplot(data = subset(Season_Dino_monthly, Code_point_Libelle != 0)) +
  geom_tile(aes(x = Year, y = Month, fill = log(prop_pos+1))) +
  scale_fill_viridis() +
  facet_wrap(facets = 'Code_point_Libelle') +
  scale_y_continuous(limits = c(0, 13), 
                     breaks = c(1,2,3,4,5,6,7,8,9,10,11,12)) +
  labs(title = 'Dinophysis seasonality over the study period, 
by sampling site',
       fill = 'Dinophysis presence 
(proportion of samples, log scale)') +
  theme_classic() +
  theme(legend.position = 'bottom')

# This one's pretty awful.
# ggsave('Plots/REPHY/REPHY_Hovmoller_proportion.tiff', height = 225, width = 250, units = 'mm',
#        dpi = 300, compression = 'lzw')

### Finally, Hovmoller diagram of Dinophysis abundance

ggplot(data = subset(Season_Dino_monthly, Code_point_Libelle != 0)) +
  geom_tile(aes(x = Year, y = Month, fill = log10(mean_Dino+1))) +
  scale_fill_viridis() +
  facet_wrap(facets = 'Code_point_Libelle') +
  scale_y_continuous(limits = c(0, 13), 
                     breaks = c(1,2,3,4,5,6,7,8,9,10,11,12)) +
  labs(title = 'Dinophysis seasonality over the study period, 
by sampling site',
       fill = 'Mean Dinophysis abundance (log10)') +
  theme_classic() +
  theme(legend.position = 'bottom')

# This one's a bit better but still not great.
# ggsave('Plots/REPHY/REPHY_Hovmoller_abundance.tiff', height = 225, width = 250, units = 'mm',
#        dpi = 300, compression = 'lzw')

## Hovmoller diagrams - CPR timeline ####

# We elongate the x-axis so we have the same time range than the CPR data

# Antifer
ggplot(data = subset(Season_Dino_monthly, Code_point_Libelle == 'Antifer ponton pétrolier')) +
  geom_tile(aes(x = Year, y = Month, fill = log10(mean_Dino+1))) +
  scale_fill_viridis() +
  scale_x_continuous(limits = c(1958, 2023), 
                     breaks = c(1960, 2000, 2007, 2020)) +
  scale_y_continuous(limits = c(0, 13), 
                     breaks = c(1,2,3,4,5,6,7,8,9,10,11,12)) +
  labs(title = 'Antifer ponton pétrolier',
       fill = 'Mean Dinophysis abundance (log10)') +
  theme_classic() +
  theme(legend.position = 'bottom')

# Cabourg
ggplot(data = subset(Season_Dino_monthly, Code_point_Libelle == 'Cabourg')) +
  geom_tile(aes(x = Year, y = Month, fill = log10(mean_Dino+1))) +
  scale_fill_viridis() +
  scale_x_continuous(limits = c(1958, 2023), 
                     breaks = c(1960, 2000, 2010, 2020)) +
  scale_y_continuous(limits = c(0, 13), 
                     breaks = c(1,2,3,4,5,6,7,8,9,10,11,12)) +
  labs(title = 'Cabourg',
       fill = 'Mean Dinophysis abundance (log10)') +
  theme_classic() +
  theme(legend.position = 'bottom')

# Men er Roue
ggplot(data = subset(Season_Dino_monthly, Code_point_Libelle == 'Men er Roue')) +
  geom_tile(aes(x = Year, y = Month, fill = log10(mean_Dino+1))) +
  scale_fill_viridis() +
  scale_x_continuous(limits = c(1958, 2023), 
                     breaks = c(1960, 2000, 2010, 2020)) +
  scale_y_continuous(limits = c(0, 13), 
                     breaks = c(1,2,3,4,5,6,7,8,9,10,11,12)) +
  labs(title = 'Men er Roue',
       fill = 'Mean Dinophysis abundance (log10)') +
  theme_classic() +
  theme(legend.position = 'bottom')

# Ouest Loscolo
ggplot(data = subset(Season_Dino_monthly, Code_point_Libelle == 'Ouest Loscolo')) +
  geom_tile(aes(x = Year, y = Month, fill = log10(mean_Dino+1))) +
  scale_fill_viridis() +
  scale_x_continuous(limits = c(1958, 2023), 
                     breaks = c(1960, 2000, 2010, 2020)) +
  scale_y_continuous(limits = c(0, 13), 
                     breaks = c(1,2,3,4,5,6,7,8,9,10,11,12)) +
  labs(title = 'Ouest Loscolo',
       fill = 'Mean Dinophysis abundance (log10)') +
  theme_classic() +
  theme(legend.position = 'bottom')

## That's all very nice, let's group them by region to obtain a broader picture

Season_Dino_region <- Season_Dino %>%
  group_by(Year, Month, Code.Region) %>%
  # Correct an error in the Code.Region
  mutate(Code.Region = ifelse(Code_point_Libelle == 'Ouest Loscolo',
                              21,
                              Code.Region)) %>%
  # Create a dummy variable for Dinophysis presence
  mutate(Dinodummy = ifelse(Dinophysis_genus > 0, 1, 0)) %>%
  # Summarise as desired values
  summarise(mean_Dino = mean(Dinophysis_genus),
            # Total sum of Dinophysis counted
            sumDino = sum(Dinophysis_genus),
            # Number of samples with Dinophysis
            nsamples_Dino = sum(Dinodummy),
            # Total number of samples
            nsamples = n(),
            .groups = 'keep') %>%
  # Compute the proportion of Dinophysis
  mutate(prop_pos = nsamples_Dino/nsamples) %>%
  # Add the name of regions
  mutate(Region_name = ifelse(
    Code.Region == 11, 'Pas de Calais',
    ifelse(Code.Region == 12, 'Seine Bay',
           ifelse(Code.Region == 13, 'Western Channel',
                  ifelse(Code.Region == 21, 'Southern Brittany',
                         ifelse(Code.Region == 22, 'Pertuis',
                                ifelse(Code.Region == 23, 'Arcachon',
                                       ifelse(Code.Region == 31, 'Golfe du Lion',
                                              'Corsica'))))))
  )) %>%
  # And recode it as factor
  mutate(Region_name = as_factor(Region_name)) %>%
  mutate(Region_name = fct_relevel(Region_name,
                                   'Pas de Calais', 'Seine Bay', 'Western Channel',
                                   'Southern Brittany', 'Pertuis', 'Arcachon',
                                   'Golfe du Lion', 'Corsica'))

# Hovmoller diagram : all regions
ggplot(data = Season_Dino_region) +
  geom_tile(aes(x = Year, y = Month, fill = log10(mean_Dino+1))) +
  scale_fill_viridis(breaks = c(log10(0+1), log10(10+1), log10(100+1),
                                log10(1000+1)), labels = c('0','10','100','1000')) +
  scale_x_continuous(limits = c(1958, 2023), 
                     breaks = c(1960, 2007, 2020)) +
  scale_y_continuous(limits = c(0, 13), 
                     breaks = c(1,2,3,4,5,6,7,8,9,10,11,12)) +
  facet_wrap(facets = 'Region_name') +
  labs(title = 'Dinophysis seasonality, REPHY data',
       fill = 'Mean Dinophysis abundance (log10)') +
  theme_classic() +
  theme(legend.position = 'bottom')

# save plot
# ggsave('Plots/REPHY/Hovmoller_abundance_all.tiff', height = 150, width = 164,
#        units = 'mm', compression = 'lzw')

# Only Atlantic regions
ggplot(data = subset(Season_Dino_region, Code.Region < 30)) +
  geom_tile(aes(x = Year, y = Month, fill = log10(mean_Dino+1))) +
  scale_fill_viridis(breaks = c(log10(0+1), log10(10+1), log10(100+1),
                                log10(1000+1)), labels = c('0','10','100','1000')) +
  scale_x_continuous(limits = c(1958, 2023), 
                     breaks = c(1960, 2007, 2020)) +
  scale_y_continuous(limits = c(0, 13), 
                     breaks = c(1,2,3,4,5,6,7,8,9,10,11,12)) +
  facet_wrap(facets = 'Region_name') +
  labs(title = 'Dinophysis seasonality, REPHY data',
       fill = 'Mean Dinophysis abundance (log10)') +
  theme_classic() +
  theme(legend.position = 'bottom')

# save plot
# ggsave('Plots/REPHY/Hovmoller_abundance_atlantic.tiff', height = 125, width = 164,
#        units = 'mm', compression = 'lzw')

### Environmental data ####

### Import data -----------

## Here we will call a table of environmental data to plot a heatmap
# underneath our Dinophysis data
# This table was created by an earlier script

Table_hydro_daily <- read.csv2('Data/REPHY_outputs/Table_hydro_daily_20240821.csv', 
                               header = TRUE,
                               fileEncoding = 'ISO-8859-1') %>%
  # Erase the chl a values that are improper
  select(-c('CHLOROA.med', 'CHLOROA.mean'))

# Now call the fortinightly table with the right values
Table_hydro_fortnightly <- read.csv2('Data/REPHY_outputs/Table_hydro_fortnightly_20250812.csv', 
                                     header = TRUE,
                                     fileEncoding = 'ISO-8859-1')

# Aaaaand... left_join
Table_hydro_daily <- left_join(Table_hydro_daily, Table_hydro_fortnightly)

# Great!

## Now, we will import the stratification data from the GAMAR model (summarised
# by fortnight)
Table_stratif_fortnightly <- read.csv2('Data/Models/GAMAR/Outputs/Table_stratif_fortnightly.csv', 
          header = TRUE, fileEncoding = 'ISO-8859-1')

# Create a daily table
Table_daily <- select(Table_hydro_daily, c('Day', 'Fortnight', 
                                           'Code_point_Libelle'))

Table_stratif_daily <- left_join(Table_daily, Table_stratif_fortnightly,
                                 by = c('Fortnight', 'Code_point_Libelle'))

# Excellent

# Now, we plot! First, what we want to do is to have an idea of the distribution
# of Dinophysis maxima for each site. We will focus on Channel/Atlantic sites,
# and only those where we observe Dinophysis
Table_hydro_daily_select <- filter(Table_hydro_daily,
                                   Code_point_Libelle %in%
                                     c('Antifer ponton pétrolier', 'Cabourg',
                                       'Men er Roue', 'Ouest Loscolo',
                                       'Le Cornard', 'Auger',
                                       'Arcachon - Bouée 7', 'Teychan bis'))

Table_stratif_daily_select <- filter(Table_stratif_daily,
                                   Code_point_Libelle %in%
                                     c('Antifer ponton pétrolier', 'Cabourg',
                                       'Men er Roue', 'Ouest Loscolo',
                                       'Le Cornard', 'Auger',
                                       'Arcachon - Bouée 7', 'Teychan bis'))

Maxima_Dino_2_stats_select <- filter(Maxima_Dino_2_stats,
                                   Code_point_Libelle %in%
                                     c('Antifer ponton pétrolier', 'Cabourg',
                                       'Men er Roue', 'Ouest Loscolo',
                                       'Le Cornard', 'Auger',
                                       'Arcachon - Bouée 7', 'Teychan bis'))
Maxima_Dino_select <- filter(Maxima_Dino_2,
                             Code_point_Libelle %in%
                               c('Antifer ponton pétrolier', 'Cabourg',
                                 'Men er Roue', 'Ouest Loscolo',
                                 'Le Cornard', 'Auger',
                                 'Arcachon - Bouée 7', 'Teychan bis'))

### Plots ------

# New color scale
pheno_palette8 <- c('red3', 'orangered', '#2156A1', '#5995E3', 
                    '#1F3700', '#649003','#F7B41D', '#FBB646')

# First, let's plot the distribution of the maxima (without the heatmap)
ggplot(Maxima_Dino_2_stats_select, aes(color = Code_point_Libelle)) +
  
  # All maxima as smaller translucid points
  geom_point(data = Maxima_Dino_select, aes(x = Day, y = Code_point_Libelle, 
                                            color = Code_point_Libelle), 
             size = 2.5, alpha = .5) +
  
  ## Mean +- sd of maxima (by site and by period)
  # error bars
  geom_errorbar(aes(y = Code_point_Libelle, xmin = mean_daymax - stdev_daymax,
                    xmax = mean_daymax + stdev_daymax), linewidth = .5,
                color = 'grey10', width = .2) +
  
  # points for means
  geom_point(aes(y = Code_point_Libelle, x = mean_daymax, 
                 fill = Code_point_Libelle), size = 4, color = 'grey10',
             shape = 21, stroke = .5) +
  
  # Color scales
  scale_color_discrete(type = pheno_palette8, guide = 'none') +
  scale_fill_discrete(type = pheno_palette8, guide = 'none') +
  
  ## Axis stuff
  scale_x_continuous(limits = c(1,365), breaks = c(1, 100, 200, 300, 365)) +
  # reverse y axis to get sites in the desired order
  scale_y_discrete(limits = rev) +
  
  # labels
  labs(x = 'Day of the year', y = NULL) +
  theme_classic()

# Now we try to underlie a heatmap of temperature/chl a/stratification/ssr seasonality 
# (4 versions of the plot)
ggplot(Maxima_Dino_2_stats_select) +
  geom_tile(data = Table_stratif_daily_select, # Table_hydro_daily_select
            aes(x = Day, y = Code_point_Libelle, fill = SI_date.med), 
            # fill = TEMP.med ; fill = CHLOROA.med ; fill = ssr.med
            alpha = .8) +
  # color palette for fill
  # Stratification version
  scale_fill_cmocean(name = 'deep', direction = 1) +
  # Temperature version
  # scale_fill_distiller(palette = 'RdBu', direction = -1) +
  # Chl a version
  # scale_fill_cmocean(name = 'algae') +
  
  # Add labels here so the name of the 'fill' legend is correct
  # Labels (for TEMP plot fill = 'Median SST (°C)',
  # Chl a: fill = 'Median [chl a] (mg/m3)',
  # Stratif: fill = 'Median Stratification Index (-)')
  labs(x = 'Day of the year', y = NULL, fill = 'Median Stratification Index (-)',
       title = 'Maxima of Dinophysis counts (2007-2022)') +
  
  # New fill scale for points
  new_scale_fill() +
  
  # All maxima as smaller translucid points
  geom_point(data = Maxima_Dino_select, aes(x = Day, y = Code_point_Libelle, 
                                            color = Code_point_Libelle), 
             size = 2.5, alpha = .5) +
  
  ## Mean +- sd of maxima (by site and by period)
  # error bars
  geom_errorbar(aes(y = Code_point_Libelle, xmin = mean_daymax - stdev_daymax,
                    xmax = mean_daymax + stdev_daymax), linewidth = .5,
                color = 'grey10', width = .2) +
  
  # points for means
  geom_point(aes(y = Code_point_Libelle, x = mean_daymax, 
                 fill = Code_point_Libelle), size = 4, color = 'grey10',
             shape = 21, stroke = .5) +
  
  # Color scales
  scale_color_discrete(type = pheno_palette8, guide = 'none') +
  scale_fill_discrete(type = pheno_palette8, guide = 'none') +
  
  # Axis stuff
  # reverse y axis to get sites in the desired order
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(1,365), breaks = c(1, 100, 200, 300, 365)) +
  
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

# Nice one. Let's save that
# ggsave('Plots/REPHY/Dino_max_chla.tiff', dpi = 300, height = 180, width = 150,
#        units = 'mm', compression = 'lzw')

####

## What if instead of mean+-sd, we plot the peaks of the GAM? ##

# Summarise GAM data to remove the year
GAM_peak_plot <- GAM_peak_select2 %>%
  filter(Day < 300) %>%
  group_by(Code_point_Libelle, period) %>%
  summarise(Daymax = mean(Day), .groups = 'keep')

GAM_peak_plot <- left_join(GAM_peak_plot, Maxima_Dino_2_stats_select,
                           by = c('Code_point_Libelle', 'period')) %>%
  filter(Code_point_Libelle %in%
           c('Antifer ponton pétrolier', 'Cabourg',
             'Men er Roue', 'Ouest Loscolo',
             'Le Cornard', 'Auger',
             'Arcachon - Bouée 7', 'Teychan bis'))

# And new version of the plot
ggplot(GAM_peak_plot) +
  geom_tile(data = Table_hydro_daily_select, 
            aes(x = Day, y = Code_point_Libelle, fill = TEMP.med), #fill = CHLOROA.med
            alpha = .8) +
  # color palette for fill
  # Temperature version
  scale_fill_distiller(palette = 'RdBu', direction = -1) +
  # Chl a version
  # scale_fill_cmocean(name = 'algae') +
  
  # Add labels here so the name of the 'fill' legend is correct
  # Labels (for chla plot fill = 'Median [chl a] (mg/m3)')
  labs(x = 'Day of the year', y = NULL, fill = 'Median SST (°C)',
       title = 'Maxima of Dinophysis counts (2007-2022)') +
  
  # New fill scale for points
  new_scale_fill() +
  
  # All maxima as smaller translucid points
  geom_point(data = Maxima_Dino_select, aes(x = Day, y = Code_point_Libelle, 
                                            color = Code_point_Libelle), 
             size = 2.5, alpha = .5) +
  
  # points for GAM peak
  geom_point(aes(y = Code_point_Libelle, x = Daymax, 
                 fill = Code_point_Libelle), size = 4, color = 'grey10',
             shape = 21, stroke = .5) +
  
  # Color scales
  scale_color_discrete(type = pheno_palette8, guide = 'none') +
  scale_fill_discrete(type = pheno_palette8, guide = 'none') +
  
  # Axis stuff
  # reverse y axis to get sites in the desired order
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(1,365), breaks = c(1, 100, 200, 300, 365)) +
  
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

# Great! Let's save that
# ggsave('Plots/REPHY/Dino_gam_max_temp.tiff', dpi = 300, height = 180, width = 150,
#        units = 'mm', compression = 'lzw')

### Alternative : plot only Seine Bay and Southern Brittany sites ###
ggplot(subset(GAM_peak_plot, Code_point_Libelle %in% 
                c('Antifer ponton pétrolier', 'Cabourg', 
                  'Men er Roue', 'Ouest Loscolo'))) +
  geom_tile(data = subset(Table_hydro_daily_select, Code_point_Libelle %in% 
                            c('Antifer ponton pétrolier', 'Cabourg', 
                              'Men er Roue', 'Ouest Loscolo')), 
            aes(x = Day, y = Code_point_Libelle, fill = TEMP.med), #fill = CHLOROA.med
            alpha = .8) +
  # color palette for fill
  # Temperature version
  scale_fill_distiller(palette = 'RdBu', direction = -1) +
  # Chl a version
  # scale_fill_cmocean(name = 'algae') +
  
  # Add labels here so the name of the 'fill' legend is correct
  # Labels (for chla plot fill = 'Median [chl a] (mg/m3)')
  labs(x = 'Day of the year', y = NULL, fill = 'Median SST (°C)',
       title = 'Maxima of Dinophysis counts (2007-2022)') +
  
  # New fill scale for points
  new_scale_fill() +
  
  # All maxima as smaller translucid points
  geom_point(data = subset(Maxima_Dino_select, Code_point_Libelle %in%
                             c('Antifer ponton pétrolier', 'Cabourg', 
                               'Men er Roue', 'Ouest Loscolo')), 
             aes(x = Day, y = Code_point_Libelle, 
                 color = Code_point_Libelle), 
             size = 2.5, alpha = .5) +
  
  # points for GAM peak
  geom_point(aes(y = Code_point_Libelle, x = Daymax, 
                 fill = Code_point_Libelle), size = 4, color = 'grey10',
             shape = 21, stroke = .5) +
  
  # Color scales
  scale_color_discrete(type = pheno_palette8, guide = 'none') +
  scale_fill_discrete(type = pheno_palette8, guide = 'none') +
  
  # Axis stuff
  # reverse y axis to get sites in the desired order
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(1,365), breaks = c(1, 100, 200, 300, 365)) +
  
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

# ggsave('Plots/REPHY/Dino_gam_max_temp_4sites.tiff', dpi = 300, height = 135, width = 150,
#        units = 'mm', compression = 'lzw')


### Plotting the sites used in the Random Forest model ####

Table_hydro_daily_RF <- filter(Table_hydro_daily,
                                   Code_point_Libelle %in%
                                     c('Point 1 Boulogne', 'At so',
                                       'les Hébihens', 'Loguivy',
                                       'Antifer ponton pétrolier', 'Cabourg',
                                       'Men er Roue', 'Ouest Loscolo',
                                       'Le Cornard',
                                       'Arcachon - Bouée 7')) %>%
  # Recode the Code_point_Libelle
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard',
                                          'Arcachon - Bouée 7'))

# Stratification
Table_stratif_daily_RF <- filter(Table_stratif_daily,
                               Code_point_Libelle %in%
                                 c('Point 1 Boulogne', 'At so',
                                   'les Hébihens', 'Loguivy',
                                   'Antifer ponton pétrolier', 'Cabourg',
                                   'Men er Roue', 'Ouest Loscolo',
                                   'Le Cornard',
                                   'Arcachon - Bouée 7')) %>%
  # Recode the Code_point_Libelle
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard',
                                          'Arcachon - Bouée 7'))

# Dinophysis maxima (REPHY data)
Maxima_Dino_2_stats_RF <- filter(Maxima_Dino_2_stats,
                                     Code_point_Libelle %in%
                                       c('Point 1 Boulogne', 'At so',
                                         'les Hébihens', 'Loguivy',
                                         'Antifer ponton pétrolier', 'Cabourg',
                                         'Men er Roue', 'Ouest Loscolo',
                                         'Le Cornard',
                                         'Arcachon - Bouée 7'))
Maxima_Dino_RF <- filter(Maxima_Dino_2,
                             Code_point_Libelle %in%
                               c('Point 1 Boulogne', 'At so',
                                 'les Hébihens', 'Loguivy',
                                 'Antifer ponton pétrolier', 'Cabourg',
                                 'Men er Roue', 'Ouest Loscolo',
                                 'Le Cornard', 'Arcachon - Bouée 7')) %>%
  # Recode the Code_point_Libelle
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'les Hébihens', 'Loguivy',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard',
                                          'Arcachon - Bouée 7'))

# New color scales
pheno_palette10 <- c('sienna4', 'tan3', 'red3', 'orangered', 
                     '#0A1635', '#2B4561', '#2156A1', '#5995E3', 
                     '#1F3700', '#F7B41D')
pheno_palette6 <- c('red3', 'orangered', '#2156A1', '#5995E3', 
                     '#1F3700', '#F7B41D')

# New GAM peak table (without 2 sites)
GAM_peak_plot_RF <- filter(GAM_peak_plot, Code_point_Libelle %in% 
                            c('Antifer ponton pétrolier', 'Cabourg',
                              'Men er Roue', 'Ouest Loscolo',
                              'Le Cornard',
                              'Arcachon - Bouée 7'))

## Plotting the heatmap, maxima and GAM peaks
# 4 versions of the plot : temperature/chl a/stratification/ssr heatmaps
ggplot(Maxima_Dino_2_stats_select) +
  geom_tile(data = Table_stratif_daily_RF, # Table_hydro_daily_RF
            aes(x = Day, y = Code_point_Libelle, fill = SI_date.med), 
            # fill = TEMP.med ; fill = CHLOROA.med ; fill = ssr.med
            alpha = .8) +
  # color palette for fill
  # Stratification version
  scale_fill_cmocean(name = 'deep', direction = 1) +
  # Temperature version
  # scale_fill_distiller(palette = 'RdBu', direction = -1) +
  # Chl a version
  # scale_fill_cmocean(name = 'algae') +
  
  # Add labels here so the name of the 'fill' legend is correct
  # Labels (for TEMP plot fill = 'Median SST (°C)',
  # Chl a: fill = 'Median [chl a] (mg/m3)',
  # Stratif: fill = 'Median Stratification Index (-)')
  labs(x = 'Day of the year', y = NULL, fill = 'Median Stratification Index (-)',
       title = NULL) +
  
  # New fill scale for points
  new_scale_fill() +
  
  # All maxima as smaller translucid points
  geom_point(data = Maxima_Dino_RF, aes(x = Day, y = Code_point_Libelle, 
                                            color = Code_point_Libelle), 
             size = 2.5, alpha = .5) +
  # Color scales (Maxima)
  scale_color_discrete(type = pheno_palette10, guide = 'none') +
  scale_fill_discrete(type = pheno_palette10, guide = 'none') +
  # New color scale
  new_scale_fill() +
  
  # points for GAM peak
  geom_point(data = GAM_peak_plot_RF, aes(y = Code_point_Libelle, x = Daymax, 
                 fill = Code_point_Libelle), size = 4, color = 'grey10',
             shape = 21, stroke = .5) +
  # Color scales (GAM peak)
  scale_fill_discrete(type = pheno_palette6, guide = 'none') +
  
  
  # Axis stuff
  # reverse y axis to get sites in the desired order
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(1,365), breaks = c(1, 100, 200, 300, 365)) +
  
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

# Great! Let's save that
# ggsave('Plots/REPHY/Dino_gam_max_stratif_RF.tiff', dpi = 300, height = 180, width = 100,
#        units = 'mm', compression = 'lzw')

### Plotting maxima by latitude ####

ggplot() +
  
  # All maxima as smaller translucid points
  geom_point(data = Maxima_Dino_select, aes(x = Day, y = Latitude, 
                                            fill = Code_point_Libelle,
                                            shape = as.factor(period),
                                            alpha = log10(true_count)),
             stroke = 0, color = 'transparent',
             size = 2.5) +
  
  ## Mean +- sd of maxima (by site and by period)
  # error bars
  geom_errorbar(data = Maxima_Dino_2_stats_select,
                aes(y = median.lat, xmin = mean_daymax - stdev_daymax,
                    xmax = mean_daymax + stdev_daymax),
                linewidth = .5,
                color = 'grey10', width = .1) +
  
  # points for means
  geom_point(data = Maxima_Dino_2_stats_select,
             aes(y = median.lat, x = mean_daymax, 
                 fill = Code_point_Libelle, shape = as.factor(period)),
             size = 4, color = 'grey10', stroke = .5) +
  
  # Color scales
  scale_color_discrete(type = pheno_palette8, guide = 'none') +
  scale_fill_discrete(type = pheno_palette8, guide = 'none') +
  # Shape
  scale_shape_manual(values = c(21, 24), guide = 'none') +
  # Transparency
  scale_alpha_continuous(guide = 'none') +
  
  ## Axis stuff
  scale_x_continuous(limits = c(1,365), breaks = c(1, 100, 200, 300, 365)) +
  # reverse y axis to get sites in the desired order
  # scale_y_discrete(limits = rev) +
  
  # labels
  labs(x = 'Day of the year', y = 'Latitude') +
  theme_classic() +
  theme(legend.position = 'right')

# Save this plot: maxima depending on latitude
# ggsave('Dino_max_lat.tiff', dpi = 300, height = 180, width = 150,
#        units = 'mm', compression = 'lzw')

## What if instead of mean+-sd, we plot the peaks of the GAM?

# First, let's add the latitude to the GAM data

# Summarise GAM data to remove the year
GAM_peak_plot <- GAM_peak_select2 %>%
  filter(Day < 300) %>%
  group_by(Code_point_Libelle, period) %>%
  summarise(Daymax = mean(Day), .groups = 'keep')

GAM_peak_plot <- left_join(GAM_peak_plot, Maxima_Dino_2_stats_select,
                           by = c('Code_point_Libelle', 'period'))

# Plot this
ggplot() +
  
  # All maxima as smaller translucid points
  geom_point(data = Maxima_Dino_select, aes(x = Day, y = Latitude, 
                                            fill = Code_point_Libelle,
                                            shape = as.factor(period),
                                            alpha = log10(true_count)),
             stroke = 0, color = 'transparent',
             size = 2.5) +
  
  ## Peaks from the GAM (no error bars)
  
  # points for GAM peaks
  geom_point(data = GAM_peak_plot,
             aes(y = median.lat, x = Daymax, 
                 fill = Code_point_Libelle, shape = as.factor(period)),
             size = 4, color = 'grey10', stroke = .5) +
  
  # Color scales
  scale_color_discrete(type = pheno_palette8, guide = 'none') +
  scale_fill_discrete(type = pheno_palette8, guide = 'none') +
  # Shape
  scale_shape_manual(values = c(21, 24), guide = 'none') +
  # Transparency
  scale_alpha_continuous(guide = 'none') +
  
  ## Axis stuff
  scale_x_continuous(limits = c(1,365), breaks = c(1, 100, 200, 300, 365)) +
  # reverse y axis to get sites in the desired order
  # scale_y_discrete(limits = rev) +
  
  # labels
  labs(x = 'Day of the year', y = 'Latitude') +
  theme_classic() +
  theme(legend.position = 'right')

# Save this plot: maxima depending on latitude
# ggsave('Dino_max_lat_GAM_peak.tiff', dpi = 300, height = 180, width = 150,
#        units = 'mm', compression = 'lzw')

### Focusing on Ouest Loscolo - different years ####

# Select only Ouest Loscolo
Season_Dino_OL <- filter(Season_Dino, Code_point_Libelle == 'Ouest Loscolo')

#### Plotting the GAM on the different years ##
## The logic here is to see how well the GAM fits the data year by year
# Can we see some years where the phenology significantly deviates from
# what the GAM predicts?

rpp_OL <- filter(Dino_response_pred, Code_point_Libelle == 'Ouest Loscolo')

ggplot(Season_Dino_OL, aes(x = Day, y = true_count)) +
  geom_point(color = 'dodgerblue3', alpha = .8) +
  geom_line(data = rpp_OL, aes(x = Day, y = fit), color = 'dodgerblue3') +
  geom_ribbon(data = rpp_OL, aes(x = Day, y = fit,
                                 ymin = lwrS, ymax = uprS), 
              color = 'dodgerblue3', alpha = 0.2) +
  facet_wrap(facets = c('Year')) +
  # scale_y_continuous(limits = c(0,30)) +
  theme_classic()

# Ok nice. Now we'll try to add dots for satellite observations of Mesodinium
# blooms in the region.

Satellite_survey <- read.csv2('Data/Satellite/Blooms_list_satellite_Vilaine-Gironde.csv',
                              header = TRUE, fileEncoding = 'ISO-8859-1')

# Let's focus on Mesodinium in the Loire-Atlantique region
Satellite_survey_Meso <- filter(Satellite_survey,
                                (grepl('Mesodinium', Main.species) |
                                  grepl('Mesodinium', Comment)) &
                                  Region == 'Loire Atlantique') %>%
  # Add the date information that we will need for plotting
  mutate(Date = ymd(Date)) %>%
  mutate(Year = year(Date)) %>%
  mutate(Day = yday(Date)) %>%
  # Finally, add a false "true_count" value for plotting alongside the Dinophysis
  # phenology
  mutate(true_count = 50)

# Let's plot again with satellite observations
ggplot(Season_Dino_OL, aes(x = Day, y = true_count)) +
  # Dinophysis GAM +- confidence interval
  geom_line(data = rpp_OL, aes(x = Day, y = fit), color = '#2156A1') +
  geom_ribbon(data = rpp_OL, aes(x = Day, y = fit,
                                 ymin = lwrS, ymax = uprS), 
              color = '#2156A1',
              fill = '#BBD4F2', alpha = 0.2) +
  # Dinophysis counts
  geom_point(color = '#2156A1', alpha = .8) +
  
  # Satellite observations of Mesodinium blooms
  geom_point(data = Satellite_survey_Meso, aes(x = Day, y = true_count),
             color = 'firebrick4', shape = 8, size = 3.5, stroke = .45) +
  
  # 1 subplot for each year
  facet_wrap(facets = c('Year')) +
  # labels
  labs(y = 'Dinophysis cells observed in 10 mL',
       x = 'Julian day', title = 'Ouest Loscolo') +
  
  theme_classic()

# Seems quite nice. let's save this plot.
# ggsave('Ouest_Loscolo_multiyear.tiff', dpi = 300, height = 164, width = 164,
#                units = 'mm', compression = 'lzw')


### Focusing on Basse Michaud - different years ####

# Import Season Dino at basse Michaud
Season_Dino_BM <- read.csv2('Data/REPHY_outputs/Season_Dino_BM_20250401.csv',
                            header = TRUE, fileEncoding = 'ISO-8859-1')

# Import response pred of GAM Dino at Basse Michaud (generated by script
# 'GAM_Dino_OLBM_comparison')
rpp_BM <- read.csv2('Data/GAM_outputs/response_pred_GAMDino_BM_20250401.csv',
                    header = TRUE, fileEncoding = 'ISO-8859-1')

#### Plotting the GAM on the different years ##
## The logic here is to see how well the GAM fits the data year by year
# Can we see some years where the phenology significantly deviates from
# what the GAM predicts?

ggplot(Season_Dino_BM, aes(x = Day, y = true_count)) +
  geom_point(color = '#1E3523', alpha = .8) +
  geom_line(data = rpp_BM, aes(x = Day, y = fit), color = '#1E3523') +
  geom_ribbon(data = rpp_BM, aes(x = Day, y = fit,
                                 ymin = lwrS, ymax = uprS), 
              color = '#1E3523', fill = '#377185', alpha = 0.2) +
  facet_wrap(facets = c('Year')) +
  # scale_y_continuous(limits = c(0,30)) +
  theme_classic()

# Ok nice. Now we'll try to add dots for satellite observations of Mesodinium
# blooms in the region.

Satellite_survey <- read.csv2('Data/Satellite/Blooms_list_satellite_Vilaine-Gironde.csv',
                              header = TRUE, fileEncoding = 'ISO-8859-1')

# Let's focus on Mesodinium in the Loire-Atlantique region
Satellite_survey_Meso <- filter(Satellite_survey,
                                (grepl('Mesodinium', Main.species) |
                                   grepl('Mesodinium', Comment)) &
                                  Region == 'Loire Atlantique') %>%
  # Add the date information that we will need for plotting
  mutate(Date = ymd(Date)) %>%
  mutate(Year = year(Date)) %>%
  mutate(Day = yday(Date)) %>%
  # Finally, add a false "true_count" value for plotting alongside the Dinophysis
  # phenology
  mutate(true_count = 25)

# Let's plot again with satellite observations
ggplot(Season_Dino_BM, aes(x = Day, y = true_count)) +
  # Dinophysis GAM +- confidence interval
  geom_line(data = rpp_BM, aes(x = Day, y = fit), color = '#1E3523') +
  geom_ribbon(data = rpp_BM, aes(x = Day, y = fit,
                                 ymin = lwrS, ymax = uprS), 
              color = '#1E3523', fill = '#377185', alpha = 0.2) +
  # Dinophysis counts
  geom_point(color = '#1E3523', alpha = .8) +
  
  # Satellite observations of Mesodinium blooms
  geom_point(data = Satellite_survey_Meso, aes(x = Day, y = true_count),
             color = 'firebrick4', shape = 8, size = 3.5, stroke = .45) +
  
  # 1 subplot for each year
  facet_wrap(facets = c('Year'), nrow = 2) +
  # labels
  labs(y = 'Dinophysis cells observed in 10 mL',
       x = 'Julian day', title = 'Basse Michaud') +
  
  theme_classic()

# Seems quite nice. let's save this plot.
# ggsave('Basse_Michaud_multiyear.tiff', dpi = 300, height = 135, width = 164,
#                units = 'mm', compression = 'lzw')


## Cabourg ####

rpp_Cab <- filter(Dino_response_pred, Code_point_Libelle == 'Cabourg')
Season_Dino_Cab <- filter(Season_Dino, Code_point_Libelle == 'Cabourg')

# Let's plot that
ggplot(Season_Dino_Cab, aes(x = Day, y = true_count)) +
  # Dinophysis GAM +- confidence interval
  geom_line(data = rpp_Cab, aes(x = Day, y = fit), color = 'orangered') +
  geom_ribbon(data = rpp_Cab, aes(x = Day, y = fit,
                                 ymin = lwrS, ymax = uprS), 
              color = 'orangered',
              fill = '#FF6448', alpha = 0.2) +
  # Dinophysis counts
  geom_point(color = 'orangered', alpha = .8) +
  
  # # Satellite observations of Mesodinium blooms (don't have it (yet))
  # geom_point(data = Satellite_survey_Meso, aes(x = Day, y = true_count),
  #            color = 'firebrick4', shape = 8, size = 3.5, stroke = .45) +
  
  # 1 subplot for each year
  facet_wrap(facets = c('Year'), scales = 'free_y') +
  # labels
  labs(y = 'Dinophysis cells observed in 10 mL',
       x = 'Julian day', title = 'Cabourg') +
  
  theme_classic()

# There may be an issue with the very high count in 2007.
# ggsave('Cabourg_multiyear.tiff', dpi = 300, height = 164, width = 164,
#                units = 'mm', compression = 'lzw')
