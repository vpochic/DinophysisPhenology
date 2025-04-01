###### Dinophysis phenology: visualisations on multiple years ###
## V. POCHIC
# 2025-03-31

# Additional visualisations using the results of GAM, REPHY data and satellite
# observations, to observe Dinophysis dynamics over several years

#### Packages and functions ####
library(tidyverse)
library(ggplot2)
library(ggnewscale)
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
                         'Arcachon - Bouée 7', 'Teychan bis')) %>%
# We see there are 2 pits for many sites: 1 corresponds to the pit between the 
# 2 peaks, and 1 to the "renewal" after the winter. We will select only the 
# former (Day > 70 and < 300)
  filter(Day > 70 & Day < 300) %>%
  # Reorder the site as a factor
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis'))

GAM_peak_select <- filter(GAM_peak, Code_point_Libelle %in% 
                           c('Antifer ponton pétrolier', 'Cabourg',
                             'Men er Roue', 'Ouest Loscolo',
                             'Le Cornard', 'Auger',
                             'Arcachon - Bouée 7', 'Teychan bis')) %>%
  
  # Reorder the site as a factor
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis'))

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
           # All other sites
           1))))))))))))) %>%
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
                                          'Sète mer', 'Diana centre'))

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
# ggsave('Maxima_plot2_8sites.tiff', dpi = 300, height = 125, width = 164,
#                units = 'mm', compression = 'lzw')


### Environmental data ####

## Here we will call a table of environmental data to plot a heatmap
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

# Now we try to underlie a heatmap of temperature/chl a seasonality 
# (2 versions of the plot)
ggplot(Maxima_Dino_2_stats_select) +
  geom_tile(data = Table_hydro_daily_select, 
            aes(x = Day, y = Code_point_Libelle, fill = TEMP.med), #fill = CHLOROA.med
            alpha = .8) +
  # color palette for fill
  # Temperature version
  scale_fill_distiller(palette = 'RdBu', direction = -1) +
  # Chl a version
  # scale_fill_cmocean(name = 'algae') +
  
  # Add labels here so the name of the 'fill' legend is correct
  # Labels (for chl a plot fill = 'Median [chl a] (mg/m3)')
  labs(x = 'Day of the year', y = NULL, fill = 'Median SST (°C)',
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
# ggsave('Dino_max_temp.tiff', dpi = 300, height = 180, width = 150,
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
