### GAM of Dinophysis variation rate as function of temperature ###
# V. POCHIC, 2025/11/18

# The idea here is to compute a GAM of Dinophysis variation rate (obtained from
# the derivative of the Dinophysis GAM), as a function of temperature, for each
# sampling site

#### Packages ####
library(tidyverse)
library(ggplot2)
library(mgcv)
library(gratia)
library(RColorBrewer)
library(ggnewscale)

#### Import data ####

### True Dinophysis sampling data
Season_Dino <- read.csv2('Data/REPHY_outputs/Season_Dino_20250604.csv', header = TRUE, 
                         fileEncoding = 'ISO-8859-1')


Season_Dino_12sites <- filter(Season_Dino, Code_point_Libelle %in% 
                                c(# Baie de Seine
                                  'Antifer ponton pétrolier', 'Cabourg',
                                  # Bretagne Sud
                                  'Men er Roue', 'Ouest Loscolo',
                                  # Pertuis charentais
                                  'Auger', 'Le Cornard',
                                  # Arcachon
                                  'Arcachon - Bouée 7', 'Teychan bis')) %>%
  # Get rid of the one case in which Heure is ""
  # (Antifer, 2020/05/27, 0 Dinophysis)
  filter(Heure != "") %>%
  # Get a date-time variable by pasting together date and time
  mutate(DateTime = paste0(Date, ' ', Heure)) %>%
  # Get the Target date in the right format for corresponding with the stratif
  # dataset
  mutate(Target_Date = ymd_hms(DateTime)) %>%
  # And get the Date variable in the right format
  mutate(Date = ymd(Date))

### Dinophysis variation rate
Dino_variation_rate <- read.csv2('Data/GAM_outputs/Dinophysis/All_sites/gam_Dino_multiyear_deriv_20251117.csv',
                                 fileEncoding = 'ISO-8859-1', header = TRUE) %>%
  mutate(Date = as_date(Date))

### Hydrological parameters form the REPHY dataset
Table_hydro <- read.csv2('Data/REPHY_outputs/Table1_hydro_select.csv', header = TRUE, 
                         fileEncoding = "ISO-8859-1")

# Select only stations of interest from 2007 --> 2022
Table_hydro_select <- filter(Table_hydro, Code_point_Libelle %in% 
                               c(# Baie de Seine
                                 'Antifer ponton pétrolier', 'Cabourg',
                                 # Bretagne Sud
                                 'Men er Roue', 'Ouest Loscolo',
                                 # Pertuis charentais
                                 'Auger', 'Le Cornard',
                                 # Arcachon
                                 'Arcachon - Bouée 7', 'Teychan bis')
) %>%
  # Only after 2007
  filter(Year >= 2007) %>%
  # Only events for which Temperature and Salinity were properly recorded
  filter(is.na(TEMP) == FALSE) %>%
  filter(is.na(SALI) == FALSE) %>%
  # Change to date format
  mutate(Date = ymd(Date)) %>%
  mutate(Month = month(Date)) %>%
  mutate(Day = as.numeric(yday(Date))) %>%
  # Create a 'fortnight' variable to match the sampling frequency
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2))

# Checking for NAs
# Temperature
isNA_TEMP <- filter(Table_hydro_select, is.na(TEMP) == TRUE)
# Clean
# Salinity
isNA_SALI <- filter(Table_hydro_select, is.na(SALI) == TRUE)
# Clean
# Chlorophyll a
isNA_CHLOROA <- filter(Table_hydro_select, is.na(CHLOROA) == TRUE)
# Not clean. 733 obs with CHLOROA = NA, from 2007 to 2022, in several sites

#### Join datasets ####

# Join Season_Dino and variation rate
Table_data_multiyear <- left_join(Season_Dino_12sites, Dino_variation_rate, 
                                  by = c('Code_point_Libelle', 'Date'), 
                                  suffix = c('','')) %>%
  # change the name of the derivative for consistency
  mutate(.derivative = deriv)

# And then hydrology data
Table_data_multiyear <- left_join(Table_data_multiyear, Table_hydro_select, 
                                  by = c('Year', 'Day', 'Code_point_Libelle'),
                                  suffix = c('','')) %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis')) %>%
  filter(is.na(TEMP) == FALSE)

## Nice :)

#### GAM formulation ####

gam_deriv <- gam(data = Table_data_multiyear, 
                # Only a spline for the temperature
                # We'll use a thin plate basis spline
                # 'k = -1' allows the model to fix the 'best' number of basic
                # functions (= knots)
                formula = deriv~s(TEMP, bs = 'tp', k = -1,
                                       # separate each site
                                       by = Code_point_Libelle),
                # Using a gaussian distribution for the derivative data
                family = gaussian(),
                # Restricted maximum likelihood estimation (recommended method)
                method = 'REML')

summary(gam_deriv)
gam.check(gam_deriv)

#### Model checks ####

# Create a new 'empty' dataset for storing model prediction (fit)
gam_deriv_newdata <- expand_grid(TEMP=seq(min(Table_data_multiyear$TEMP), max(Table_data_multiyear$TEMP),
                                          # We will predict on a range of temperature that is precise to
                                          # 0.1°C
                                          by = 0.1),
                                # We also add the site data
                                Code_point_Libelle = unique(Table_data_multiyear$Code_point_Libelle))

## Get the inverse link function of the model
# With this function, we can transform the prediction we make on the link
#scale to the response scale
ilink <- gam_deriv$family$linkinv

## Predict : add fit and se.fit on the **link** scale
gam_deriv_newdata <- bind_cols(gam_deriv_newdata, 
                              setNames(as_tibble(
                                predict(gam_deriv, gam_deriv_newdata, 
                                        se.fit = TRUE, type = 'link',
                                        re.form = ~ 1|Year)[1:2]),
                                nm = c('fit_link','se_link')))

## Create the 95% confidence interval (2*standard error fit) and backtransform 
# to response variable using the inverse link function
gam_deriv_newdata <- mutate(gam_deriv_newdata,
                           fit_resp  = ilink(fit_link),
                           right_upr = ilink(fit_link + (2 * se_link)),
                           right_lwr = ilink(fit_link - (2 * se_link)))
# Check the confidence interval. It should not extend below 0 (negative counts
# are impossible)
min(gam_deriv_newdata$right_lwr) # Nice :)
min(gam_deriv_newdata$fit_resp)
max(gam_deriv_newdata$right_upr) # Nice too
max(gam_deriv_newdata$fit_resp)

### Plotting ####

# Color palette
pheno_palette8 <- c('red3', 'orangered', '#2156A1', '#5995E3', 
                     '#1F3700', '#649003','#F7B41D', '#FBB646')

ggplot(gam_deriv_newdata, aes(x = TEMP, y = fit_resp, color = Code_point_Libelle)) +
  geom_ribbon(aes(x=TEMP, ymin=right_lwr, ymax=right_upr,
                  color = Code_point_Libelle), fill = 'grey70', alpha=0.7) +
  geom_line(linewidth = 1) +
  geom_point(data = Table_data_multiyear,
             aes(x = TEMP, y = deriv, color = Code_point_Libelle), shape = 21) +
  scale_color_discrete(type = pheno_palette8) +
  facet_wrap(facets = 'Code_point_Libelle', scale = 'free_y', nrow = 2) +
  theme_classic()

# Saving plot
# ggsave('Plots/GAMs/Dinophysis/gam_deriv_8_sites.tiff', dpi = 300, height = 120, width = 160, 
#        units = 'mm', compression = 'lzw')

## Bonus : rug plot of temperature ####

# Create a separate data table for temperature
temp_data_plot <- gam_deriv_newdata %>%
  # We select 1 sampling point at random (all of them have the whole temperature
  # range in the 'newdata' table)
  filter(Code_point_Libelle == 'Auger') %>%
  select(c('TEMP'))

ggplot(gam_deriv_newdata) +
  # rug plot for temperature range (purely visual)
  geom_rug(data = temp_data_plot, aes(x = TEMP, color = TEMP),
           linewidth = .3,
           length = unit(0.5, 'cm')) +
  scale_color_distiller(palette = 'RdBu', direction = -1, guide = 'none') +
  new_scale_color() +
  # Points for true data
  geom_point(data = Table_data_multiyear,
             aes(x = TEMP, y = deriv, color = Code_point_Libelle), shape = 21,
             alpha = .3) +
  # Ribbon for uncertainty
  geom_ribbon(aes(x=TEMP, ymin=right_lwr, ymax=right_upr,
                  color = Code_point_Libelle, fill = Code_point_Libelle),
              alpha=0.2) +
  # Line for model fit
  geom_line(aes(x = TEMP, y = fit_resp, color = Code_point_Libelle),
            linewidth = 1) +
  # color scales
  scale_color_discrete(type = pheno_palette8, guide = 'none') +
  scale_fill_discrete(type = pheno_palette8, guide = 'none') +
  facet_wrap(facets = 'Code_point_Libelle', scale = 'free_y', nrow = 2) +
  # labels
  labs(x = 'Temperature (°C)', y = 'Dinophysis variation rate (cells/10mL/d-1)') +
  theme_classic() +
  theme(strip.text.x = element_text(size = 7.2, colour = "black"))

# Save plot
# ggsave('Plots/GAMs/Dinophysis/gam_deriv_8_sites_rug.tiff', dpi = 300, height = 120, width = 164,
#        units = 'mm', compression = 'lzw')
