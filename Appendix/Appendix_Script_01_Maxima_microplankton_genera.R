#### Additional script: broad phenology of other microplankton genera ###
# Author: V. POCHIC
# Last modif: 2026/07/06

## General information ####

## Description: 
# The goal of this script is to check if the regularity in the timing of the
# yearly maximum is a Dinophysis-only thing or if it's broadly shared among
# other microplankton genera (Thomas' idea)
# We will do that for 4 genera: Dinophysis (to compare with the others),
# Protoperidinium (dino, heterotrophic), Alexandrium (dino, photo/mixotrophic),
# Pseudo-nitzschia (diatom, autotrophic)
# We'll do this for the 4 sampling sites where Dinophysis regularity in yearly
# maximum is the most pronounced, i.e. in Seine Bay and Southern Brittany.

## Files required:
# Data/REPHY_outputs/Table_phyto_taxon.csv
## --> REPHY data to go look into for other phytoplankton genera
# Data/GAM_outputs/Dinophysis/All_sites/Dino_GAM_response_pred.csv
## --> GAM fit of Dinophysis phenology (function of DOY)
# Data/Models/ERA5/Outputs/...
# Data/Models/GAMAR/Outputs/Stratif_index_GAMAR_12sites.csv
## --> Data from oceanographic/atmospheric physical models

## Outputs:
## Figures:
# Figure A1 - Distribution of occurrences for different phyto genera
# Figure A2 - Timing of Dinophysis (+ other genera) maxima over the study 
# period (+ GLM fit)

#### Required packages ####

library(tidyverse)
library(tidymodels)
library(mgcv)
library(ggplot2)
library(ggpubr)
library(grid)
library(deeptime)

#### Import data ####

Table_phyto_taxon <- read.csv2('Data/REPHY_outputs/Table_phyto_taxon.csv',
                        header = TRUE, fileEncoding = 'ISO-8859-1')

# Select only sampling sites of interest from 2007 --> 2022
Table_phyto_select <- filter(Table_phyto_taxon, Code_point_Libelle %in% 
                               c(# Baie de Seine
                                 'Antifer ponton pétrolier', 'Cabourg',
                                 # Bretagne Sud
                                 'Men er Roue', 'Ouest Loscolo')
) %>%
  # Only years in the period [2007-2022]
  filter(Year >= 2007 & Year <= 2022)

## More zeros

### The goal here is to create a table where all absences of detection are 
# listed as zeros
Table_phyto_zeros <- Table_phyto_select %>%
  # 'Comptage' is the cell count
  pivot_wider(names_from = Taxon, values_from = Comptage, values_fill = 0) %>%
  # Create Latitude and Longitude variables for selection (because lon and lat
  # introduce a bunch of unwanted taxa )
  mutate(Longitude = as.numeric(lon)) %>%
  mutate(Latitude = as.numeric(lat)) %>%
  # We keep only desired genera and the metadata columns
  select(starts_with('Dinophysis') | starts_with('Pseudo-nitzschia') |
           starts_with('Protoperidinium') | starts_with('Alexandrium') |
           contains(c('Code_point_Libelle', 'Year', 'Month',
                      'Date', 'Heure', 'Code.parametre', 
                      'Longitude', 'Latitude', 
                      'ID.interne.passage'))) %>%
  # create a day of the year variable and get the Date in a 'Date' format
  mutate(Date = ymd(Date)) %>%
  mutate(Day = as.numeric(yday(Date)))

#### Compute smaller, taxon-specific datasets ####

# Dinophysis
Season_Dino <- Table_phyto_zeros %>%
  # We create a count variable for Dinophysis as a genus
  mutate(Dinophysis_genus = rowSums(across(contains('Dinophysis')))) %>%
  #--- Create a 'true count' variable ---#
  # this variable corresponds to the number of cells that were actually
  # counted by the operator (in 10 mL). This variable will follow a Poisson
  # distribution, contrary to 'Dinophysis_genus' because the conversion from
  # 10 mL to 1L (*100) in the original data prevents some intermediate values 
  # (e.g., 150 cells.L-1)
  mutate(true_count = Dinophysis_genus/100) %>%
  # Filter out exceptionnaly high counts (>500 cells observed in 10 mL)
  # this represents 3 events (2 in Antifer and 1 in Cabourg)
  filter(true_count < 500)

# Protoperidinium
Season_Protop <- Table_phyto_zeros %>%
  # We create a count variable for Protoperidinium as a genus
  mutate(Protop_genus = rowSums(across(contains('Protoperidinium')))) %>%
  #--- Create a 'true count' variable ---#
  # this variable corresponds to the number of cells that were actually
  # counted by the operator (in 10 mL). This variable will follow a Poisson
  # distribution, contrary to 'Dinophysis_genus' because the conversion from
  # 10 mL to 1L (*100) in the original data prevents some intermediate values 
  # (e.g., 150 cells.L-1)
  mutate(true_count = Protop_genus/100)

# Alexandrium
Season_Alex <- Table_phyto_zeros %>%
  # We create a count variable for Alexandrium as a genus
  mutate(Alex_genus = rowSums(across(contains('Alexandrium')))) %>%
  #--- Create a 'true count' variable ---#
  # this variable corresponds to the number of cells that were actually
  # counted by the operator (in 10 mL). This variable will follow a Poisson
  # distribution, contrary to 'Dinophysis_genus' because the conversion from
  # 10 mL to 1L (*100) in the original data prevents some intermediate values 
  # (e.g., 150 cells.L-1)
  mutate(true_count = Alex_genus/100)

# Pseudo-nitzschia
Season_PN <- Table_phyto_zeros %>%
  # We create a count variable for Pseudo-nitzschia as a genus
  mutate(PN_genus = rowSums(across(contains('Pseudo-nitzschia')))) %>%
  #--- Create a 'true count' variable ---#
  # this variable corresponds to the number of cells that were actually
  # counted by the operator (in 10 mL). This variable will follow a Poisson
  # distribution, contrary to 'Dinophysis_genus' because the conversion from
  # 10 mL to 1L (*100) in the original data prevents some intermediate values 
  # (e.g., 150 cells.L-1)
  mutate(true_count = PN_genus/100)

#---- Plots: barplot of presence (seasonality) ----####

# For each of our 4 genera, we will make a barplot of how much each taxon is
# present (= detected) by fortnight. That way, we can assess if the distribution
# is unimodal or bimodal (or else), and inform how we analyse the timing of
# yearly maximum.

# We will keep the information of the different taxa present, to give us an
# idea of the variability in seasonality across species within a genus.

# (Each plot script is a bit long)

## Dinophysis ####

Season_Dino_longer <- Season_Dino %>%
  # Pivot longer with values to 'Comptage' (cell count)
  pivot_longer(cols = contains(c('Dinophysis')), 
               names_to = 'Taxon', 
               values_to = 'Comptage') %>%
  # Remove Dinophysis_genus (else we count everything twice)
  filter(Taxon != 'Dinophysis_genus') %>%
  # Remove these values that were calculated on Dinophysis_genus
  select(-c('true_count')) %>%
  # Create a 'Fortnight' variable as sampling occurs fortnightly
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2))

# A plot-ready version of the table
Season_Dino_fortnight <- Season_Dino_longer %>%
  # Group data by fortnight to get a 1-year seasonality
  group_by(Code_point_Libelle, Fortnight, Taxon) %>%
  # Summarise by fortnight
  summarise(Comptage = sum(Comptage), .groups = 'keep')


#--- Now plot ---#

## Color palette (for eight taxa)
taxo_palette8 <- c(# Dinophysis  
  '#11203E',
  # Dinophysis acuminata-complex (acuminata, sacculus, fortii)
  # Though I'm not really sure anymore than D. fortii is part of the complex...
  '#435E7B', '#7F96B6', '#BBD4F2',
  # Dinophysis acuta
  '#8D6456',
  # Dinophysis caudata and tripos
  '#FBA823', '#FBB646',
  # Dinophysis + phalacroma
  '#FF6448'
)

# Changing factor order so sampling sites appear in the order we want
Season_Dino_fortnight <- Season_Dino_fortnight %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# And we do the same for the Dinophysis taxa
Season_Dino_fortnight <- Season_Dino_fortnight %>%
  mutate(Taxon = as_factor(Taxon)) %>%
  mutate(Taxon = fct_relevel(Taxon,
                             'Dinophysis', 'Dinophysis acuminata',
                             'Dinophysis sacculus', 'Dinophysis fortii',
                             'Dinophysis acuta', 'Dinophysis caudata',
                             'Dinophysis tripos',
                             'Dinophysis + phalacroma')) %>%
  # Add a 'true count' variable that identifies the number of cells actually 
  # observed and counted
  mutate(true_count = Comptage/100)

# Little fancy trick: we want to color the top labels of each subplot by
# phenoregion. We will create a dataframe with the names of the sampling sites
# that make our facets, and the colors we want associated to them.
# The name of the columns in the dataframe is super important here!
# (i.e., you have to name them 'name' and 'color')
df_color <- data.frame(
  name = c('Antifer ponton pétrolier', 'Cabourg',
           'Men er Roue', 'Ouest Loscolo'),
  color = c('red3', 'orangered', '#2156A1', '#5995E3'))

## Plot command
Dino_barplot <- ggplot(Season_Dino_fortnight, aes(x = Fortnight, y = true_count, fill = Taxon)) +
  geom_col() +
  # Facet wrap
  facet_wrap_color(facets = c('Code_point_Libelle'), scales = 'free_y',
                   colors = df_color, nrow = 1) +
  # axis
  scale_x_continuous(breaks = c(1, 8, 16, 20, 26),
                     labels = c('1', '8', '16', '20', '26')) +
  # Labels
  labs(title = "Dinophysis (2007-2022)",
       y = "Sum of observed Dinophysis cells",
       x = 'Fortnight',
       fill = NULL
  ) +
  # Legend labels
  scale_fill_discrete(# Color palette
                      type = taxo_palette8,
                      # guide legend
                      guide = 'none') +
  # Theme
  theme_classic() +
  theme(# title and facets text
    title = element_text(size = 11),
    strip.text = element_text(size = 7.5),
    # axes text
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7.5),
    # legend
    legend.key.height = unit(.75, 'cm'),
    legend.key.width = unit(.6, 'cm'),
    legend.title = element_text(size=9),
    legend.text = element_text(size=8),
    legend.position = 'bottom',
    legend.background = element_rect(color = 'black', linewidth = .25),
    strip.background = element_rect(color = 'white')
  )

# Ok great we have our first plot, which is just an adaptation of Figure S2.
# Now let's see how it goes with the 3 other genera.
## Pseudo-nitzschia ####

Season_PN_longer <- Season_PN %>%
  # Pivot longer with values to 'Comptage' (cell count)
  pivot_longer(cols = contains(c('Pseudo-nitzschia')), 
               names_to = 'Taxon', 
               values_to = 'Comptage') %>%
  # Remove PN_genus (else we count everything twice)
  filter(Taxon != 'PN_genus') %>%
  # Remove these values that were calculated on Pseudo-nitzschia_genus
  select(-c('true_count')) %>%
  # Create a 'Fortnight' variable as sampling occurs fortnightly
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2))

# A plot-ready version of the table
Season_PN_fortnight <- Season_PN_longer %>%
  # Group data by fortnight to get a 1-year seasonality
  group_by(Code_point_Libelle, Fortnight, Taxon) %>%
  # Summarise by fortnight
  summarise(Comptage = sum(Comptage), .groups = 'keep')


#--- Now plot ---#

## Color palette (for 10 taxa)
taxo_palette10 <- c(# Pseudo-nitzschia  
  '#9DA51E', 
  # Pseudo-nitzschia americana-complex
  '#B47E24', '#F7B41D',
  # Pseudo-nitzschia seriata-complex
  '#CCA377', '#869B7B', '#041503', 
  # Pseudo-nitzschia delicatissima-complex
  '#649003',
  # Pseudo-nitzschia large-group
  '#879BBD', '#C2D4E9',
  # Pseudo-nitzschia sigmoid-group
  '#292F06'
)

# Changing factor order so sampling sites appear in the order we want
Season_PN_fortnight <- Season_PN_fortnight %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# And we do the same for the Pseudo-nitzschia taxa
Season_PN_fortnight <- Season_PN_fortnight %>%
  mutate(Taxon = as_factor(Taxon)) %>%
  mutate(Taxon = fct_relevel(Taxon,
                             'Pseudo-nitzschia', 'Pseudo-nitzschia americana',
                             'Pseudo-nitzschia, complexe americana (americana + brasiliana)',
                             'Pseudo-nitzschia fraudulenta', 
                             'Pseudo-nitzschia, complexe seriata, groupe des effilées (multiseries + pungens)',
                             'Pseudo-nitzschia, complexe seriata, groupe des larges (australis + fraudulenta + seriata + subpacifica)',
                             'Pseudo-nitzschia, complexe delicatissima, groupe des fines (calliantha + delicatissima + pseudodelicatissima + subcurvata)', 
                             'Pseudo-nitzschia, groupe des larges asymétriques (australis + seriata + subpacifica)',
                             'Pseudo-nitzschia, groupe des larges symétriques (fraudulenta)',
                             'Pseudo-nitzschia, groupe des sigmoïdes (multistriata)')) %>%
  # Add a 'true count' variable that identifies the number of cells actually 
  # observed and counted
  mutate(true_count = Comptage/100)

# Little fancy trick: we want to color the top labels of each subplot by
# phenoregion. We will create a dataframe with the names of the sampling sites
# that make our facets, and the colors we want associated to them.
# The name of the columns in the dataframe is super important here!
# (i.e., you have to name them 'name' and 'color')
df_color <- data.frame(
  name = c('Antifer ponton pétrolier', 'Cabourg',
           'Men er Roue', 'Ouest Loscolo'),
  color = c('red3', 'orangered', '#2156A1', '#5995E3'))

## Plot command
PN_barplot <- ggplot(Season_PN_fortnight, aes(x = Fortnight, y = true_count, fill = Taxon)) +
  geom_col() +
  # Facet wrap
  facet_wrap_color(facets = c('Code_point_Libelle'), scales = 'free_y',
                   colors = df_color, nrow = 1) +
  # axis
  scale_x_continuous(breaks = c(1, 8, 16, 20, 26),
                     labels = c('1', '8', '16', '20', '26')) +
  # Labels
  labs(title = "Pseudo-nitzschia (2007-2022)",
       y = "Sum of observed Pseudo-nitzschia cells",
       x = 'Fortnight',
       fill = NULL
  ) +
  # Legend labels
  scale_fill_discrete(# Color palette
    type = taxo_palette10,
    # guide legend
    guide = 'none') +
  # Theme
  theme_classic() +
  theme(# title and facets text
    title = element_text(size = 11),
    strip.text = element_text(size = 7.5),
    # axes text
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7.5),
    # legend
    legend.key.height = unit(.75, 'cm'),
    legend.key.width = unit(.6, 'cm'),
    legend.title = element_text(size=9),
    legend.text = element_text(size=8),
    legend.position = 'bottom',
    legend.background = element_rect(color = 'black', linewidth = .25),
    strip.background = element_rect(color = 'white')
  )

# Nice. Let's see the next one.
## Alexandrium ####

Season_Alex_longer <- Season_Alex %>%
  # Pivot longer with values to 'Comptage' (cell count)
  pivot_longer(cols = contains(c('Alexandrium')), 
               names_to = 'Taxon', 
               values_to = 'Comptage') %>%
  # Remove Alex_genus (else we count everything twice)
  filter(Taxon != 'Alex_genus') %>%
  # Remove these values that were calculated on Alex_genus
  select(-c('true_count')) %>%
  # Create a 'Fortnight' variable as sampling occurs fortnightly
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2))

# A plot-ready version of the table
Season_Alex_fortnight <- Season_Alex_longer %>%
  # Group data by fortnight to get a 1-year seasonality
  group_by(Code_point_Libelle, Fortnight, Taxon) %>%
  # Summarise by fortnight
  summarise(Comptage = sum(Comptage), .groups = 'keep')


#--- Now plot ---#

## Color palette (for eight taxa)
taxo_palette6 <- c(# Alexandrium  
  '#711412',
  # Alexandrium affine
  '#8B064B',
  # Alexandrium minutum
  '#DEB1CC',
  # Alexandrium ostenfeldii
  '#FC4D6B',
  # Alexandrium pseudogonyaulax
  '#FF6448',
  # Alexandrium tamutum
  '#642C3A'
)

# Changing factor order so sampling sites appear in the order we want
Season_Alex_fortnight <- Season_Alex_fortnight %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# And we do the same for the Dinophysis taxa
Season_Alex_fortnight <- Season_Alex_fortnight %>%
  mutate(Taxon = as_factor(Taxon)) %>%
  mutate(Taxon = fct_relevel(Taxon,
                             'Alexandrium', 'Alexandrium affine',
                             'Alexandrium minutum',
                             'Alexandrium ostenfeldii',
                             'Alexandrium pseudogonyaulax',
                             'Alexandrium tamutum')) %>%
  # Add a 'true count' variable that identifies the number of cells actually 
  # observed and counted
  mutate(true_count = Comptage/100)

# Little fancy trick: we want to color the top labels of each subplot by
# phenoregion. We will create a dataframe with the names of the sampling sites
# that make our facets, and the colors we want associated to them.
# The name of the columns in the dataframe is super important here!
# (i.e., you have to name them 'name' and 'color')
df_color <- data.frame(
  name = c('Antifer ponton pétrolier', 'Cabourg',
           'Men er Roue', 'Ouest Loscolo'),
  color = c('red3', 'orangered', '#2156A1', '#5995E3'))

## Plot command
Alex_barplot <- ggplot(Season_Alex_fortnight, aes(x = Fortnight, y = true_count, fill = Taxon)) +
  geom_col() +
  # Facet wrap
  facet_wrap_color(facets = c('Code_point_Libelle'), scales = 'free_y',
                   colors = df_color, nrow = 1) +
  # axis
  scale_x_continuous(breaks = c(1, 8, 16, 20, 26),
                     labels = c('1', '8', '16', '20', '26')) +
  # Labels
  labs(title = "Alexandrium (2007-2022)",
       y = "Sum of observed Alexandrium cells",
       x = 'Fortnight',
       fill = NULL
  ) +
  # Legend labels
  scale_fill_discrete(# Color palette
    type = taxo_palette6,
    # guide legend
    guide = 'none') +
  # Theme
  theme_classic() +
  theme(# title and facets text
    title = element_text(size = 8),
    strip.text = element_text(size = 7.5),
    # axes text
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7.5),
    # legend
    legend.key.height = unit(.75, 'cm'),
    legend.key.width = unit(.6, 'cm'),
    legend.title = element_text(size=9),
    legend.text = element_text(size=8),
    legend.position = 'bottom',
    legend.background = element_rect(color = 'black', linewidth = .25),
    strip.background = element_rect(color = 'white')
  )

# Ok. Next one.


## Protoperidinium ####

Season_Protop_longer <- Season_Protop %>%
  # Pivot longer with values to 'Comptage' (cell count)
  pivot_longer(cols = contains(c('Protoperidinium')), 
               names_to = 'Taxon', 
               values_to = 'Comptage') %>%
  # Remove Protop_genus (else we count everything twice)
  filter(Taxon != 'Protop_genus') %>%
  # Remove these values that were calculated on Protop_genus
  select(-c('true_count')) %>%
  # Create a 'Fortnight' variable as sampling occurs fortnightly
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2))

# A plot-ready version of the table
Season_Protop_fortnight <- Season_Protop_longer %>%
  # Group data by fortnight to get a 1-year seasonality
  group_by(Code_point_Libelle, Fortnight, Taxon) %>%
  # Summarise by fortnight
  summarise(Comptage = sum(Comptage), .groups = 'keep')


#--- Now plot ---#

## Color palette (for eight taxa)
taxo_palette8b <- c(# Protoperidinium  
  'sienna',
  # Protoperidinium + Peridinium
  'tan3',
  # Protoperidinium bipes
  'sienna1',
  # Protoperidinium conicum
  'tomato2',
  # Protoperidinium depressum
  'red2',
  # Protoperidinium divergens
  'red4',
  # Protoperidinium steinii
  'orange',
  # Protoperidinium steinii + pyriforme
  'khaki3'
)

# Changing factor order so sampling sites appear in the order we want
Season_Protop_fortnight <- Season_Protop_fortnight %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# And we do the same for the Dinophysis taxa
Season_Protop_fortnight <- Season_Protop_fortnight %>%
  mutate(Taxon = as_factor(Taxon)) %>%
  mutate(Taxon = fct_relevel(Taxon,
                             'Protoperidinium', 'Protoperidinium + Peridinium',
                             'Protoperidinium bipes', 'Protoperidinium conicum',
                             'Protoperidinium depressum', 'Protoperidinium divergens',
                             'Protoperidinium steinii',
                             'Protoperidinium steinii + pyriforme')) %>%
  # Add a 'true count' variable that identifies the number of cells actually 
  # observed and counted
  mutate(true_count = Comptage/100)

# Little fancy trick: we want to color the top labels of each subplot by
# phenoregion. We will create a dataframe with the names of the sampling sites
# that make our facets, and the colors we want associated to them.
# The name of the columns in the dataframe is super important here!
# (i.e., you have to name them 'name' and 'color')
df_color <- data.frame(
  name = c('Antifer ponton pétrolier', 'Cabourg',
           'Men er Roue', 'Ouest Loscolo'),
  color = c('red3', 'orangered', '#2156A1', '#5995E3'))

## Plot command
Protop_barplot <- ggplot(Season_Protop_fortnight, aes(x = Fortnight, y = true_count, fill = Taxon)) +
  geom_col() +
  # Facet wrap
  facet_wrap_color(facets = c('Code_point_Libelle'), scales = 'free_y',
                   colors = df_color, nrow = 1) +
  # axis
  scale_x_continuous(breaks = c(1, 8, 16, 20, 26),
                     labels = c('1', '8', '16', '20', '26')) +
  # Labels
  labs(title = "Protoperidinium (2007-2022)",
       y = "Sum of observed Protoperidinium cells",
       x = 'Fortnight',
       fill = NULL
  ) +
  # Legend labels
  scale_fill_discrete(# Color palette
    type = taxo_palette8b,
    # guide legend
    guide = 'none') +
  # Theme
  theme_classic() +
  theme(# title and facets text
    title = element_text(size = 8),
    strip.text = element_text(size = 7.5),
    # axes text
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7.5),
    # legend
    legend.key.height = unit(.75, 'cm'),
    legend.key.width = unit(.6, 'cm'),
    legend.title = element_text(size=9),
    legend.text = element_text(size=8),
    legend.position = 'bottom',
    legend.background = element_rect(color = 'black', linewidth = .25),
    strip.background = element_rect(color = 'white')
  )

# Ok, done!


### Arrange all 4 plots ####

# Now let's arrange all 4 plots so they fit all in 1 figure.
barplots4 <- ggarrange(Dino_barplot + rremove("ylab") + rremove("xlab"),
                       Alex_barplot + rremove("ylab") + rremove("xlab"),
                       Protop_barplot + rremove("ylab") + rremove("xlab"),
                       PN_barplot + rremove("ylab") + rremove("xlab"),
                       nrow = 4, align = 'v')

annotate_figure(barplots4, 
               left = textGrob('Total number of cells observed', 
                               rot = 90, vjust = 1,
                               gp=gpar(fontsize=10, col="black")),
               bottom = textGrob('Fortnight of sampling',
                                 gp=gpar(fontsize=10, col="black")))

# Save the plot
# ggsave('Appendix/Plots/Microplankton_genera/FigA01_Barplots_4genera.tiff',
#        height = 180, width = 164, units = 'mm', compression = 'lzw',
#        dpi = 500)

# Ok it's quite nice! Now let's do the hard part: the yearly maxima.

####-----------------------------------------------------------------------####
#### Yearly maxima over the study period ####

# The objective is to look at the yearly maximum of each genus over the study
# period, and see how it compares to Dinophysis in terms of regularity.
# (and, as a bonus, to see if there's a trend)

# For Dinophysis in the paper, we used the GAM to determine if we needed to
# split the year into 2 periods. We can't do that because we won't make a GAM
# for each genus, so instead we'll decide based on the looks of the distribution
# we just plotted: unimodal = 1 period, bimodal = 2 periods.

### Dinophysis ####

## Yearly maxima ####

# Let's compute the yearly maximum for each year of the study period

# First, get rid of zeros
Season_Dino_nozeros <- filter(Season_Dino, true_count != 0)

# Then, establish periods. Based on our plots of the distribution of cells
# observed, we have a unimodal distribution in Antifer and Cabourg (i.e., 1
# period) and a bimodal distribution in Men er Roue and Ouest Loscolo (i.e., 
# 2 periods)

Maxima_Dino <- Season_Dino_nozeros %>%
  # The first mode of the distribution in MeR and OL ends around fortnight 16,
  # so let's cut the year like that:
  mutate(period = ifelse(Code_point_Libelle %in% c('Men er Roue', 'Ouest Loscolo')
                         & Day > (16*14), 2, 1))

# Let's compute the number of observations of Dinophysis by period
Maxima_Dino_n <- Maxima_Dino %>%
  # groups
  group_by(Code_point_Libelle, Year, period) %>%
  # summarise to get n
  summarise(n_max = n(), max_Dino = max(true_count), .groups = 'keep') %>%
  # Create a "shape" variable for plotting
  mutate(shape = ifelse(period == 1 & n_max > 1 & max_Dino > 1, 1, 
                        ifelse(period == 2 & n_max > 1 & max_Dino > 1, 2,
                               ifelse(period == 1 & (n_max == 1 || max_Dino == 1), 3,
                                      4))))

# Now extract maxima. Careful, this step has to be done AFTER the previous one.
Maxima_Dino <- Maxima_Dino %>%
  # Add period as a grouping variable
  group_by(Code_point_Libelle, Year, period) %>%
  # extract rows with yearly maxima of 'true_count' (slice() retains groups)
  slice(which.max(true_count)) %>%
  # Keep only Dinophysis cell densities >= 100 cells/L (more than 1 cell counted)
  filter(true_count >= 1)

# Join both datasets to get n() info on maxima
Maxima_Dino <- left_join(Maxima_Dino, Maxima_Dino_n) %>%
  # Code_point_Libelle as factor
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# Seems perfect.

# Plot this

# The df color for facet labels
df_color <- data.frame(
  name = c('Antifer ponton pétrolier', 'Cabourg',
           'Men er Roue', 'Ouest Loscolo'),
  color = c('red3', 'orangered', '#2156A1', '#5995E3'))

# plot
ggplot(Maxima_Dino, 
       aes(x = Year, y = Day, 
           shape = as.factor(shape))) +
  # then the points
  geom_point(size = 4, alpha = .8,
             color = '#11203E') +
  # aesthetics
  facet_wrap_color(facets = 'Code_point_Libelle', nrow = 1,
                   colors = df_color, axes = 'all_y',
                   axis.labels = 'margins') +
  scale_y_continuous(limits = c(1, 365), 
                     breaks = c(1, 100, 200, 300, 365)) +
  # shape scale
  scale_shape_manual(values = c(16, 17, 1, 2), guide = 'none') +
  # scale_shape_manual(values = c(1, 3), guide = 'none') +
  # labels
  labs(y = 'Day of maximum Dinophysis count') +
  theme_classic()

# Ok that's nice. Now let's model that properly with some good ol' GLMs

## GLMs ####

# Very simple: our response variable is day of Dino max (DOY), our predictor is
# the Year

# We'll do GLMs, with the help of the tutorial by Dr. Bede Davies
# https://bedeffinianrowedavies.com/statisticstutorials/poissonglms

Maxima_Dino_trends <- Maxima_Dino %>%
  # We create a category variable that is the sampling site + the period
  mutate(category = paste(Code_point_Libelle, period, sep = "_")) %>%
  # Relevel as a factor (category = Sampling site + period)
  mutate(category = as_factor(category)) %>%
  mutate(category = fct_relevel(category,
                                'Antifer ponton pétrolier_1', 'Cabourg_1',
                                'Men er Roue_1', 'Men er Roue_2',
                                'Ouest Loscolo_1', 'Ouest Loscolo_2'))

# Let's transform our response variable so it follows a beta distribution
# (bounded between 0 and 1). This allows us to fit a GLM properly
Maxima_Dino_trends <- Maxima_Dino_trends %>%
  # For this we divide Day by the maximum value of this variable (365)
  mutate(Day_beta = Day/365)

# Now, we'll fit a GLM for every series of maxima (1 for each period in each 
# sampling site)
# With the beta-transformed variable (we'll use mgcv's gam() function)
glm_trends_beta <- gam(Day_beta~Year*category,
                       data = Maxima_Dino_trends,
                       # Specifying a beta distribution of the variable
                       family = betar(link="logit"))

summary(glm_trends_beta)

### Diagnostic plots ###

### Checking the model
ModelOutputs<-data.frame(Fitted=fitted(glm_trends_beta),
                         Residuals=resid(glm_trends_beta))

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- glm_trends_beta$model
# then we add the values of fitted and residuals
qq_data <- bind_cols(qq_data, ModelOutputs)

# Reorder data
qq_data_reordered <- qq_data %>%
  mutate(category = fct_relevel(category,
                                'Antifer ponton pétrolier_1', 'Cabourg_1',
                                'Men er Roue_1', 'Men er Roue_2',
                                'Ouest Loscolo_1', 'Ouest Loscolo_2')) %>%
  # we ungroup the data frame to produce only 1 qq-plot
  ungroup()

# And (qq-)plot
qqplot_custom <- ggplot(qq_data_reordered) +
  stat_qq(aes(sample=Residuals, color = category), alpha = .7) +
  stat_qq_line(aes(sample=Residuals, color = category)) +
  facet_wrap(facets = c('category')) +
  scale_color_discrete(guide = 'none') +
  theme_classic() +
  labs(y="Sample Quantiles",x="Theoretical Quantiles")

qqplot_custom

# Ok!

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data_reordered)+
  geom_point(aes(x=Fitted,y=Residuals, color = category), 
             alpha = .7) +
  facet_wrap(facets = c('category'), scales = 'free') +
  scale_color_discrete(guide = 'none') +
  theme_classic() +
  labs(y="Residuals",x="Fitted Values")

RvFplot_custom

# Okish

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data_reordered, aes(x = Residuals, 
                                                fill = category))+
  geom_histogram(binwidth = 1)+
  facet_wrap(facets = c('category'), scales = 'free') +
  scale_fill_discrete(guide = 'none') +
  theme_classic() +
  labs(x='Residuals', y = 'Count')

HistRes_custom

# It's approximately ok

### Let's predict with the GLM
NewData_1 <- expand_grid(# Year variable
  Year=seq(min(Maxima_Dino_trends$Year), 
           max(Maxima_Dino_trends$Year)),
  # We also add the site data
  category = unique(Maxima_Dino_trends$category))

Pred <- predict(glm_trends_beta, NewData_1, se.fit=TRUE, type="response")

NewData<-NewData_1 %>% 
  mutate(response=Pred$fit,
         se.fit=Pred$se.fit,
         Upr=response+(se.fit*1.96),
         Lwr=response-(se.fit*1.96))

# Join datasets of prediction and true data
Data_trends_Dino <- left_join(Maxima_Dino_trends, NewData)

#---- Plot with the GLM fits ----####
# Plot (we multiply the predictions made on the beta-transformed variable
# by 365 so that it fits the true data)
plot_trends_Dino <- ggplot(Data_trends_Dino)+
  # confidence interval of model
  geom_ribbon(aes(x = Year,
                  ymax = Upr*365,
                  ymin = Lwr*365, 
                  group = category),
              alpha=0.35, fill = '#11203E') +
  # fit of model
  geom_line(aes(x = Year,
                y = response*365,
                group = category),
            color = '#11203E') +
  # true data
  geom_point(aes(x = Year, y = Day,
                 shape = as_factor(shape)), color = '#11203E') +
  # shape scale
  scale_shape_manual(values = c(16, 17, 1, 2), guide = 'none') +
  # y-axis scale
  scale_y_continuous(limits = c(1, 410), 
                     breaks = c(1, 100, 200, 300, 365)) +
  # Labels
  labs(title = "Dinophysis (2007-2022)",
       x = "Year", 
       y = "Day of maximum Dinophysis count") +
  # Facet wrap
  facet_wrap_color(facets = c('Code_point_Libelle'),
                   colors = df_color, nrow = 1, axes = 'all_y', 
                   axis.labels = 'margins') +
  theme_classic() +
  theme(title = element_text(size = 8),
    strip.text = element_text(
    size = 7, color = 'black'),
    strip.background = element_rect(colour = 'white'))

# Excellent stuff. Now to the other genera.

### Alexandrium ####

## Yearly maxima ####

# Let's compute the yearly maximum for each year of the study period

# First, get rid of zeros
Season_Alex_nozeros <- filter(Season_Alex, true_count != 0)

# Then, establish periods. Based on our plots of the distribution of cells
# observed, we have a unimodal distribution in Antifer and Cabourg (i.e., 1
# period) and a somewhat bimodal distribution in Men er Roue and Ouest Loscolo 
# (i.e., 2 periods), with a very small second mode

Maxima_Alex <- Season_Alex_nozeros %>%
  # Once again, the first mode of the distribution in MeR and OL ends around 
  # fortnight 16, so let's cut the year like that:
  mutate(period = ifelse(Code_point_Libelle %in% c('Men er Roue', 'Ouest Loscolo')
                         & Day > (16*14), 2, 1))

# Let's compute the number of observations of Alexandrium by period
Maxima_Alex_n <- Maxima_Alex %>%
  # groups
  group_by(Code_point_Libelle, Year, period) %>%
  # summarise to get n
  summarise(n_max = n(), max_Alex = max(true_count), .groups = 'keep') %>%
  # Create a "shape" variable for plotting
  mutate(shape = ifelse(period == 1 & n_max > 1 & max_Alex > 1, 1, 
                        ifelse(period == 2 & n_max > 1 & max_Alex > 1, 2,
                               ifelse(period == 1 & (n_max == 1 || max_Alex == 1), 3,
                                      4))))

# Now extract maxima. Careful, this step has to be done AFTER the previous one.
Maxima_Alex <- Maxima_Alex %>%
  # Add period as a grouping variable
  group_by(Code_point_Libelle, Year, period) %>%
  # extract rows with yearly maxima of 'true_count' (slice() retains groups)
  slice(which.max(true_count)) %>%
  # Keep only Alexandrium cell densities >= 100 cells/L (more than 1 cell counted)
  filter(true_count >= 1)

# Join both datasets to get n() info on maxima
Maxima_Alex <- left_join(Maxima_Alex, Maxima_Alex_n) %>%
  # Code_point_Libelle as factor
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# Seems perfect.

# Plot this

# The df color for facet labels
df_color <- data.frame(
  name = c('Antifer ponton pétrolier', 'Cabourg',
           'Men er Roue', 'Ouest Loscolo'),
  color = c('red3', 'orangered', '#2156A1', '#5995E3'))

# plot
ggplot(Maxima_Alex, 
       aes(x = Year, y = Day, 
           shape = as.factor(shape))) +
  # then the points
  geom_point(size = 4, alpha = .8,
             color = '#711412') +
  # aesthetics
  facet_wrap_color(facets = 'Code_point_Libelle', nrow = 1,
                   colors = df_color, axes = 'all_y',
                   axis.labels = 'margins') +
  scale_y_continuous(limits = c(1, 365), 
                     breaks = c(1, 100, 200, 300, 365)) +
  # shape scale
  scale_shape_manual(values = c(16, 17, 1, 2), guide = 'none') +
  # scale_shape_manual(values = c(1, 3), guide = 'none') +
  # labels
  labs(y = 'Day of maximum Alexandrium count') +
  theme_classic()

# Ok that's nice. Now let's model that properly with some good ol' GLMs

## GLMs ####

# Very simple: our response variable is day of Alex max (DOY), our predictor is
# the Year

# We'll do GLMs, with the help of the tutorial by Dr. Bede Davies
# https://bedeffinianrowedavies.com/statisticstutorials/poissonglms

Maxima_Alex_trends <- Maxima_Alex %>%
  # We create a category variable that is the sampling site + the period
  mutate(category = paste(Code_point_Libelle, period, sep = "_")) %>%
  # Relevel as a factor (category = Sampling site + period)
  mutate(category = as_factor(category)) %>%
  mutate(category = fct_relevel(category,
                                'Antifer ponton pétrolier_1', 'Cabourg_1',
                                'Men er Roue_1', 'Men er Roue_2',
                                'Ouest Loscolo_1', 'Ouest Loscolo_2'))

# Let's transform our response variable so it follows a beta distribution
# (bounded between 0 and 1). This allows us to fit a GLM properly
Maxima_Alex_trends <- Maxima_Alex_trends %>%
  # For this we divide Day by the maximum value of this variable (365)
  mutate(Day_beta = Day/365)

# Now, we'll fit a GLM for every series of maxima (1 for each period in each 
# sampling site)
# With the beta-transformed variable (we'll use mgcv's gam() function)
glm_trends_beta <- gam(Day_beta~Year*category,
                       data = Maxima_Alex_trends,
                       # Specifying a beta distribution of the variable
                       family = betar(link="logit"))

summary(glm_trends_beta)

### Diagnostic plots ###

### Checking the model
ModelOutputs<-data.frame(Fitted=fitted(glm_trends_beta),
                         Residuals=resid(glm_trends_beta))

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- glm_trends_beta$model
# then we add the values of fitted and residuals
qq_data <- bind_cols(qq_data, ModelOutputs)

# Reorder data
qq_data_reordered <- qq_data %>%
  mutate(category = fct_relevel(category,
                                'Antifer ponton pétrolier_1', 'Cabourg_1',
                                'Men er Roue_1', 'Men er Roue_2',
                                'Ouest Loscolo_1', 'Ouest Loscolo_2')) %>%
  # we ungroup the data frame to produce only 1 qq-plot
  ungroup()

# And (qq-)plot
qqplot_custom <- ggplot(qq_data_reordered) +
  stat_qq(aes(sample=Residuals, color = category), alpha = .7) +
  stat_qq_line(aes(sample=Residuals, color = category)) +
  facet_wrap(facets = c('category')) +
  scale_color_discrete(guide = 'none') +
  theme_classic() +
  labs(y="Sample Quantiles",x="Theoretical Quantiles")

qqplot_custom

# Ok!

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data_reordered)+
  geom_point(aes(x=Fitted,y=Residuals, color = category), 
             alpha = .7) +
  facet_wrap(facets = c('category'), scales = 'free') +
  scale_color_discrete(guide = 'none') +
  theme_classic() +
  labs(y="Residuals",x="Fitted Values")

RvFplot_custom

# Okish

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data_reordered, aes(x = Residuals, 
                                                fill = category))+
  geom_histogram(binwidth = 1)+
  facet_wrap(facets = c('category'), scales = 'free') +
  scale_fill_discrete(guide = 'none') +
  theme_classic() +
  labs(x='Residuals', y = 'Count')

HistRes_custom

# It's approximately ok

### Let's predict with the GLM
NewData_1 <- expand_grid(# Year variable
  Year=seq(min(Maxima_Alex_trends$Year), 
           max(Maxima_Alex_trends$Year)),
  # We also add the site data
  category = unique(Maxima_Alex_trends$category))

Pred <- predict(glm_trends_beta, NewData_1, se.fit=TRUE, type="response")

NewData<-NewData_1 %>% 
  mutate(response=Pred$fit,
         se.fit=Pred$se.fit,
         Upr=response+(se.fit*1.96),
         Lwr=response-(se.fit*1.96))

# Join datasets of prediction and true data
Data_trends_Alex <- left_join(Maxima_Alex_trends, NewData)

#---- Plot with the GLM fits ----####
# Plot (we multiply the predictions made on the beta-transformed variable
# by 365 so that it fits the true data)
plot_trends_Alex <- ggplot(Data_trends_Alex)+
  # confidence interval of model
  geom_ribbon(aes(x = Year,
                  ymax = Upr*365,
                  ymin = Lwr*365, 
                  group = category),
              alpha=0.35, fill = '#711412') +
  # fit of model
  geom_line(aes(x = Year,
                y = response*365,
                group = category),
            color = '#711412') +
  # true data
  geom_point(aes(x = Year, y = Day,
                 shape = as_factor(shape)), color = '#711412') +
  # shape scale
  scale_shape_manual(values = c(16, 17, 1, 2), guide = 'none') +
  # y-axis scale
  scale_y_continuous(limits = c(1, 410), 
                     breaks = c(1, 100, 200, 300, 365)) +
  # Labels
  labs(title = "Alexandrium (2007-2022)",
       x = "Year", 
       y = "Day of maximum Alexandrium count") +
  # Facet wrap
  facet_wrap_color(facets = c('Code_point_Libelle'),
                   colors = df_color, nrow = 1, axes = 'all_y', 
                   axis.labels = 'margins') +
  theme_classic() +
  theme(title = element_text(size = 8),
    strip.text = element_text(
    size = 7, color = 'black'),
    strip.background = element_rect(colour = 'white'))

# Good. 2 more to go.

### Protoperidinium ####

## Yearly maxima ####

# Let's compute the yearly maximum for each year of the study period

# First, get rid of zeros
Season_Protop_nozeros <- filter(Season_Protop, true_count != 0)

# For Protoperidinium, the distribution is basically unimodal in every site. So,
# We'll not cut the year and keep only 1 period.

Maxima_Protop <- Season_Protop_nozeros %>%
  mutate(period = 1)

# Let's compute the number of observations of Protoperidinium by period
Maxima_Protop_n <- Maxima_Protop %>%
  # groups
  group_by(Code_point_Libelle, Year, period) %>%
  # summarise to get n
  summarise(n_max = n(), max_Protop = max(true_count), .groups = 'keep') %>%
  # Create a "shape" variable for plotting
  mutate(shape = ifelse(period == 1 & n_max > 1 & max_Protop > 1, 1, 
                        ifelse(period == 2 & n_max > 1 & max_Protop > 1, 2,
                               ifelse(period == 1 & (n_max == 1 || max_Protop == 1), 3,
                                      4))))

# Now extract maxima. Careful, this step has to be done AFTER the previous one.
Maxima_Protop <- Maxima_Protop %>%
  # Add period as a grouping variable
  group_by(Code_point_Libelle, Year, period) %>%
  # extract rows with yearly maxima of 'true_count' (slice() retains groups)
  slice(which.max(true_count)) %>%
  # Keep only Protoperidinium cell densities >= 100 cells/L (more than 1 cell counted)
  filter(true_count >= 1)

# Join both datasets to get n() info on maxima
Maxima_Protop <- left_join(Maxima_Protop, Maxima_Protop_n) %>%
  # Code_point_Libelle as factor
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# Seems perfect.

# Plot this

# The df color for facet labels
df_color <- data.frame(
  name = c('Antifer ponton pétrolier', 'Cabourg',
           'Men er Roue', 'Ouest Loscolo'),
  color = c('red3', 'orangered', '#2156A1', '#5995E3'))

# plot
ggplot(Maxima_Protop, 
       aes(x = Year, y = Day, 
           shape = as.factor(shape))) +
  # then the points
  geom_point(size = 4, alpha = .8,
             color = 'sienna') +
  # aesthetics
  facet_wrap_color(facets = 'Code_point_Libelle', nrow = 1,
                   colors = df_color, axes = 'all_y',
                   axis.labels = 'margins') +
  scale_y_continuous(limits = c(1, 365), 
                     breaks = c(1, 100, 200, 300, 365)) +
  # shape scale
  scale_shape_manual(values = c(16, 17, 1, 2), guide = 'none') +
  # scale_shape_manual(values = c(1, 3), guide = 'none') +
  # labels
  labs(y = 'Day of maximum Protoperidinium count') +
  theme_classic()

# Ok that's nice. Now let's model that properly with some good ol' GLMs

## GLMs ####

# Very simple: our response variable is day of Protop max (DOY), our predictor is
# the Year

# We'll do GLMs, with the help of the tutorial by Dr. Bede Davies
# https://bedeffinianrowedavies.com/statisticstutorials/poissonglms

Maxima_Protop_trends <- Maxima_Protop %>%
  # We create a category variable that is the sampling site + the period
  mutate(category = paste(Code_point_Libelle, period, sep = "_")) %>%
  # Relevel as a factor (category = Sampling site + period)
  mutate(category = as_factor(category)) %>%
  mutate(category = fct_relevel(category,
                                'Antifer ponton pétrolier_1', 'Cabourg_1',
                                'Men er Roue_1', 'Ouest Loscolo_1'))

# Let's transform our response variable so it follows a beta distribution
# (bounded between 0 and 1). This allows us to fit a GLM properly
Maxima_Protop_trends <- Maxima_Protop_trends %>%
  # For this we divide Day by the maximum value of this variable (365)
  mutate(Day_beta = Day/365)

# Now, we'll fit a GLM for every series of maxima (1 for each period in each 
# sampling site)
# With the beta-transformed variable (we'll use mgcv's gam() function)
glm_trends_beta <- gam(Day_beta~Year*category,
                       data = Maxima_Protop_trends,
                       # Specifying a beta distribution of the variable
                       family = betar(link="logit"))

summary(glm_trends_beta)

### Diagnostic plots ###

### Checking the model
ModelOutputs<-data.frame(Fitted=fitted(glm_trends_beta),
                         Residuals=resid(glm_trends_beta))

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- glm_trends_beta$model
# then we add the values of fitted and residuals
qq_data <- bind_cols(qq_data, ModelOutputs)

# Reorder data
qq_data_reordered <- qq_data %>%
  mutate(category = fct_relevel(category,
                                'Antifer ponton pétrolier_1', 'Cabourg_1',
                                'Men er Roue_1', 'Ouest Loscolo_1')) %>%
  # we ungroup the data frame to produce only 1 qq-plot
  ungroup()

# And (qq-)plot
qqplot_custom <- ggplot(qq_data_reordered) +
  stat_qq(aes(sample=Residuals, color = category), alpha = .7) +
  stat_qq_line(aes(sample=Residuals, color = category)) +
  facet_wrap(facets = c('category')) +
  scale_color_discrete(guide = 'none') +
  theme_classic() +
  labs(y="Sample Quantiles",x="Theoretical Quantiles")

qqplot_custom

# Ok!

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data_reordered)+
  geom_point(aes(x=Fitted,y=Residuals, color = category), 
             alpha = .7) +
  facet_wrap(facets = c('category'), scales = 'free') +
  scale_color_discrete(guide = 'none') +
  theme_classic() +
  labs(y="Residuals",x="Fitted Values")

RvFplot_custom

# Ok

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data_reordered, aes(x = Residuals, 
                                                fill = category))+
  geom_histogram(binwidth = 1)+
  facet_wrap(facets = c('category'), scales = 'free') +
  scale_fill_discrete(guide = 'none') +
  theme_classic() +
  labs(x='Residuals', y = 'Count')

HistRes_custom

# It's ok

### Let's predict with the GLM
NewData_1 <- expand_grid(# Year variable
  Year=seq(min(Maxima_Protop_trends$Year), 
           max(Maxima_Protop_trends$Year)),
  # We also add the site data
  category = unique(Maxima_Protop_trends$category))

Pred <- predict(glm_trends_beta, NewData_1, se.fit=TRUE, type="response")

NewData<-NewData_1 %>% 
  mutate(response=Pred$fit,
         se.fit=Pred$se.fit,
         Upr=response+(se.fit*1.96),
         Lwr=response-(se.fit*1.96))

# Join datasets of prediction and true data
Data_trends_Protop <- left_join(Maxima_Protop_trends, NewData)

#---- Plot with the GLM fits ----####
# Plot (we multiply the predictions made on the beta-transformed variable
# by 365 so that it fits the true data)
plot_trends_Protop <- ggplot(Data_trends_Protop)+
  # confidence interval of model
  geom_ribbon(aes(x = Year,
                  ymax = Upr*365,
                  ymin = Lwr*365, 
                  group = category),
              alpha=0.35, fill = 'sienna') +
  # fit of model
  geom_line(aes(x = Year,
                y = response*365,
                group = category),
            color = 'sienna') +
  # true data
  geom_point(aes(x = Year, y = Day,
                 shape = as_factor(shape)), color = 'sienna') +
  # shape scale
  scale_shape_manual(values = c(16, 17, 1, 2), guide = 'none') +
  # y-axis scale
  scale_y_continuous(limits = c(1, 410), 
                     breaks = c(1, 100, 200, 300, 365)) +
  # Labels
  labs(title = "Protoperidinium (2007-2022)",
       x = "Year", 
       y = "Day of maximum Protoperidinium count") +
  # Facet wrap
  facet_wrap_color(facets = c('Code_point_Libelle'),
                   colors = df_color, nrow = 1, axes = 'all_y', 
                   axis.labels = 'margins') +
  theme_classic() +
  theme(title = element_text(size = 8),
    strip.text = element_text(
    size = 7, color = 'black'),
    strip.background = element_rect(colour = 'white'))

# Good. Only Pseudo-nitzschia left.

### Pseudo-nitzschia ####

## Yearly maxima ####

# Let's compute the yearly maximum for each year of the study period

# First, get rid of zeros
Season_PN_nozeros <- filter(Season_PN, true_count != 0)

# Then, establish periods. Based on our plots of the distribution of cells
# observed, we have a roughly unimodal distribution in Antifer, Cabourg and Men
# er Roue, and a somewhat bimodal distribution in Ouest Loscolo 
# (i.e., 2 periods), with a very small second mode

Maxima_PN <- Season_PN_nozeros %>%
  # Once again, the first mode of the distribution in MeR and OL ends around 
  # fortnight 16, so let's cut the year like that:
  mutate(period = ifelse(Code_point_Libelle %in% c('Ouest Loscolo')
                         & Day > (16*14), 2, 1))

# Let's compute the number of observations of Pseudo-nitzschia by period
Maxima_PN_n <- Maxima_PN %>%
  # groups
  group_by(Code_point_Libelle, Year, period) %>%
  # summarise to get n
  summarise(n_max = n(), max_PN = max(true_count), .groups = 'keep') %>%
  # Create a "shape" variable for plotting
  mutate(shape = ifelse(period == 1 & n_max > 1 & max_PN > 1, 1, 
                        ifelse(period == 2 & n_max > 1 & max_PN > 1, 2,
                               ifelse(period == 1 & (n_max == 1 || max_PN == 1), 3,
                                      4))))

# Now extract maxima. Careful, this step has to be done AFTER the previous one.
Maxima_PN <- Maxima_PN %>%
  # Add period as a grouping variable
  group_by(Code_point_Libelle, Year, period) %>%
  # extract rows with yearly maxima of 'true_count' (slice() retains groups)
  slice(which.max(true_count)) %>%
  # Keep only Pseudo-nitzschia cell densities >= 100 cells/L (more than 1 cell counted)
  filter(true_count >= 1)

# Join both datasets to get n() info on maxima
Maxima_PN <- left_join(Maxima_PN, Maxima_PN_n) %>%
  # Code_point_Libelle as factor
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Men er Roue', 'Ouest Loscolo'))

# Seems perfect.

# Plot this

# The df color for facet labels
df_color <- data.frame(
  name = c('Antifer ponton pétrolier', 'Cabourg',
           'Men er Roue', 'Ouest Loscolo'),
  color = c('red3', 'orangered', '#2156A1', '#5995E3'))

# plot
ggplot(Maxima_PN, 
       aes(x = Year, y = Day, 
           shape = as.factor(shape))) +
  # then the points
  geom_point(size = 4, alpha = .8,
             color = '#649003') +
  # aesthetics
  facet_wrap_color(facets = 'Code_point_Libelle', nrow = 1,
                   colors = df_color, axes = 'all_y',
                   axis.labels = 'margins') +
  scale_y_continuous(limits = c(1, 365), 
                     breaks = c(1, 100, 200, 300, 365)) +
  # shape scale
  scale_shape_manual(values = c(16, 17, 1, 2), guide = 'none') +
  # scale_shape_manual(values = c(1, 3), guide = 'none') +
  # labels
  labs(y = 'Day of maximum Pseudo-nitzschia count') +
  theme_classic()

# Ok that's nice. Now let's model that properly with some good ol' GLMs

## GLMs ####

# Very simple: our response variable is day of PN max (DOY), our predictor is
# the Year

# We'll do GLMs, with the help of the tutorial by Dr. Bede Davies
# https://bedeffinianrowedavies.com/statisticstutorials/poissonglms

Maxima_PN_trends <- Maxima_PN %>%
  # We create a category variable that is the sampling site + the period
  mutate(category = paste(Code_point_Libelle, period, sep = "_")) %>%
  # Relevel as a factor (category = Sampling site + period)
  mutate(category = as_factor(category)) %>%
  mutate(category = fct_relevel(category,
                                'Antifer ponton pétrolier_1', 'Cabourg_1',
                                'Men er Roue_1',
                                'Ouest Loscolo_1', 'Ouest Loscolo_2'))

# Let's transform our response variable so it follows a beta distribution
# (bounded between 0 and 1). This allows us to fit a GLM properly
Maxima_PN_trends <- Maxima_PN_trends %>%
  # For this we divide Day by the maximum value of this variable (365)
  mutate(Day_beta = Day/365)

# Now, we'll fit a GLM for every series of maxima (1 for each period in each 
# sampling site)
# With the beta-transformed variable (we'll use mgcv's gam() function)
glm_trends_beta <- gam(Day_beta~Year*category,
                       data = Maxima_PN_trends,
                       # Specifying a beta distribution of the variable
                       family = betar(link="logit"))

summary(glm_trends_beta)

### Diagnostic plots ###

### Checking the model
ModelOutputs<-data.frame(Fitted=fitted(glm_trends_beta),
                         Residuals=resid(glm_trends_beta))

# We're gonna make our own qq plot with colors identifying sites
# We base it on the structure of the model
qq_data <- glm_trends_beta$model
# then we add the values of fitted and residuals
qq_data <- bind_cols(qq_data, ModelOutputs)

# Reorder data
qq_data_reordered <- qq_data %>%
  mutate(category = fct_relevel(category,
                                'Antifer ponton pétrolier_1', 'Cabourg_1',
                                'Men er Roue_1',
                                'Ouest Loscolo_1', 'Ouest Loscolo_2')) %>%
  # we ungroup the data frame to produce only 1 qq-plot
  ungroup()

# And (qq-)plot
qqplot_custom <- ggplot(qq_data_reordered) +
  stat_qq(aes(sample=Residuals, color = category), alpha = .7) +
  stat_qq_line(aes(sample=Residuals, color = category)) +
  facet_wrap(facets = c('category')) +
  scale_color_discrete(guide = 'none') +
  theme_classic() +
  labs(y="Sample Quantiles",x="Theoretical Quantiles")

qqplot_custom

# Ok

# We can do the same for residuals vs fitted
RvFplot_custom <- ggplot(qq_data_reordered)+
  geom_point(aes(x=Fitted,y=Residuals, color = category), 
             alpha = .7) +
  facet_wrap(facets = c('category'), scales = 'free') +
  scale_color_discrete(guide = 'none') +
  theme_classic() +
  labs(y="Residuals",x="Fitted Values")

RvFplot_custom

# Ok

# And let's do one last diagnostic plot with histogram of residuals
HistRes_custom <- ggplot(qq_data_reordered, aes(x = Residuals, 
                                                fill = category))+
  geom_histogram(binwidth = 1)+
  facet_wrap(facets = c('category'), scales = 'free') +
  scale_fill_discrete(guide = 'none') +
  theme_classic() +
  labs(x='Residuals', y = 'Count')

HistRes_custom

# It's quite ok

### Let's predict with the GLM
NewData_1 <- expand_grid(# Year variable
  Year=seq(min(Maxima_PN_trends$Year), 
           max(Maxima_PN_trends$Year)),
  # We also add the site data
  category = unique(Maxima_PN_trends$category))

Pred <- predict(glm_trends_beta, NewData_1, se.fit=TRUE, type="response")

NewData<-NewData_1 %>% 
  mutate(response=Pred$fit,
         se.fit=Pred$se.fit,
         Upr=response+(se.fit*1.96),
         Lwr=response-(se.fit*1.96))

# Join datasets of prediction and true data
Data_trends_PN <- left_join(Maxima_PN_trends, NewData)

#---- Plot with the GLM fits ----####
# Plot (we multiply the predictions made on the beta-transformed variable
# by 365 so that it fits the true data)
plot_trends_PN <- ggplot(Data_trends_PN)+
  # confidence interval of model
  geom_ribbon(aes(x = Year,
                  ymax = Upr*365,
                  ymin = Lwr*365, 
                  group = category),
              alpha=0.35, fill = '#649003') +
  # fit of model
  geom_line(aes(x = Year,
                y = response*365,
                group = category),
            color = '#649003') +
  # true data
  geom_point(aes(x = Year, y = Day,
                 shape = as_factor(shape)), color = '#649003') +
  # shape scale
  scale_shape_manual(values = c(16, 17, 1, 2), guide = 'none') +
  # y-axis scale
  scale_y_continuous(limits = c(1, 410), 
                     breaks = c(1, 100, 200, 300, 365)) +
  # Labels
  labs(title = "Pseudo-nitzschia (2007-2022)",
       x = "Year", 
       y = "Day of maximum Pseudo-nitzschia count") +
  # Facet wrap
  facet_wrap_color(facets = c('Code_point_Libelle'),
                   colors = df_color, nrow = 1, axes = 'all_y', 
                   axis.labels = 'margins') +
  theme_classic() +
  theme(title = element_text(size = 8),
    strip.text = element_text(
    size = 7, color = 'black'),
    strip.background = element_rect(colour = 'white'))

# Good! How do these 4 plots look like side by side?

### Arrange all 4 plots ####

# Now let's arrange all 4 plots so they fit all in 1 figure.
trendplots4 <- ggarrange(plot_trends_Dino + rremove("ylab") + rremove("xlab"),
                       plot_trends_Alex + rremove("ylab") + rremove("xlab"),
                       plot_trends_Protop + rremove("ylab") + rremove("xlab"),
                       plot_trends_PN + rremove("ylab") + rremove("xlab"),
                       nrow = 4, align = 'v')

annotate_figure(trendplots4, 
                left = textGrob('Day of maximum cells observed', 
                                rot = 90, vjust = 1,
                                gp=gpar(fontsize=10, col="black")),
                bottom = textGrob('Year',
                                  gp=gpar(fontsize=10, col="black")))

# Save the plot
# ggsave('Appendix/Plots/Microplankton_genera/FigA02_Trend_plots_4genera.tiff',
#        height = 180, width = 164, units = 'mm', compression = 'lzw',
#        dpi = 500)

# Great stuff. Mission accomplished.

####-------------------------- End of script ------------------------------####
