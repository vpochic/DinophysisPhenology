# Dinophysis phenology - species plot #
# V. POCHIC - 2025/11/21

# Just a little script to plot the distribution of the different Dinophysis
# taxa in our dataset

## Packages ####
library(tidyverse)

## Import data ####

# Import Season_Dino (all observations of Dinophysis in our dataset of 16 sites)
Season_Dino <- read.csv2('Data/REPHY_outputs/Season_Dino_20250604.csv', header = TRUE,
                         fileEncoding = 'ISO-8859-1')

# And add the fortnight variable
Season_Dino <- Season_Dino %>%
  mutate(Date = ymd(Date)) %>%
  mutate(Week = week(Date)) %>%
  mutate(Fortnight = ceiling(Week/2))

#### Re-introducing species in the dataset ####



Season_Dino_longer <- Season_Dino %>%
  pivot_longer(cols = contains(c('Dinophysis')), 
               names_to = 'Taxon', 
               values_to = 'Comptage') %>%
  # Remove Dinophysis_genus
  filter(Taxon != 'Dinophysis_genus') %>%
  # Remove these values that were calculated on Dinophysis_genus
  select(-c('true_count', 'log_c')) %>%
  # compute the genus from the first word of 'Taxon' (here, only 'Dinophysis')
  mutate(Genus = str_extract(Taxon, 'Pseudo-nitzschia|[:alnum:]+'))

### Ok

### Now, we need to group the dataset by sampling site, fortnight, genus, and
# finally species ('taxon', more accurately)

Season_Dino_fortnight <- Season_Dino_longer %>%
  group_by(Code_point_Libelle, Fortnight, Genus, Taxon)
Season_Dino_month <- Season_Dino_longer %>%
  group_by(Code_point_Libelle, Month, Genus, Taxon)

# And now, summarise by sampling site, fortnight, genus, and taxon
Season_Dino_species_fortnight <- Season_Dino_fortnight %>%
  summarise(Comptage = sum(Comptage), .groups = 'keep')
Season_Dino_species_month <- Season_Dino_month %>%
  summarise(Comptage = sum(Comptage), .groups = 'keep')


#### A nice plot ####

# A color palette for the nine taxa, with help from our homemade Bretagne 
# color palettes :)
taxo_palette9 <- c(# Dinophysis  
  '#11203E',
  # Dinophysis acuminata-complex (acuminata, sacculus, fortii)
  '#435E7B', '#7F96B6', '#BBD4F2',
  # Dinophysis acuta
  '#8D6456',
  # Dinophysis caudata and tripos
  '#FBA823', '#FBB646',
  # Dinophysis hastata + odiosa
  '#8B064B',
  # Dinophysis + phalacroma
  '#FF6448'
)

# Changing factor order so sampling sites appear in the order we want
Season_Dino_species_fortnight <- Season_Dino_species_fortnight %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Loguivy', 'les Hébihens',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))
Season_Dino_species_month <- Season_Dino_species_month %>%
  mutate(Code_point_Libelle = as_factor(Code_point_Libelle)) %>%
  mutate(Code_point_Libelle = fct_relevel(Code_point_Libelle,
                                          'Point 1 Boulogne', 'At so',
                                          'Antifer ponton pétrolier', 'Cabourg',
                                          'Loguivy', 'les Hébihens',
                                          'Men er Roue', 'Ouest Loscolo',
                                          'Le Cornard', 'Auger',
                                          'Arcachon - Bouée 7', 'Teychan bis',
                                          'Parc Leucate 2', 'Bouzigues (a)',
                                          'Sète mer', 'Diana centre'))

# And we do the same for the Dinophysis taxa
Season_Dino_species_fortnight <- Season_Dino_species_fortnight %>%
  mutate(Taxon = as_factor(Taxon)) %>%
  mutate(Taxon = fct_relevel(Taxon,
                             'Dinophysis', 'Dinophysis.acuminata',
                             'Dinophysis.sacculus', 'Dinophysis.fortii',
                             'Dinophysis.acuta', 'Dinophysis.caudata',
                             'Dinophysis.tripos', 'Dinophysis.hastata...odiosa',
                             'Dinophysis...phalacroma')) %>%
  # Add a 'true count' variable that identifies the number of cells actually 
  # observed and counted
  mutate(true_count = Comptage/100)
Season_Dino_species_month <- Season_Dino_species_month %>%
  mutate(Taxon = as_factor(Taxon)) %>%
  mutate(Taxon = fct_relevel(Taxon,
                             'Dinophysis', 'Dinophysis.acuminata',
                             'Dinophysis.sacculus', 'Dinophysis.fortii',
                             'Dinophysis.acuta', 'Dinophysis.caudata',
                             'Dinophysis.tripos', 'Dinophysis.hastata...odiosa',
                             'Dinophysis...phalacroma')) %>%
  # Add a 'true count' variable that identifies the number of cells actually 
  # observed and counted
  mutate(true_count = Comptage/100)

### Plotting (by fortnight)
ggplot(Season_Dino_species_fortnight, aes(x = Fortnight, y = true_count, fill = Taxon)) +
  geom_col() +
  facet_wrap(facets = 'Code_point_Libelle', scales = 'free_y') +
  # axis
  scale_x_continuous(breaks = c(1, 8, 16, 20, 26),
                   labels = c('1', '8', '16', '20', '26')) +
  # Labels
  labs(title = "Observations of Dinophysis taxa (2007-2022)",
       y = "Sum of observed Dinophysis cells",
       x = 'Fortnight',
       fill = 'Taxon:'
  ) +
  # Legend labels
  scale_fill_discrete(labels = c('Dinophysis', 'D. acuminata',
                                 'D. sacculus', 'D. fortii',
                                 'D. acuta', 'D. caudata',
                                 'D. tripos', 'D. hastata + odiosa',
                                 'Dinophysis + Phalacroma'),
                      # Color palette
                      type = taxo_palette9,
                      # guide legend
                      guide = guide_legend(nrow = 3)) +
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
        legend.background = element_rect(color = 'black', linewidth = .25)
  )

# Save this awesome plot
ggsave('Plots/REPHY/Taxa_16sites_Dinophenology_fortnight.tiff', dpi = 300, width = 164, height = 164,
       units = 'mm', compression = 'lzw')

### Plotting (by month)
ggplot(Season_Dino_species_month, aes(x = Month, y = true_count, fill = Taxon)) +
  geom_col() +
  facet_wrap(facets = 'Code_point_Libelle', scales = 'free_y') +
  # axis
  scale_x_discrete(limits = c(1,2,3,4,5,6,7,8,9,10,11,12),
                   labels = c('1','2','3','4','5','6','7','8','9','10',
                              '11','12')) +
  # Labels
  labs(title = "Observations of Dinophysis taxa",
       y = "Sum of observed Dinophysis cells",
       x = 'Month',
       fill = 'Taxon'
  ) +
  # Legend labels
  scale_fill_discrete(labels = c('Dinophysis', 'D. acuminata',
                                 'D. sacculus', 'D. fortii',
                                 'D. acuta', 'D. caudata',
                                 'D. tripos', 'D. hastata + odiosa',
                                 'Dinophysis + Phalacroma'),
                      # Color palette
                      type = taxo_palette9) +
  # Theme
  theme_classic() +
  theme(legend.position = 'bottom',
        strip.text = element_text(size = 7.5),
        axis.text = element_text(size = 6),
        title = element_text(size = 13),
        axis.title = element_text(size = 7.5))

# Save this awesome plot
# ggsave('Taxa_16sites_Dinophenology_month.tiff', dpi = 300, width = 300, height = 200,
#        units = 'mm', compression = 'lzw')
