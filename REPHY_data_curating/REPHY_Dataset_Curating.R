#### Script REPHY : curating data
# V. POCHIC 2024/08/19

# Packages

library(tidyverse)
library(stringr)
library(FactoMineR)
library(factoextra)
library(vegan)
library(paletteer)
library(mgcv)
library(gratia)
library(scatterpie)


#### Import data ####

# Open with the correct file encoding allows to preserve accents
DataREPHY_MA <- read.csv2('REPHY_Manche_Atlantique_1987-2022.csv', fileEncoding = "ISO-8859-1")
DataREPHY_Med <- read.csv2('REPHY_Med_1987-2022.csv', fileEncoding = "ISO-8859-1")

# Creating a vector to import all columns as characters
classes_char <- rep('character', 56)
# CAREFUL! ALL 2023 DATA LACK TEMPERATURE! Problem seen with Maud on 2024/01/09. To be continued.
DataREPHY_2023 <- read.csv2('Extraction SEANOE_REPHY_phyto-Hydro_Manche_Atl-Med validé 19122023.csv', 
                            fileEncoding = "ISO-8859-1", colClasses = classes_char)

# Merge the first 2 datasets
DataREPHY <- bind_rows(DataREPHY_MA, DataREPHY_Med) %>%
  # This line allows to remove the problematic row that separates the hydro and phyto datasets
  filter(Passage...Mois != 'Passage : Mois')

# Binding the 2 datasets
DataREPHY_8723 <- bind_rows(DataREPHY, DataREPHY_2023)


### Load tables used to enrich the dataset
# Load table "Zones_marines"
ZM <- read.csv('Zones_marines.csv', sep = ';', header = TRUE)

# Load table "Liste_phylum.classe"
PhyClasse <- read.csv('Liste_phylum.classe_REPHY.csv', sep =';', header = TRUE, fileEncoding = 'ISO-8859-1')

#### Formatting the dataframe with better column names ####

# Extracting the numeric code for the ZM
Table1 <- DataREPHY_8723 %>%
  mutate(ZM_Quadrige_Numero = as.numeric(str_extract(Lieu.de.surveillance...Entité.de.classement...Libellé, '[:alnum:]+')))

# Change column names (to match ZM)
colnames(Table1)[which(names(Table1) == "Lieu.de.surveillance...Mnémonique")] <- "Code_point_Mnemonique"
colnames(Table1)[which(names(Table1) == "Lieu.de.surveillance...Libellé")] <- "Code_point_Libelle"
colnames(Table1)[which(names(Table1) == "Passage...Date")] <- "Date"
colnames(Table1)[which(names(Table1) == "Coordonnées.passage...Coordonnées.minx")] <- "lon"
colnames(Table1)[which(names(Table1) == "Coordonnées.passage...Coordonnées.miny")] <- "lat"
colnames(Table1)[which(names(Table1) == "Résultat...Libellé.unité.de.mesure.associé.au.quintuplet")] <- "Mesure_Unite"
colnames(Table1)[which(names(Table1) == "Résultat...Symbole.unité.de.mesure.associé.au.quintuplet")] <- "Mesure_Symbole"
colnames(Table1)[which(names(Table1) == "Résultat...Nom.du.taxon.référent")] <- "Taxon"
colnames(Table1)[which(names(Table1) == "Résultat...Libellé.du.groupe.de.taxon")] <- "Groupe_Taxon"
colnames(Table1)[which(names(Table1) == "Résultat...Valeur.de.la.mesure")] <- "Valeur_mesure"
colnames(Table1)[which(names(Table1) == "Prélèvement...Immersion")] <- "Profondeur.metre"
colnames(Table1)[which(names(Table1) == "Prélèvement...Niveau")] <- "Prelevement.niveau"
colnames(Table1)[which(names(Table1) == "Résultat...Code.paramètre")] <- "Code.parametre"
colnames(Table1)[which(names(Table1) == "Résultat...Libellé.paramètre")] <- "Parametre"
colnames(Table1)[which(names(Table1) == "Résultat...Niveau.de.qualité")] <- "Qualite.resultat"
colnames(Table1)[which(names(Table1) == "Prélèvement...Niveau.de.qualité")] <- "Qualite.prelevement"
colnames(Table1)[which(names(Table1) == "Passage...Service.saisisseur...Libellé")] <- "Service.saisie"
colnames(Table1)[which(names(Table1) == "Passage...Heure")] <- "Heure"
colnames(Table1)[which(names(Table1) == "Résultat...Service.analyste...Code")] <- "Service.analyse"
colnames(Table1)[which(names(Table1) == "Prélèvement...Service.préleveur...Code")] <- "Service.prelevement"
colnames(Table1)[which(names(Table1) == "Prélèvement...Identifiant.interne")] <- "ID.interne.prelevement"
colnames(Table1)[which(names(Table1) == "Passage...Identifiant.interne")] <- "ID.interne.passage"

#### Curate table to keep only desired variables ####
Table1 <- Table1 %>%
  dplyr::select(c('ZM_Quadrige_Numero', 'Code_point_Mnemonique', 'Code_point_Libelle', 'Date', 
                  'Heure', 'lon', 'lat', 'Mesure_Unite', 'Mesure_Symbole', 'Taxon', 'Valeur_mesure', 
                  'Prelevement.niveau', 'Profondeur.metre', 'Code.parametre', 'Parametre', 
                  'Qualite.prelevement', 'Qualite.resultat', 'ID.interne.prelevement', 'ID.interne.passage'))

# Modifying date format so that it gives the year and month, and getting rid of rows with no Year value
Table1 <- Table1 %>%
  mutate(Date = dmy(Date)) %>%
  # modifies the date format
  mutate(Day = day(Date)) %>%
  mutate(Month = month(Date, label = F)) %>%
  mutate(Year = year(Date)) %>%
  filter(!is.na(Year))

# Transform the measured values into numbers
# Caution! The decimal separator is a comma, 
# we need first to transform it to a full stop to avoid fuckery
Table1$Valeur_mesure <-  str_replace_all(Table1$Valeur_mesure, ',', '.')
Table1 <- Table1 %>%
  mutate(Valeur_mesure = as.numeric(Valeur_mesure))

#### Tidying table structure ####

## Associate a region with each ZM code
Table1 <- left_join(Table1, ZM, by='ZM_Quadrige_Numero', suffix=c('',''))

# Basically, we want the hydrological measurements as columns and the phytoplankton taxa as rows

# Separate the table into 2 : 1 for hydrology and the other for phytoplankton
Table1_hydro <- Table1 %>%
  filter(Taxon == "")

Table1_phyto <- Table1 %>%
  filter(Taxon != "")

### Manipulating the phyto table ###
# Associate the Phylum/Class to the taxon
Table1_phyto <- left_join(Table1_phyto, PhyClasse, by='Taxon', suffix=c('',''))

# Compute genus from Taxon
# We find the genus by selecting only the first word in the string 'Taxon'
Table1_phyto <- Table1_phyto %>%
  mutate(Genus = str_extract(Taxon, 'Pseudo-nitzschia|[:alnum:]+'))

# Applying quality control and selecting only FLORTOT
Table1_phyto_select <- Table1_phyto %>%
  #filter out Code.Region = 0 (only a few mysterious events in 2011-2013)
  filter(Code.Region != 0) %>%
  #filter(Qualite.prelevement == 'Bon') %>%
  filter(Code.parametre == 'FLORTOT' | Code.parametre == 'FLORIND' | Code.parametre == 'FLORPAR') %>%
  # Keeping only surface (or near surface) sampling
  filter(Prelevement.niveau == 'Surface (0-1m)' | Prelevement.niveau == '2 metres' |
           Prelevement.niveau == 'de 3 à 5 metres')
# We're cutting the pipe here so we can create 3 tables: 1 grouped by genus, 1 by class and 1 by taxon

Table1_phyto_genus <- Table1_phyto_select %>%
  group_by(Code.Region, Code_point_Libelle, lon, lat, Year, Month, Date, Heure,
           Code.parametre, ID.interne.passage, Genus) %>%
  summarise(Comptage = sum(Valeur_mesure), .groups = 'keep') #%>%
# pivot_wider(names_from = Genus, values_from = Valeur_mesure)

Table1_phyto_class <- left_join(Table1_phyto_select, PhyClasse, by = 'Taxon', suffix = c('','')) %>%
  group_by(Code.Region, Code_point_Libelle, lon, lat, Year, Month, Date, Heure,
           Code.parametre, ID.interne.passage, Phylum.Classe) %>%
  summarise(Comptage = sum(Valeur_mesure), .groups = 'keep') #%>%
# pivot_wider(names_from = Phylum.Classe, values_from = Valeur_mesure)

Table1_phyto_taxon <- Table1_phyto_select %>%
  group_by(Code.Region, Code_point_Libelle, lon, lat, Year, Month, Date, Heure,
           Code.parametre, ID.interne.passage, Taxon) %>%
  summarise(Comptage = sum(Valeur_mesure), .groups = 'keep') #%>%

# write the tables to free some memory space
# write.csv2(Table1_phyto_genus, 'Table1_phyto_genus.csv', row.names = FALSE, fileEncoding = "ISO-8859-1")
# write.csv2(Table1_phyto_class, 'Table1_phyto_class.csv', row.names = FALSE, fileEncoding = "ISO-8859-1")
# write.csv2(Table1_phyto_taxon, 'Table1_phyto_taxon.csv', row.names = FALSE, fileEncoding = "ISO-8859-1")

### Manipulating the hydro table ###
# Spread the hydrology measurements

Table1_hydro_select <- Table1_hydro %>%
  select(c('ID.interne.passage', 'Qualite.resultat', 'Code.parametre', 'Valeur_mesure', 'Prelevement.niveau',
           'Code.Region', 'Region', 'Date', 'Day', 'Month', 'Year', 'Code_point_Libelle', 
           'Code_point_Mnemonique', 'lon', 'lat')) %>%
  #filter(Code_point_Libelle == 'Ouest Loscolo')
  #filter out Code.Region = 0 (only a few mysterious events in 2011-2013)
  filter(Code.Region != 0) %>%
  #filter(Qualite.resultat == 'Bon') %>%
  filter(Prelevement.niveau == 'Surface (0-1m)' | Prelevement.niveau == '2 metres' |
           Prelevement.niveau == 'de 3 à 5 metres') %>%
  group_by(Code.Region, Code_point_Libelle, lon, lat, Year, Month, Date, ID.interne.passage, Code.parametre) %>%
  # There is probably something shifty here, regarding multiple CHLOROA measurements at certain stations,
  # made with different methods. For now we average everything but this is quite bad.
  summarise(Valeur_mesure = mean(Valeur_mesure), .groups = 'keep') %>%
  pivot_wider(names_from = Code.parametre, values_from = Valeur_mesure)

# write the table to free some memory space
write.csv2(Table1_hydro_select, 'Table1_hydro_select.csv', row.names = FALSE, fileEncoding = "ISO-8859-1")


#### Clear the environment from objects to free some memory space ####
#### Now reload the needed tables ####

# Table with hydrology measurements
Table1_hydro_select <- read.csv2('Table1_hydro_select.csv', header = TRUE)

# Table with phytoplankton counts
Table1_phyto_genus <- read.csv2('Table1_phyto_genus.csv', header = TRUE)
Table1_phyto_class <- read.csv2('Table1_phyto_class.csv', header = TRUE)
Table1_phyto_taxon <- read.csv2('Table1_phyto_taxon.csv', header = TRUE)


# Apparently it's alright
Table1_rejoin_genus <- left_join(Table1_phyto_genus, Table1_hydro_select, by = 'ID.interne.passage', suffix = c('',''))
Table1_rejoin_class <- left_join(Table1_phyto_class, Table1_hydro_select, by = 'ID.interne.passage', suffix = c('',''))
Table1_rejoin_taxon <- left_join(Table1_phyto_taxon, Table1_hydro_select, by = 'ID.interne.passage', suffix = c('',''))

# Save table
write.csv2(Table1_rejoin_taxon, 'Table1_rejoin_taxon.csv', row.names = FALSE, 
           fileEncoding = "ISO-8859-1")

#### Chl a time series ####

# Color palettes
viridisdiverging <- c('#440154', '#443B84', '#2F6B8E', '#1F948C', '#3FBC73', '#A2DA37','#E8E419',
                      '#75D054', '#26AD81', '#25838E', '#38598C', '#482576')
zissoudiverging <- c('#3B99B1', '#59B29D', '#A2C194', '#EAC728', '#E89E16', '#E96E00','#F13D09',
                     '#E78300', '#E9AF1D', '#C9C989', '#87BC95', '#3EAAA6')

# Let's first try this at Ouest Loscolo
Table_OL_hydro <- Table1_hydro_select %>%
  filter(Code_point_Libelle == 'Ouest Loscolo') %>%
  mutate(Day = day(Date)) %>%
  mutate(month_day = as.numeric(make_date(month = Month, day = Day))) %>%
  filter(is.na(CHLOROA) == FALSE)

gam_chla_OL <- gam(CHLOROA ~ s(month_day), data = Table_OL_hydro)

ggplot(Table_OL_hydro, aes(x = month_day, y = CHLOROA)) +
  geom_smooth(method = 'gam') +
  geom_point(aes(x = month_day, y = CHLOROA, color = as.factor(Month))) +
  scale_color_manual(values = viridisdiverging) +
  geom_point() +
  labs(x = 'Date', y = 'Chl a (mg.m-3)')+
  scale_y_continuous(limits = c(0,30))+
  theme(axis.title.x = element_text(size=10), axis.text = element_text(size=10, color = 'black'), axis.title.y =element_text(size=10),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = 'bold'), 
        legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.position = 'bottom',
        panel.grid.major = element_line(color = 'gray70', linewidth = .5))+
  guides(color = guide_legend(override.aes = list(size = 5)))

# By month 
Table_OL_month_chla <- Table_OL_hydro %>%
  group_by(Month) %>%
  # Temporary solution: remove rows with CHLOROA = NA (TEMPORARY!!!)
  filter(is.na(CHLOROA) == FALSE) %>%
  summarise(CHLOROA = median(CHLOROA), .groups = 'keep')

ggplot(Table_OL_month_chla, aes(x = Month, y = CHLOROA, color = as.factor(Month))) +
  scale_color_manual(values = viridisdiverging) +
  geom_point() +
  labs(x = 'Month', y = 'Chl a (mg.m-3)')+
  theme(axis.title.x = element_text(size=10), axis.text = element_text(size=10, color = 'black'), axis.title.y =element_text(size=10),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = 'bold'), 
        legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.position = 'bottom',
        panel.grid.major = element_line(color = 'gray70', linewidth = .5))+
  guides(color = guide_legend(override.aes = list(size = 5)))

#### Class ####

Table_OL_class <- Table1_rejoin_class %>%
  filter(Code_point_Libelle == 'Le Cornard') %>%
  filter(Phylum.Classe == 'Dinophyceae' | Phylum.Classe == 'Bacillariophyceae') %>%
  filter(Comptage != 0) %>%
  mutate(Day = day(Date)) %>%
  mutate(month_day = make_date(month = Month, day = Day)) %>%
  mutate(month_day = as.numeric(make_date(month = Month, day = Day))) %>%
  mutate(Phylum.Classe = as.factor(Phylum.Classe))
  

gam_class_OL <- gam(Comptage ~ s(month_day, by=Phylum.Classe), data = Table_OL_class, method = 'REML')
draw(gam_class_OL)

sm <- smooth_estimates(gam_class_OL, data = Table_OL_class) %>%
  add_confint()

ggplot(sm, aes(x = month_day, y = est, colour = Phylum.Classe)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, colour = NULL, fill = Phylum.Classe),
              alpha = 0.2) +
  geom_line() +
  facet_wrap(~ Phylum.Classe)

ggplot(Table_OL_class, aes(x = month_day, y = Comptage, color = Phylum.Classe)) +
  geom_smooth(method='gam') +
  geom_point(alpha = .3) +
  scale_color_manual(values = c('#440154','#E96E00','#2F6B8E')) +
  #scale_shape_manual(values = c(16,17)) +
 scale_y_continuous(limits = c(0,150000)) +
 labs(x = 'Day', y = 'Comptage (cells.L-1)')+
 theme(axis.title.x = element_text(size=10), axis.text = element_text(size=10, color = 'black'), axis.title.y =element_text(size=10),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = 'bold'), 
        legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.position = 'bottom',
        panel.background = element_blank(),
        panel.grid.major = element_line(color = 'black', linewidth = .5))+
  guides(color = guide_legend(override.aes = list(size = 5)))

# By month
Table_OL_class_month <- Table_OL_class %>%
  group_by(Month, Phylum.Classe) %>%
  summarise(Comptage = median(Comptage),.groups = 'keep')

ggplot(Table_OL_class_month, aes(x = Month, y = Comptage, color = Phylum.Classe),
                               shape = Phylum.Classe) +
  geom_point(size = 4) +
  geom_line(linewidth = .3) +
  scale_color_manual(values = c('#440154','#E96E00'))+
  scale_shape_manual(values = c(16,17)) +
  #scale_y_continuous(limits = c(0,50000)) +
  labs(x = 'Month', y = 'Cells.L-1')+
  theme(axis.title.x = element_text(size=10), axis.text = element_text(size=10, color = 'black'), axis.title.y =element_text(size=10),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = 'bold'), 
        legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.position = 'bottom',
        panel.background = element_blank(),
        panel.grid.major = element_line(color = 'black', linewidth = .5))+
  guides(color = guide_legend(override.aes = list(size = 5)))

# ggsave('phenology_class_Le_Cornard.tiff', dpi = 600, device = 'tiff',
#       width = 164, height = 147.6, units = 'mm', compression = 'lzw')

# By year + month ---> This isn't that useful imo
Table_OL_class_year_month <- Table_OL_class %>%
  group_by(Year, Month, Phylum.Classe) %>%
  summarise(Comptage = median(Comptage),.groups = 'keep') %>%
  mutate(Date = make_date(month = Month, year = Year))

ggplot(Table_OL_class_year_month, aes(x = Date, y = Comptage, color = interaction(as.factor(Month), Phylum.Classe),
                                 shape = Phylum.Classe)) +
  geom_point(size = 4) +
  scale_color_manual(values = c(viridisdiverging,zissoudiverging))+
  scale_shape_manual(values = c(16,17)) +
  scale_y_continuous(limits = c(0,100000)) +
  labs(x = 'Month', y = 'Cells.L-1')+
  theme(axis.title.x = element_text(size=10), axis.text = element_text(size=10, color = 'black'), axis.title.y =element_text(size=10),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = 'bold'), 
        legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.position = 'bottom',
        panel.grid.major = element_line(color = 'gray70', linewidth = .5))+
  guides(color = guide_legend(override.aes = list(size = 5)))

# By month but only from 2011 on (last ten years)
Table_OL_class_month <- Table_OL_class %>%
  filter(Year > 2011) %>%
  group_by(Month, Phylum.Classe) %>%
  summarise(Comptage = median(Comptage), .groups = 'keep')

ggplot(Table_OL_class_month, aes(x = Month, y = Comptage, color = interaction(as.factor(Month), Phylum.Classe),
                                 shape = Phylum.Classe)) +
  geom_point(size = 4) +
  scale_color_manual(values = c(viridisdiverging,zissoudiverging))+
  scale_shape_manual(values = c(16,17)) +
  labs(x = 'Month', y = 'Cells.L-1')+
  theme(axis.title.x = element_text(size=10), axis.text = element_text(size=10, color = 'black'), axis.title.y =element_text(size=10),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = 'bold'), 
        legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.position = 'bottom',
        panel.grid.major = element_line(color = 'gray70', linewidth = .5))+
  guides(color = guide_legend(override.aes = list(size = 5)))

#### Genus ####

Table_OL_genus <- Table1_rejoin_genus %>%
  filter(Code_point_Libelle == 'Ouest Loscolo') %>%
  filter(Genus == 'Alexandrium' | Genus == 'Lepidodinium'
         | Genus == 'Dinophysis') %>%
  filter(Comptage != 0) %>%
  mutate(Day = day(Date)) %>%
  mutate(month_day = make_date(month = Month, day = Day)) %>%
  mutate(month_day = as.numeric(make_date(month = Month, day = Day)))
  

ggplot(Table_OL_genus, aes(x = month_day, y = Comptage, color = Genus,
                           shape = Genus)) +
  geom_point() +
  #scale_color_manual(values = c(viridisdiverging,zissoudiverging))+
  scale_shape_manual(values = c(16,17,18,19,20)) +
  #scale_y_continuous(limits = c(0,2000)) +
  labs(x = 'Day of year', y = 'Dinophyceae (cells.L)')+
  theme(axis.title.x = element_text(size=10), axis.text = element_text(size=10, color = 'black'), 
        axis.title.y =element_text(size=10),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = 'bold'), 
        legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.position = 'bottom',
        panel.grid.major = element_line(color = 'gray70', linewidth = .5))+
  guides(color = guide_legend(override.aes = list(size = 5)))

# By month
Table_OL_genus_month <- Table_OL_genus %>%
  group_by(Month, Genus) %>%
  summarise(Comptage = median(Comptage),.groups = 'keep') %>%
  mutate(Genus = factor(Genus, levels=c('Alexandrium', 'Lepidodinium', 'Dinophysis')))

ggplot(Table_OL_genus_month, aes(x = Month, y = log(Comptage), color = Genus)) +
  geom_point(size = 4) +
  geom_line(linewidth = .3) +
  scale_color_manual(values = c('#59B29D', '#E89E16', '#F13D09' )) + #,'#2F6B8E', '#3FBC73','#440154'
  scale_shape_manual(values = c(16,17)) +
  labs(x = 'Month', y = 'log(Cells.L-1)')+
  theme(axis.title.x = element_text(size=10), axis.text = element_text(size=10, color = 'black'), axis.title.y =element_text(size=10),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = 'bold'), 
        legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.position = 'bottom',
        panel.background = element_blank(),
        panel.grid.major = element_line(color = 'black', linewidth = .5))+
  guides(color = guide_legend(override.aes = list(size = 5)))

# ggsave('phenology_dinos_Ouest_Loscolo.tiff', device = 'tiff', dpi = 600, 
#       width = 164, height = 147.6, units = 'mm', compression = 'lzw')

# By year + month ---> This isn't that useful imo
Table_OL_genus_year_month <- Table_OL_genus %>%
  group_by(Year, Month, Genus) %>%
  summarise(Comptage = median(Comptage),.groups = 'keep') %>%
  mutate(Date = make_date(month = Month, year = Year))

ggplot(Table_OL_genus_year_month, aes(x = Date, y = Comptage, color = interaction(as.factor(Month), Genus),
                                      shape = Genus)) +
  geom_point(size = 4) +
  scale_color_manual(values = c(viridisdiverging,zissoudiverging))+
  scale_shape_manual(values = c(16,17)) +
  scale_y_continuous(limits = c(0,100000)) +
  labs(x = 'Date', y = 'Cells.L-1')+
  theme(axis.title.x = element_text(size=10), axis.text = element_text(size=10, color = 'black'), axis.title.y =element_text(size=10),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = 'bold'), 
        legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.position = 'bottom',
        panel.grid.major = element_line(color = 'gray70', linewidth = .5))+
  guides(color = guide_legend(override.aes = list(size = 5)))

# By month but only from 2011 on (last ten years)
Table_OL_genus_month <- Table_OL_genus %>%
  filter(Year > 2011) %>%
  group_by(Month, Genus) %>%
  summarise(Comptage = median(Comptage), .groups = 'keep')

ggplot(Table_OL_genus_month, aes(x = Month, y = Comptage, color = interaction(as.factor(Month), Genus),
                                 shape = Genus)) +
  geom_point(size = 4) +
  scale_color_manual(values = c(viridisdiverging,zissoudiverging))+
  scale_shape_manual(values = c(16,17)) +
  labs(x = 'Month', y = 'Cells.L-1')+
  theme(axis.title.x = element_text(size=10), axis.text = element_text(size=10, color = 'black'), axis.title.y =element_text(size=10),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = 'bold'), 
        legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.position = 'bottom',
        panel.grid.major = element_line(color = 'gray70', linewidth = .5))+
  guides(color = guide_legend(override.aes = list(size = 5)))

###### A scatter plot of T/S conditions for detection of several genera ####

### Computing the T/S "envelope" for all the taxa of interest
# Transform the data into something we can work with
Table1_rejoin_genus_select <- Table1_rejoin_genus %>%
  # Only the events for which we have Temperature AND Salinity values
  filter(is.na(TEMP) == FALSE & is.na(SALI) == FALSE & Comptage != 0) %>%
  # Select Atlantic region
  filter(Code.Region == 21 | Code.Region == 22 | Code.Region == 23) %>%
  # 2021 and 2022 out (these are the 2 years we want to analyse)
  filter(Year != 2021 & Year != 2022) %>%
  filter(Genus == 'Leptocylindrus' | Genus == 'Thalassiosira' | Genus == 'Skeletonema'
         | Genus == 'Lingulodinium' | Genus == 'Alexandrium' | Genus == 'Lepidodinium'
         | Genus == 'Protoceratium') %>%
  # Use a mutate to specify the order the genera will appear (to sort dinos from diatoms here)
  mutate(Genus = factor(Genus, levels=c('Leptocylindrus', 'Skeletonema', 'Thalassiosira', 'Alexandrium',
                                        'Lepidodinium', 'Lingulodinium', 'Protoceratium')))

# Compute a table with centroids for each taxon (from a weighted mean with abundance as the weight factor)
Table_centroids <- Table1_rejoin_genus_select %>%
  group_by(Genus) %>%
  summarise(TEMP = weighted.mean(TEMP, w = Comptage), SALI = weighted.mean(SALI, w = Comptage),
            .groups = 'keep')
# And standard deviations
Table_sd <- Table1_rejoin_genus_select %>%
  group_by(Genus) %>%
  summarise(TEMP.sd = sd(TEMP), SALI.sd = sd(SALI), .groups = 'keep')
# Joining both tables
Table_centroids <- left_join(Table_centroids, Table_sd, by = 'Genus', suffix = c('',''))


#### Computing tables for the events we are interested in (2021 and 2022) ####

# Defining a function to scale phyto cell densities between 0 and 1
# We will use this function to scale our data and produce a nice size vector
# and make fancy pie charts
scale_values <- function(x) { (x-min(x))/(max(x)-min(x))}

### Setting up a table of T, S and taxon abundance at points of interest in 2021
Table_POI_2021 <- Table1_rejoin_genus %>%
  filter(Code_point_Libelle == 'Ouest Loscolo' | Code_point_Libelle == 'Basse Michaud' |
           Code_point_Libelle == 'Bois de la Chaise large') %>%
  filter(Year == 2021) %>%
  pivot_wider(names_from = Genus, values_from = Comptage, values_fill = 0) %>%
  # Creating a 'Protoceratium' column for code consistency, because there is none observed in 2021
  # Needed for the plot (even if in the end it does not appear)
  mutate(Protoceratium = as.numeric(Lepidodinium - Lepidodinium)) %>%
  
  mutate(Other_phyto = rowSums(across(-c('Code.Region', 'Code_point_Libelle', 'lon', 'lat', 'Year', 'Month', 'Date',
                        'ID.interne.passage', 'SALI', 'TEMP', 'CHLOROA',	'TURB',	'NH4',	'PO4',	
                        'SIOH',	'OXYGENE',	'NO3.NO2',
                        'Allo',	'Anth',	'Asta',	'But.fuco',	'C.neo',	'CHLOROB',	'CHLOROC2',
                        'CHLOROC3',	'Cantha',	'Chlide.a',	'Chlide.b',	'DVCHLOROA',	'Diadino',
                        'Diato',	'Dino',	'Fuco',	'Hex.fuco',	'Lut',	'Peri',	'Phe.a',	'Pheide.a',
                        'Pras',	'Viola',	'Zea',	'beta.beta.Car',	'TURB.FNU',	'TCHLOROB',	'NANOSUP3',
                        'PEUKINF3',	'PCYAN', 'Leptocylindrus', 'Skeletonema', 'Thalassiosira', 'Alexandrium',
                        'Lepidodinium', 'Lingulodinium', 'Protoceratium')))) %>%
  mutate(Total_phyto = rowSums(across(-c('Code.Region', 'Code_point_Libelle', 'lon', 'lat', 'Year', 'Month', 'Date',
                               'ID.interne.passage', 'SALI', 'TEMP', 'CHLOROA',	'TURB',	'NH4',	'PO4',	
                               'SIOH',	'OXYGENE',	'NO3.NO2',
                               'Allo',	'Anth',	'Asta',	'But.fuco',	'C.neo',	'CHLOROB',	'CHLOROC2',
                               'CHLOROC3',	'Cantha',	'Chlide.a',	'Chlide.b',	'DVCHLOROA',	'Diadino',
                               'Diato',	'Dino',	'Fuco',	'Hex.fuco',	'Lut',	'Peri',	'Phe.a',	'Pheide.a',
                               'Pras',	'Viola',	'Zea',	'beta.beta.Car',	'TURB.FNU',	'TCHLOROB',	'NANOSUP3',
                               'PEUKINF3',	'PCYAN', 'Other_phyto')))) %>%
  mutate(Check_sum = rowSums(across(c('Leptocylindrus', 'Skeletonema', 'Thalassiosira', 'Alexandrium',
                                'Lepidodinium', 'Lingulodinium', 'Protoceratium', 'Other_phyto')))) %>%
  
  # Select only relevant variables
  select(c('Code.Region', 'Code_point_Libelle', 'lon', 'lat', 'Year', 'Month', 'Date',
                'ID.interne.passage', 'SALI', 'TEMP','CHLOROA', 'Leptocylindrus', 'Skeletonema', 'Thalassiosira', 'Alexandrium',
                'Lepidodinium', 'Lingulodinium', 'Protoceratium', 'Other_phyto', 'Total_phyto', 'Check_sum')) %>%
  
  # Additional step : to make the plot clearer, only select events where at least 1 of the 7 taxa is present
  filter(Leptocylindrus != 0 | Skeletonema != 0 | Thalassiosira != 0 | Alexandrium != 0 | Lepidodinium != 0
         | Lingulodinium != 0 | Protoceratium != 0)

# Creating a radius column for the scatterpies
Table_POI_2021$radius <- 2*(0.1+scale_values(Table_POI_2021$CHLOROA))
  

### Setting up a table of T, S and taxon abundance at points of interest in 2022
Table_POI_2022 <- Table1_rejoin_genus %>%
  filter(Code_point_Libelle == 'Ouest Loscolo' | Code_point_Libelle == 'Basse Michaud' |
           Code_point_Libelle == 'Bois de la Chaise large') %>%
  filter(Year == 2022) %>%
  pivot_wider(names_from = Genus, values_from = Comptage, values_fill = 0) %>%
  mutate(Other_phyto = rowSums(across(-c('Code.Region', 'Code_point_Libelle', 'lon', 'lat', 'Year', 'Month', 'Date',
                                         'ID.interne.passage', 'SALI', 'TEMP', 'CHLOROA',	'TURB',	'NH4',	'PO4',	
                                         'SIOH',	'OXYGENE',	'NO3.NO2',
                                         'Allo',	'Anth',	'Asta',	'But.fuco',	'C.neo',	'CHLOROB',	'CHLOROC2',
                                         'CHLOROC3',	'Cantha',	'Chlide.a',	'Chlide.b',	'DVCHLOROA',	'Diadino',
                                         'Diato',	'Dino',	'Fuco',	'Hex.fuco',	'Lut',	'Peri',	'Phe.a',	'Pheide.a',
                                         'Pras',	'Viola',	'Zea',	'beta.beta.Car',	'TURB.FNU',	'TCHLOROB',	'NANOSUP3',
                                         'PEUKINF3',	'PCYAN', 'Leptocylindrus', 'Skeletonema', 'Thalassiosira', 'Alexandrium',
                                         'Lepidodinium', 'Lingulodinium', 'Protoceratium')))) %>%
  mutate(Total_phyto = rowSums(across(-c('Code.Region', 'Code_point_Libelle', 'lon', 'lat', 'Year', 'Month', 'Date',
                                         'ID.interne.passage', 'SALI', 'TEMP', 'CHLOROA',	'TURB',	'NH4',	'PO4',	
                                         'SIOH',	'OXYGENE',	'NO3.NO2',
                                         'Allo',	'Anth',	'Asta',	'But.fuco',	'C.neo',	'CHLOROB',	'CHLOROC2',
                                         'CHLOROC3',	'Cantha',	'Chlide.a',	'Chlide.b',	'DVCHLOROA',	'Diadino',
                                         'Diato',	'Dino',	'Fuco',	'Hex.fuco',	'Lut',	'Peri',	'Phe.a',	'Pheide.a',
                                         'Pras',	'Viola',	'Zea',	'beta.beta.Car',	'TURB.FNU',	'TCHLOROB',	'NANOSUP3',
                                         'PEUKINF3',	'PCYAN', 'Other_phyto')))) %>%
  mutate(Check_sum = rowSums(across(c('Leptocylindrus', 'Skeletonema', 'Thalassiosira', 'Alexandrium',
                                      'Lepidodinium', 'Lingulodinium', 'Protoceratium', 'Other_phyto')))) %>%
  
  # Select only relevant variables
  select(c('Code.Region', 'Code_point_Libelle', 'lon', 'lat', 'Year', 'Month', 'Date',
           'ID.interne.passage', 'SALI', 'TEMP', 'CHLOROA', 'Leptocylindrus', 'Skeletonema', 'Thalassiosira', 'Alexandrium',
           'Lepidodinium', 'Lingulodinium', 'Protoceratium', 'Other_phyto', 'Total_phyto', 'Check_sum')) %>%

  # Additional step : to make the plot clearer, only select events where at least 1 of the 7 taxa is present
  filter(Leptocylindrus != 0 | Skeletonema != 0 | Thalassiosira != 0 | Alexandrium != 0 | Lepidodinium != 0
         | Lingulodinium != 0 | Protoceratium != 0)

# Creating a radius column for the scatterpies !!! Not good yet !!!
Table_POI_2022$radius <- 2*(0.1+scale_values(Table_POI_2022$CHLOROA))


##### The plot(s) #####


#### Scatterpies ####

### For 2021

ggplot(Table1_rejoin_genus_select,aes(x = SALI, y = TEMP, color = Genus)) +
  
  # Scatter plot
  geom_point(aes(size=log10(Comptage)), alpha = 0.008) +
  
  
  # POI of 2021
  geom_scatterpie(aes(x = SALI, y = TEMP, group = Code_point_Libelle, r = radius), data = Table_POI_2021,
                  cols = c('Leptocylindrus', 'Skeletonema', 'Thalassiosira', 'Alexandrium',
                           'Lepidodinium', 'Lingulodinium', 'Protoceratium', 'Other_phyto'), alpha = .75) +
  # legend for the scatterpies
  geom_scatterpie_legend(Table_POI_2021$radius, x = 27, y = 3, labeller = function(z) 
    ((z/2)-0.1)*(max(Table_POI_2021$CHLOROA)-min(Table_POI_2021$CHLOROA))+min(Table_POI_2021$CHLOROA)) +
  
  # Centroids and error bars
  geom_point(data = Table_centroids, aes(x = SALI, y = TEMP, color = Genus), size = 5) +
  geom_errorbar(data = Table_centroids, aes(ymin=TEMP-TEMP.sd, ymax=TEMP+TEMP.sd, color = Genus), width = .2,
                position=position_dodge(), linewidth = .8) +
  geom_errorbarh(data = Table_centroids, aes(xmin=SALI-SALI.sd, xmax=SALI+SALI.sd, color = Genus), height = .8,
                 position=position_dodge(), linewidth = .8) +
  
  coord_cartesian() +
  # geom_point(data = Table_POI, aes(x = SALI, y = TEMP, fill = Code_point_Libelle),
  #            shape = 17, size = 3, color = 'black') +
  
  # Ellipses
  # stat_ellipse(level = .99, linewidth = 1) +
  
  # Colors
  scale_color_manual(values = c('seagreen1', 'seagreen3', 'seagreen', 'violetred4', 'violetred3', 'violetred1', 'violet')) +
  scale_fill_manual(values = c('seagreen1', 'seagreen3', 'seagreen', 'violetred4', 'violetred3', 'violetred1', 'grey45')) +
  
  #scale_x_continuous(n.breaks = 11)+
  scale_x_continuous(breaks = seq(25, 37.5, by = 5), limits = c(25,37.5))+
  # scale_y_continuous(breaks = seq(15,35, by = 2), limits = c(15,35))+
  labs(title = "Niche T/S des taxons étudiés",
       x = "Salinité", y = "Température") +
  theme(legend.position = "bottom",legend.box = "horizontal",legend.justification = c("left", "bottom")) +
  theme(axis.title.x = element_text(size=15), axis.text = element_text(size=15, color = 'black'), axis.title.y =element_text(size=15),
        legend.text = element_text(size = 15), legend.title = element_text(size = 15, face = 'bold'), 
        legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.position = 'bottom',
        panel.grid.major = element_line(color = 'black', linewidth = .5),
        panel.background = element_blank())

# ggsave('Scatterplot_REPHY_2021bis.tiff', dpi = 600, width = 600, height = 600, units = 'mm', compression = 'lzw')


### For 2022

ggplot(Table1_rejoin_genus_select,aes(x = SALI, y = TEMP, color = Genus)) +
  # Scatter plot
  geom_point(aes(size=log10(Comptage)), alpha = 0.008) +
  
  # POI of 2022
  geom_scatterpie(aes(x = SALI, y = TEMP, group = Code_point_Libelle, r = radius), data = Table_POI_2022,
                  cols = c('Leptocylindrus', 'Skeletonema', 'Thalassiosira', 'Alexandrium',
                           'Lepidodinium', 'Lingulodinium', 'Protoceratium', 'Other_phyto'), alpha = .75) +
  # legend for the scatterpies
  geom_scatterpie_legend(Table_POI_2022$radius, x = 27, y = 3, labeller = function(z) 
    ((z/2)-0.1)*(max(Table_POI_2022$Total_phyto)-min(Table_POI_2022$Total_phyto))+min(Table_POI_2022$Total_phyto)) +
  
  # Centroids and error bars
  geom_point(data = Table_centroids, aes(x = SALI, y = TEMP, color = Genus), size = 5) +
  geom_errorbar(data = Table_centroids, aes(ymin=TEMP-TEMP.sd, ymax=TEMP+TEMP.sd, color = Genus), width = .2,
                position=position_dodge(), linewidth = .8) +
  geom_errorbarh(data = Table_centroids, aes(xmin=SALI-SALI.sd, xmax=SALI+SALI.sd, color = Genus), height = .8,
                 position=position_dodge(), linewidth = .8) +
  
  coord_cartesian() +
  # geom_point(data = Table_POI, aes(x = SALI, y = TEMP, fill = Code_point_Libelle),
  #            shape = 17, size = 3, color = 'black') +
  
  # Ellipses
  # stat_ellipse(level = .99, linewidth = 1) +
  
  # Colors
  scale_color_manual(values = c('seagreen1', 'seagreen3', 'seagreen', 'violetred4', 'violetred3', 'violetred1', 'violet')) +
  scale_fill_manual(values = c('seagreen1', 'seagreen3', 'seagreen', 'violetred4', 'violetred3', 'violetred1', 'violet', 'grey45')) +
  
  #scale_x_continuous(n.breaks = 11)+
  scale_x_continuous(breaks = seq(25, 37.5, by = 5), limits = c(25,37.5))+
  # scale_y_continuous(breaks = seq(15,35, by = 2), limits = c(15,35))+
  labs(title = "Niche T/S des taxons étudiés",
       x = "Salinité", y = "Température") +
  theme(legend.position = "bottom",legend.box = "horizontal",legend.justification = c("left", "bottom")) +
  theme(axis.title.x = element_text(size=15), axis.text = element_text(size=15, color = 'black'), axis.title.y =element_text(size=15),
        legend.text = element_text(size = 15), legend.title = element_text(size = 15, face = 'bold'), 
        legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.position = 'bottom',
        panel.grid.major = element_line(color = 'black', linewidth = .5),
        panel.background = element_blank())

# ggsave('Scatterplot_REPHY_2022bis.tiff', dpi = 600, width = 600, height = 600, units = 'mm', compression = 'lzw')

#### Dots as months ####

# Color palettes
viridisdiverging <- c('#440154', '#443B84', '#2F6B8E', '#1F948C', '#3FBC73', '#A2DA37','#E8E419',
                      '#75D054', '#26AD81', '#25838E', '#38598C', '#482576')
zissoudiverging <- c('#3B99B1', '#59B29D', '#A2C194', '#EAC728', '#E89E16', '#E96E00','#F13D09',
                     '#E78300', '#E9AF1D', '#C9C989', '#87BC95', '#3EAAA6')
seasonal <- c('#440154', '#492650', '#443B25', '#EAC728', '#E89E16', '#E96E00', 
              '#F13D15', '#9C0824', '#E12020', '#C9C989', '#87BC95', '#3EAAA6')

### For 2021

ggplot(Table1_rejoin_genus_select,aes(x = SALI, y = TEMP, color = Genus)) +
  
  # Scatter plot
  geom_point(aes(size=log10(Comptage)), alpha = 0.008) +
  
  
  # POI of 2021 with fill as month
  geom_point(data = Table_POI_2021, aes(x = SALI, y = TEMP, fill = as.factor(Month), shape = Code_point_Libelle), 
                  color = 'black', size = 8, alpha = .75) +
  
  # Centroids and error bars
  geom_point(data = Table_centroids, aes(x = SALI, y = TEMP, color = Genus), size = 8) +
  geom_errorbar(data = Table_centroids, aes(ymin=TEMP-TEMP.sd, ymax=TEMP+TEMP.sd, color = Genus), width = .5,
                position=position_dodge(), linewidth = 1.2) +
  geom_errorbarh(data = Table_centroids, aes(xmin=SALI-SALI.sd, xmax=SALI+SALI.sd, color = Genus), height = 1.2,
                 position=position_dodge(), linewidth = 1.2) +
  
  coord_cartesian() +
  # geom_point(data = Table_POI, aes(x = SALI, y = TEMP, fill = Code_point_Libelle),
  #            shape = 17, size = 3, color = 'black') +
  
  # Ellipses
  # stat_ellipse(level = .99, linewidth = 1) +
  
  # Colors
  scale_color_manual(values = c('seagreen1', 'seagreen3', 'seagreen', 'violetred4', 'violetred3', 'violetred1', 'violet')) +
  scale_fill_manual(values = seasonal) +
  scale_shape_manual(values = c(21, 22, 24)) +
  
  #scale_x_continuous(n.breaks = 11)+
  scale_x_continuous(breaks = seq(20, 37.5, by = 5), limits = c(20,37.5))+
  # scale_y_continuous(breaks = seq(15,35, by = 2), limits = c(15,35))+
  labs(title = "Prélèvements REPHY - 2021 (OL, BM, BCL)",
       x = "Salinité", y = "Température") +
  theme(legend.position = "bottom",legend.box = "horizontal",legend.justification = c("left", "bottom")) +
  theme(axis.title.x = element_text(size=25), axis.text = element_text(size=20, color = 'black'), axis.title.y =element_text(size=25),
        legend.text = element_text(size = 20), legend.title = element_text(size = 20, face = 'bold'), 
        legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.position = 'bottom',
        plot.title = element_text(size = 30),
        panel.grid.major = element_line(color = 'black', linewidth = .5),
        panel.background = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 10)), 
         size = 'none',
         fill = guide_legend('Month', override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(size = 10), nrow = 3))


# ggsave('Scatterplot_REPHY_2021month.tiff', dpi = 600, width = 600, height = 600, units = 'mm', compression = 'lzw')


### For 2022

ggplot(Table1_rejoin_genus_select,aes(x = SALI, y = TEMP, color = Genus)) +
  
  # Scatter plot
  geom_point(aes(size=log10(Comptage)), alpha = 0.008) +
  
  
  # POI of 2021 with fill as month
  geom_point(data = Table_POI_2022, aes(x = SALI, y = TEMP, fill = as.factor(Month), shape = Code_point_Libelle), 
             color = 'black', size = 8, alpha = .75) +
  
  # Centroids and error bars
  geom_point(data = Table_centroids, aes(x = SALI, y = TEMP, color = Genus), size = 8) +
  geom_errorbar(data = Table_centroids, aes(ymin=TEMP-TEMP.sd, ymax=TEMP+TEMP.sd, color = Genus), width = .5,
                position=position_dodge(), linewidth = 1.2) +
  geom_errorbarh(data = Table_centroids, aes(xmin=SALI-SALI.sd, xmax=SALI+SALI.sd, color = Genus), height = 1.2,
                 position=position_dodge(), linewidth = 1.2) +
  
  coord_cartesian() +
  # geom_point(data = Table_POI, aes(x = SALI, y = TEMP, fill = Code_point_Libelle),
  #            shape = 17, size = 3, color = 'black') +
  
  # Ellipses
  # stat_ellipse(level = .99, linewidth = 1) +
  
  # Colors
  scale_color_manual(values = c('seagreen1', 'seagreen3', 'seagreen', 'violetred4', 'violetred3', 'violetred1', 'violet')) +
  scale_fill_manual(values = seasonal) +
  scale_shape_manual(values = c(21, 22, 24)) +
  
  #scale_x_continuous(n.breaks = 11)+
  scale_x_continuous(breaks = seq(20, 37.5, by = 5), limits = c(20,37.5))+
  # scale_y_continuous(breaks = seq(15,35, by = 2), limits = c(15,35))+
  labs(title = "Prélèvements REPHY - 2022 (OL, BM, BCL)",
       x = "Salinité", y = "Température") +
  theme(legend.position = "bottom",legend.box = "horizontal",legend.justification = c("left", "bottom")) +
  theme(axis.title.x = element_text(size=25), axis.text = element_text(size=20, color = 'black'), axis.title.y =element_text(size=25),
        legend.text = element_text(size = 20), legend.title = element_text(size = 20, face = 'bold'), 
        legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.position = 'bottom',
        plot.title = element_text(size = 30),
        panel.grid.major = element_line(color = 'black', linewidth = .5),
        panel.background = element_blank()) +
        guides(color = guide_legend(override.aes = list(size = 10)), 
               size = 'none',
               fill = guide_legend('Month', override.aes = list(shape = 21)),
               shape = guide_legend(override.aes = list(size = 10), nrow = 3))

# ggsave('Scatterplot_REPHY_2022month.tiff', dpi = 600, width = 600, height = 600, units = 'mm', compression = 'lzw')

#### What if we select only some regions ####

# Select only regions of the Atlantic coast (codes 21, 22, 23)
Table1_rejoin_genus_select_Atlantic <- Table1_rejoin_genus_select %>%
  filter(Code.Region < 20 & Code.Region < 30)

# Compute a table with centroids for each taxon (from a weighted mean with log10(Abundance) as the weight factor)
Table_centroids_Atlantic <- Table1_rejoin_genus_select_Atlantic %>%
  group_by(Genus) %>%
  summarise(TEMP = weighted.mean(TEMP, w = log10(Comptage)), SALI = weighted.mean(SALI, w = log10(Comptage)),
            .groups = 'keep')
# And standard deviations
Table_sd_Atlantic <- Table1_rejoin_genus_select_Atlantic %>%
  group_by(Genus) %>%
  summarise(TEMP.sd = sd(TEMP), SALI.sd = sd(SALI), .groups = 'keep')
# Joining both tables
Table_centroids_Atlantic <- left_join(Table_centroids_Atlantic, Table_sd_Atlantic, by = 'Genus', suffix = c('',''))

ggplot(Table1_rejoin_genus_select_Atlantic,aes(x = SALI, y = TEMP, color = Genus)) +
  # Scatter plot
  geom_point(aes(size=log10(Comptage)), alpha = 0.01) +
  geom_point(data = Table_centroids_Atlantic, aes(x = SALI, y = TEMP, color = Genus), size = 5) +
  # Error bars
  geom_errorbar(data = Table_centroids_Atlantic, aes(ymin=TEMP-TEMP.sd, ymax=TEMP+TEMP.sd, color = Genus), width = 0.2,
                position=position_dodge(), linewidth = .5) +
  geom_errorbarh(data = Table_centroids_Atlantic, aes(xmin=SALI-SALI.sd, xmax=SALI+SALI.sd, color = Genus), width = 0.2,
                 position=position_dodge(), linewidth = .5) +
  # Ellipses
  stat_ellipse(level = .99, linewidth = 1) +
  #scale_x_continuous(n.breaks = 11)+
  scale_x_continuous(breaks = seq(25, 40, by = 5), limits = c(25,40))+
  # scale_y_continuous(breaks = seq(15,35, by = 2), limits = c(15,35))+
  labs(title = "Niche T/S des taxons étudiés",
       x = "Salinité", y = "Température") +
  theme(legend.position = "bottom",legend.box = "horizontal",legend.justification = c("left", "bottom")) +
  theme(axis.title.x = element_text(size=15), axis.text = element_text(size=15, color = 'black'), axis.title.y =element_text(size=20),
        legend.text = element_text(size = 15), legend.title = element_text(size = 15, face = 'bold'), 
        legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.position = 'bottom',
        panel.grid.major = element_line(color = 'white', linewidth = .5))


#### Then by taking all points in region 22 ####
Tablehydro_22 <- Table1_hydro_select %>%
  filter(Code.Region == 22) %>%
  filter(is.na(CHLOROA) == FALSE)

ggplot(Tablehydro_22, aes(x = Date, y = CHLOROA, color = Month)) +
  geom_point() +
  labs(x = 'Date', y = 'Chl a (mg.m-3)')+
  scale_y_continuous(limits = c(0,30))+
  theme(axis.title.x = element_text(size=40), axis.text = element_text(size=35, color = 'black'), axis.title.y =element_text(size=40),
        legend.text = element_text(size = 35), legend.title = element_text(size = 40, face = 'bold'), 
        legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.position = 'bottom',
        panel.grid.major = element_line(color = 'gray70', linewidth = .5))+
  guides(color = guide_legend(override.aes = list(size = 10)))

Tablehydro_OL_month <- Tablehydro_OL %>%
  group_by(Month) %>%
  summarise(chl_month_median = median(CHLOROA), groups = 'keep')

Tablehydro_OL_month_sd <- Tablehydro_OL %>%
  group_by(Month) %>%
  summarise(chl_month_sd = sd(CHLOROA), groups = 'keep')

Table_chla_OL_month <- left_join(Tablehydro_OL_month, Tablehydro_OL_month_sd, by = 'Month')

ggplot(Table_chla_OL_month, aes(x = Month, y = chl_month_median, color = Month)) +
  geom_point() +
  labs(x = 'Month', y = 'Chl a (mg.m-3)')+
  theme(axis.title.x = element_text(size=40), axis.text = element_text(size=35, color = 'black'), axis.title.y =element_text(size=40),
        legend.text = element_text(size = 35), legend.title = element_text(size = 40, face = 'bold'), 
        legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
        legend.position = 'bottom',
        panel.grid.major = element_line(color = 'gray70', linewidth = .5))+
  guides(color = guide_legend(override.aes = list(size = 10)))

#### Applying a PCA ####


### ALL THE FOLLOWING IS SHIT, PURE SHIT, ABSOLUTE SHIT

ggplot(Table1_rejoin_genus,aes(x = TEMP, y = SALI, color = Genus))+
  geom_point(aes(size=log(Comptage + 1)))+
  #scale_x_continuous(n.breaks = 11)+
  # scale_x_continuous(breaks = seq(4, 24, by = 2), limits = c(4,24))+
  # scale_y_continuous(breaks = seq(15,35, by = 2), limits = c(15,35))+
  labs(title = "Niche abiotique des taxons étudiés",
       subtitle = "pour les évènements extrêmes",
       x = "Température", y = "Salinité")+
  theme(legend.position = "bottom",legend.box = "horizontal",legend.justification = c("left", "bottom")) +
  theme(axis.title.x = element_text(size=15), axis.text = element_text(size=15, color = 'black'), axis.title.y =element_text(size=20),
      legend.text = element_text(size = 15), legend.title = element_text(size = 15, face = 'bold'), 
      legend.box.margin = margin(0, 0, 0, 0, 'cm'), legend.background = element_rect(color = 'gray10'),
      legend.position = 'bottom',
      panel.grid.major = element_line(color = 'white', linewidth = .5))

REPHY.rda <- rda(Table1_rejoin$Valeur_mesure ~ Table1_rejoin$SALI + Table1_rejoin$TEMP 
                 + Table1_rejoin$TURB, data = Table1_rejoin)

REPHY.pca <- PCA(Table1_rejoin[,29:31], scale.unit=TRUE, ncp=5, graph=T)  

fviz_pca_var(REPHY.pca, col.var='cos2', gradient.cols=c('blue','green','red'),repel=T )
fviz_eig(REPHY.pca,addlabels = T,ylim=c(0,50))
fviz_cos2(REPHY.pca,choice = "ind",axes = 1:2)
fviz_contrib(REPHY.pca,choice = "ind",axes = 1, top =12)
fviz_contrib(REPHY.pca,choice = "ind",axes = 2, top =12)
fviz_contrib(REPHY.pca,choice = "ind",axes = 1:2, top =12)
var$contrib
fviz_pca_biplot(REPHY.pca, col.var="cos2", col.ind = "cos2", repel = T,gradient.cols=c('blue','green','red'))

#### Examining chl a outliers in the dataset ####

Outliers <- read.csv('Table_outliers_Phylumclasse_Genus_taxon.csv', header = TRUE, sep = ';', dec = ',', fileEncoding = "ISO-8859-1") %>%
  group_by(Code.Region, Code_point_Libelle) %>%
  arrange(.by_groups = TRUE) %>%
  filter(is.na(CHLOROA) == FALSE) #%>%
  # pivot_longer(cols = 24:267, names_to = 'Genus', values_to = 'Comptage', values_drop_na = TRUE) #%>%
  # filter(Taxon != 'Bacillariophyceae' & Taxon != 'Dinophyceae') %>%
  # select(c('Code.Region', 'Code_point_Libelle', 'cluster', 'lon', 'lat', 'Year', 'Month', 'Date',
  #          'Prelevement.niveau', 'Code.parametre', 'Methode_CHLOROA', 'CHLOROA', 'SALI', 'TEMP',
  #          'Shannon', 'Simpson', 'BergerParker', 'Pielou', 'Outlier', 'Phy' 'Taxon', 'Comptage'))

Outliers_1 <- filter(Outliers, cluster == 1)
write.csv2(Outliers_1, 'Outliers_1.csv', row.names = FALSE)

Outliers_2 <- filter(Outliers, cluster == 2)
write.csv2(Outliers_2, 'Outliers_2.csv', row.names = FALSE)

Outliers_3 <- filter(Outliers, cluster == 3)
write.csv2(Outliers_3, 'Outliers_3.csv', row.names = FALSE)

Outliers_4 <- filter(Outliers, cluster == 4)
write.csv2(Outliers_4, 'Outliers_4.csv', row.names = FALSE)

write.csv2(Outliers, 'Outliers.csv', row.names = FALSE)
