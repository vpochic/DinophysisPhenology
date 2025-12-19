#### Dinophysis in the mist ##
## Where is Dinophysis in winter?
## Looking at eDNA

# V. POCHIC
# 2025-10-30 en la Patagonia!

#### Packages ####
library(tidyverse)
library(paletteer)
library(RColorBrewer)

#### Import data ####
# Table with ASV counts
Table_ASV <- read_tsv('Data/table_with_taxo_for_stats.tsv', col_names = TRUE)

# Table with metadata
Table_metadata <- read_tsv('Data/metadata_stats2.tsv', col_names = TRUE)

#### Tidy data ####

# We need to pivot the table (to have ASVs and taxonomy as columns 
# and samples as row)
Table_ASV_pivot <- pivot_longer(Table_ASV, cols = -c('ASV_ID', 'taxonomy'), 
                                names_to = 'Site_code', values_to = 'ASV_count')

Table_ASV_pivot2 <- pivot_wider(Table_ASV_pivot, 
                                names_from = c('ASV_ID','taxonomy'),
                                values_from = 'ASV_count')

# The name combining ASV_ID and taxonomy is absolutely awful
# (But at least i's complete and accurate)

# Compute total read count for each sample
Table_summed <- Table_ASV_pivot2 %>%
  mutate(Total_eukaryote_counts = rowSums(across(-c('Site_code'))))

# Select the ASV column for Dinophysis
Table_pivot_select <- select(Table_summed, contains('Dinophysis') | 
                               contains('Site_code') |
                               contains('Total_eukaryote_counts'))

# Give a proper name to the first column
colnames(Table_pivot_select) <- c('Dinophysis_reads', 'Site_code',
                                  'Total_eukaryote_reads')

### Join the metadata
Table_ASV_metadata <- left_join(Table_pivot_select, Table_metadata,
                                by = 'Site_code')

### And select only Ouest Loscolo and Basse Michaud
Table_ASV_metadata <- filter(Table_ASV_metadata, Site == 'Ouest_Loscolo' |
                               Site == 'Basse_Michaud') %>%
  # Get a proper date value
  mutate(Date = ymd(Date)) %>%
  mutate(Year = year(Date)) %>%
  # Transform Site and Depth into factors
  mutate(Site = as_factor(Site)) %>%
  mutate(Depth = as_factor(Depth))

### Treat data to obtain percentage of reads 

MetaB_Dino <- Table_ASV_metadata %>%
  mutate(Fraction_Dinophysis_reads = Dinophysis_reads/Total_eukaryote_reads)

## This seems rather nice!
# Just need to modify the names of some variables so it is correct and in English
MetaB_Dino <- MetaB_Dino %>%
  mutate(Site = ifelse(Site == 'Ouest_Loscolo', 'Ouest Loscolo', 'Basse Michaud')) %>%
  mutate(Depth = ifelse(Depth == 'Mi-profondeur', 'Mid-depth',
                        ifelse(Depth == 'Fond', 'Bottom', 'Surface')))

# Clear tables we don't need anymore
rm(list = c('Table_ASV', 'Table_ASV_pivot', 'Table_ASV_pivot2',
            'Table_metadata', 'Table_pivot_select', 'Table_summed',
            'Table_ASV_metadata'))

#### Now let's look at that data ####

### Relevel factors

MetaB_Dino <- MetaB_Dino %>%
  mutate(Site = fct_relevel(Site, c('Basse Michaud',
                              'Ouest Loscolo'))) %>%
  mutate(Depth = fct_relevel(Depth, c('Surface',
                               'Mid-depth',
                               'Bottom'))) %>%
  # Create a dummy variable for Dinophysis presence
  mutate(Dino = ifelse(Fraction_Dinophysis_reads != 0, 'present', 'absent'))

shape_Dino <- c(4, 21)

# Plot
ggplot(MetaB_Dino, aes(x = Date, y = Depth,
                       size = 5+exp(Fraction_Dinophysis_reads),
                       shape = Dino)) +
  geom_point(fill = 'grey45', alpha = .7, stroke = .5) +
  facet_wrap(facets = 'Site') +
  scale_y_discrete(limits = rev) +
  scale_x_date(date_breaks = '3 months', date_labels = '%Y/%m',
               minor_breaks = '2 weeks') +
  scale_shape_manual(values = shape_Dino) +
  guides(size = 'none', 
         shape = guide_legend('Dinophysis eDNA detection')) +
  labs(y = 'Depth level') +
  # Theme
  theme_classic() +
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
        legend.text = element_text(size = 8, color = 'grey5'),
        legend.frame = element_rect(linewidth = .5, color = 'grey10'),
        legend.ticks = element_line(linewidth = .2, color = 'grey25'),
        legend.position = 'bottom',
        # Panel
        panel.grid = element_blank(),
        panel.background = element_rect(linewidth = .5, color = 'grey10', 
                                        fill = 'transparent'),
        # Facet labels
        strip.background = element_rect(fill = 'transparent',
                                        linewidth = 1,
                                        color = 'grey10'),
        strip.text = element_text(color = 'grey5', size = 7.5))

# Save this (approximately good-looking) plot
# ggsave('eDNA_Dinophysis_OL-BM4.tiff', height = 150, width = 250, units = 'mm',
#        compression = 'lzw', dpi = 300)
