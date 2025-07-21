#### Starvation experiment (2016) - cell counts ###
## Experiment and data collection by T. Lacour, Script by V. Pochic
# 2025-07-21

# packages ####
library(tidyverse)
library(stringr)

# Import data ####
Exp_data  <- read.csv2('Data/Starvation_Experiment/2016/Cell_counts_starv_exp_2016.csv',
                       header = TRUE, fileEncoding = 'ISO-8859-1')

# Curate data ####

# The data contains 3 conditions, each one in triplicate.
# First condition: N (Fed). The Dinophysis in this condition were fed at day 0 of
# the experiment, and were starved after.
# Second condition: NN (Unfed). The Dinophysis in this condition were not fed at
# the beginning of the experiment (i.e., they started the experiment being
# already starved)
# Third condition: R (Refed). This is a subsample of the N condition (1 for each
# replicate), for which the Dinophysis were fed again at day 43 of the experiment

# Let's get the cell counts in one column and the condition (and replicates) in
# another.

Exp_data_pivot <- pivot_longer(Exp_data, cols = c(N1, N2, N3, NN1, NN2, NN3,
                                                  R1, R2, R3),
                               values_to = 'Cell_density',
                               names_to = 'ID') %>%
  mutate(Condition = str_extract(ID, '[A-Z]+')) %>%
  mutate(Replicate = str_extract(ID, '[0-9]+'))

# Good

# Now, we'll compute the mean and standard deviation for each condition at each
# time point

Exp_data_stats <- Exp_data_pivot %>%
  group_by(Day, Condition) %>%
  summarise(celldens.mean = mean(Cell_density, na.rm = TRUE), 
            celldens.sd = sd(Cell_density, na.rm = TRUE),
            .groups = 'keep')

# Plot data ####

# Let's plot only N and R (NN is less interesting imo), and log-transform values
# so we can plot them on a log scale
Exp_data_plot <- filter(Exp_data_pivot, Condition != 'NN') %>%
  mutate(log.celldens = log(Cell_density))
Exp_data_stats_plot <- filter(Exp_data_stats, Condition != 'NN') %>%
  mutate(log.celldens.mean = log(celldens.mean))

# a color palette
palette_exp2 <- c('black', 'red2')

ggplot(data = Exp_data_stats_plot) +
  # First, the individual points for each replicate
  geom_point(data = Exp_data_plot, aes(x = Day, y = log.celldens, 
                                       color = Condition), alpha = .4,
             size = 2.5, shape = 21, fill = 'transparent') +
  # then, mean +- sd
  geom_errorbar(aes(x = Day, ymin = log(celldens.mean - celldens.sd),
                    ymax = log(celldens.mean + celldens.sd)), linewidth = .3,
                color = 'black', width = .5) +
  geom_point(aes(x = Day, y = log.celldens.mean, color = Condition), shape = 21,
             fill = 'transparent', stroke = .8, size = 3.5) +
  # color scale
  scale_color_discrete(type = palette_exp2, guide = 'none') +
  # y-axis scale
  scale_y_continuous(limits = c(log(100), log(20000)), 
                     breaks = c(log(100), log(1000), log(5000), log(10000)),
                     labels = c('100', '1000', '5000', '10000')) +
  # labels
  labs(x = 'Time (days)', y = 'Cell density (cells per mL)') +
  # theme
  theme_classic()

# Nice! Save that
# ggsave('Plots/Experiments/Starvation_exp_2016_celldens.tiff',
#        dpi = 300, height = 85, width = 164, units = 'mm', compression = 'lzw')
