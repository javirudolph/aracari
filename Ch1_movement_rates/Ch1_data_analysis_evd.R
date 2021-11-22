# Script for evd analysis

library(dplyr)
library(extRemes)
library(ggplot2)

# Bring data --------------------------------------------
load("Ch1_movement_rates/sims_backup/datagen_cp.RData")
load("Ch1_movement_rates/sims_backup/datagen_pp.RData")
load("Ch1_movement_rates/sims_backup/datagen_np.RData")
