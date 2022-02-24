#################################
## CH3 script
#################################


# Libraries -----------------------------------------------
set.seed(20220201)

library(aracari)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(purrr)
library(fitdistrplus)

source("Ch3_samplesize/Ch3_functions.R")

# ORIGINAL SCENARIO --------------------------------------
scenario <- "_orig0223"
desired_means <- c(28, 32, 40, 50)
desired_sds <- c(49.7, 39.9, 33.3, 31.1)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)

## Mixture distribution --------------------------------------
# We know that real data arise from a mix of complex processes.
# We simulate our 'truth' from a mixture distribution

# These are the weights for the distribution
pis <- c( 0.1, 0.2, 0.3, 0.4)

# The color associated to each mixture component
dens_cols <- c("#264653", "#2a9d8f", "#f4a261", "#e76f51")

# Build the densities for plotting
weighted_densities <- purrr::map(1:4, function(y) stat_function(fun = dlnorm, args = list(meanlog = pars$meanlog[y], sdlog = pars$sdlog[y]), color = dens_cols[y], size=pis[y]*5, alpha = 0.8))

# Plot density components of the mixture
ggplot() +
  weighted_densities +
  theme_minimal() +
  labs(y = "Density") +
  lims(x = c(0, 150)) -> densities_plot
densities_plot

# Plot the means and sd for each of the curves
data.frame(desired_means, desired_sds, gID = factor(1:4)) %>%
  mutate(lo = desired_means - desired_sds,
         hi = desired_means + desired_sds) %>%
  ggplot(., aes(x = gID, y = desired_means, color = gID)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.1, size = 1) +
  scale_color_manual(values = dens_cols) +
  labs(y = "Mean +/- SD") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none") -> means_sd_plot
means_sd_plot

# Align these plots to have multiple panels lates
top_row <- plot_grid(densities_plot, means_sd_plot, rel_widths = c(2,1))
top_row


## Simulate data -------------------------------------------------------

# Using purrr to get samples instead of a for loop. More efficient.

# We are getting 50k samples from the mixture
samp.n <- 50000

# These category samples are based on the weights
cat.samp.sizes <- samp.n*pis

# Using purrr to pull values from each component of the mixture according to the weights
# This instead of using a for loop
purrrsampls <- tibble(gID = c(1:length(cat.samp.sizes)), pars) %>%
  mutate(x.samps = purrr::map(gID, function(y) rlnorm(cat.samp.sizes[y], pars$meanlog[y], pars$sdlog[y])))

# Making it an easier to read data frame with only the group ID (mixture components) and the sampled data
purrrsampls %>%
  dplyr::select(., gID, x.samps) %>%
  unnest(cols = c(x.samps)) -> simplsamps

save(simplsamps, file = paste0("Ch3_samplesize/simdata/mixturesamples", scenario, ".RData"))


# This is the data that we consider our "truth"

# Extract the tail (the highest values) to visualize later in the histogram since they are too few to show in the bins
data.tail <- data.frame(values = simplsamps$x.samps, y = 100) %>% arrange(desc(values)) %>% filter(values >=100)

# Histogram plus the points in the tail.
simplsamps %>%
  ggplot(., aes(x = x.samps)) +
  geom_histogram(bins = 100) +
  geom_point(data = data.tail[1:50,], aes(x = values, y = y), color = "black", alpha = 0.5) +
  labs(y = "Frequency", x = "Distance") +
  #lims(x = c(0, 150)) +
  theme_minimal() -> truth_hist
truth_hist

# Density curve for the mixture based on the 50k samples.
simplsamps %>%
  ggplot(., aes(x = x.samps)) +
  geom_density() +
  labs(y = "Density", x = "Distance") +
  lims(x = c(0, 150), y = c(0,0.05)) +
  theme_minimal() -> truth_density
truth_density

# Visualize both in one frame
plot_grid(truth_hist, truth_density)

### Figure 1 -----------------------------

bottom_row <- plot_grid(truth_hist, truth_density)

plot_grid(top_row, bottom_row,nrow = 2)

ggsave(paste0("Ch3_samplesize/Figures/Figure1", scenario,".png"))


## Thresholds ----------------

thresh.vals <- c(50, 75, 100, 150, 250, 500, 1000)
tibble(thresh.vals) %>%
  mutate(n.overthresh = map_dbl(1:length(thresh.vals),
                                function(y) length(which(simplsamps$x.samps >= thresh.vals[y]))),
         theta_t = n.overthresh/samp.n) -> theta_ts
theta_ts
















