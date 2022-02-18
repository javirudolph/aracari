##########################################
## DATA GENERATION CH3
#########################################

# Libraries -----------------------------------------------
set.seed(20220201)

library(aracari)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(purrr)
library(fitdistrplus)

## Functions ----------------------------------------------------
source("Ch3_samplesize/Ch3_functions.R")

# MIXTURE MODEL --------------------------------------------------
## Mixture proportions ----------------------------------------
# We assume four categories of individuals with increasing movement.
pis <- c( 0.1, 0.2, 0.3, 0.4)
sum(pis)
hist(pis)

## Mixture components --------------------------------------
# Since movement lengths are only positive, we use a lognormal distribution to describe them

pars <- desired_mean_sd(mu_x = c(28, 32, 40, 48), sd_x = c(49.7, 39.9, 33.3, 31.1))

# Just checking that it makes sense now
lnorm_mean_var(pars$meanlog, pars$sdlog)
pars


mus <- pars$meanlog
sigsqs <- pars$sdlog
dens_cols <- c("#264653", "#2a9d8f", "#f4a261", "#e76f51")

### PLOT the mixture components -------------------------------

# Using I function I wrote that I use to make a lot of these density curves for a plot
# can be found in the functions script.
lnorm_densities <- lnorm_densities_fx(pars$meanlog, pars$sdlog, dens_cols)

# But choosing to use the weighted curves instead
weighted_densities <- purrr::map(1:4, function(y) stat_function(fun = dlnorm,
                                                                args = list(meanlog = pars$meanlog[y], sdlog = pars$sdlog[y]),
                                                                color = dens_cols[y], size=pis[y]*10, alpha = 0.8))

ggplot() +
  #lnorm_densities +
  weighted_densities +
  theme_minimal() +
  labs(y = "Density") +
  lims(x = c(0, 150)) -> densities_plot
densities_plot


# SAMPLING --------------------------------------------------------------

## Draw from mixture ----------------------------------------------------
# Using purrr to get samples instead of a for loop. More efficient.

# We are getting 50k samples from the mixture
samp.size <- 50000

# These category samples are based on the weights
cat.samp.sizes <- samp.size*pis

# Using purrr to pull values from each component of the mixture according to the weights
purrrsampls <- tibble(gID = c(1:4), pars) %>%
  mutate(data = purrr::map(gID, function(y) rlnorm(cat.samp.sizes[y], pars$meanlog[y], pars$sdlog[y])))

# Making it an easier to read data frame with only the group ID (mixture components) and the sampled data
purrrsampls %>%
  dplyr::select(., gID, data) %>%
  unnest(cols = c(data)) -> simplsamps

# Extract the tail (the highest values) to visualize later in the histogram since they are too few to show in the bins
data.tail <- data.frame(values = simplsamps$data, y = 100) %>% arrange(desc(values)) %>% filter(values >=100)

## Visualization -----------------------------------------------------

# Histogram plus the points in the tail.
simplsamps %>%
  ggplot(., aes(x = data)) +
  geom_histogram(bins = 100) +
  geom_point(data = data.tail[1:50,], aes(x = values, y = y), color = "black", alpha = 0.5) +
  labs(y = "Frequency", x = "Distance") +
  #lims(x = c(0, 150)) +
  theme_minimal() -> sampleshist

# Density curve for the mixture based on the 50k samples.

simplsamps %>%
  ggplot(., aes(x = data)) +
  geom_density() +
  labs(y = "Density", x = "Distance") +
  lims(x = c(0, 150), y = c(0,0.05)) +
  theme_minimal() -> samplesdens

# Visualize both in one frame
plot_grid(sampleshist, samplesdens)


# Generate Figure 1 that describes the mixture components, and the samples

bottom_row <- plot_grid(sampleshist, samplesdens)
plot_grid(densities_plot, bottom_row,nrow = 2)
ggsave("Ch3_samplesize/Figure1.png")


