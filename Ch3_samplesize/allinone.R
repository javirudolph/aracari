#
# # Libraries -----------------------------------------------
# set.seed(271367)
#
# library(dplyr)
# library(ggplot2)
# library(cowplot)
# library(tidyr)
# library(purrr)
# library(fitdistrplus)
# library(extRemes)
#
#
# source("Ch3_samplesize/Ch3_functions.R")
#
# # points_for_boxplots <- 50
# # B <- 1000
#
# points_for_boxplots <- 50
# B <- 200
#
# # ORIGINAL SCENARIO --------------------------------------
# dir_scenario <- "allinone_original"
#
# if(dir.exists(paste0("Ch3_samplesize/", dir_scenario)) == FALSE){
#   dir.create(paste0("Ch3_samplesize/", dir_scenario))}
#
# desired_means <- c(28, 32, 40, 50)
# desired_sds <- c(49.7, 39.9, 33.3, 31.1)
# pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)










# Data simulation ---------------------------

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
# densities_plot

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
# means_sd_plot

# Align these plots to have multiple panels lates
top_row <- plot_grid(densities_plot, means_sd_plot, rel_widths = c(2,1))
# top_row

# Using purrr to get samples instead of a for loop. More efficient.

# We are getting 50k samples from the mixture
tru_n <- 50000
w_samp_sizes <- pis*tru_n
# Using purrr to pull values from each component of the mixture according to the weights
# This instead of using a for loop
purrrsampls <- tibble(gID = c(1:length(w_samp_sizes)), pars) %>%
  mutate(x_samps = purrr::map(gID, function(y) rlnorm(w_samp_sizes[y], pars$meanlog[y], pars$sdlog[y])))

# Making it an easier to read data frame with only the group ID (mixture components) and the sampled data
purrrsampls %>%
  dplyr::select(., gID, x_samps) %>%
  unnest(cols = c(x_samps)) -> truth_df

save(truth_df, file = paste0("Ch3_samplesize/", dir_scenario, "/mixturesamples.RData"))

# This is the data that we consider our "truth"

# Extract the tail (the highest values) to visualize later in the histogram since they are too few to show in the bins
tru_tail <- data.frame(values = truth_df$x_samps, y = 100) %>% arrange(desc(values)) %>% filter(values >=50)

# Histogram plus the points in the tail.
truth_df %>%
  ggplot(., aes(x = x_samps)) +
  geom_histogram(bins = 100) +
  geom_point(data = tru_tail[1:50,], aes(x = values, y = y), color = "black", alpha = 0.5) +
  labs(y = "Frequency", x = "Distance") +
  #lims(x = c(0, 150)) +
  theme_minimal() -> truth_hist
# truth_hist

# Density curve for the mixture based on the 50k samples.
truth_df %>%
  ggplot(., aes(x = x_samps)) +
  geom_density() +
  labs(y = "Density", x = "Distance") +
  lims(x = c(0, 150), y = c(0,0.05)) +
  theme_minimal() -> truth_density
# truth_density

# Visualize both in one frame
# plot_grid(truth_hist, truth_density)

## Figure 1 -----------------------------

bottom_row <- plot_grid(truth_hist, truth_density)

plot_grid(top_row, bottom_row,nrow = 2)

ggsave(paste0("Ch3_samplesize/", dir_scenario, "/Figure1.png"), width = 6, height = 5)


##################################################
# The ALL IN ONE

# This means, for every experiment we do, the sample
# To the same sample we fit the glm lomax, the GP and we perfrom nonpar bootstrap correction on both
# Always using the same sample.
# We shall see how that works

# Stuff needed for plots later

# TEST VALUES --------------------------------------
## Thresholds ----------------
# Thresholds based on quantiles: because when we have a sample, we still don't know the truth
# But we can know the quantiles of our sample.
# summary(truth_df$x_samps)
# thresh_tests <- round(quantile(truth_df$x_samps, c(0.5, 0.75, 0.9, 0.99, 0.999, 0.9999), names = FALSE))
tru_n <- 50000
pis <- c( 0.1, 0.2, 0.3, 0.4)
w_samp_sizes <- pis*tru_n

thresh_tests <- c(50, 100, 150, 250, 500, 750, 1000)
tibble(thresh_tests) %>%
  mutate(n_overthresh = map_dbl(1:length(thresh_tests),
                                function(y) length(which(truth_df$x_samps >= thresh_tests[y]))),
         theta = n_overthresh/tru_n) -> thetas

## Sample sizes --------------------
samp_n_tests <- c(80, 200, 500, 800, 1000, 1600)

all_test_combos <- data.frame(expand.grid(thresh_tests = thresh_tests, samp_n_tests = samp_n_tests))


nruns <- paste0("nrun_", 1:points_for_boxplots)
combo_tests <- data.frame(expand.grid(thresh_tests = thresh_tests, samp_n_tests = samp_n_tests, nruns = nruns))
allin1df <- combo_tests
boots_df <- NULL

for(i in 1:nrow(combo_tests)){
  # Get the sample
  ith_n <- combo_tests$samp_n_tests[i]
  ith_thresh <- combo_tests$thresh_tests[i]
  ith_samps <- data.frame(x_star = sample(truth_df$x_samps, ith_n))

  # Lomax GLM
  ith_glm <- lomax.glm(formula=~1, my.dataf=ith_samps, response=ith_samps$x_star)
  ith_alpha_hat <- ith_glm$alphas.hat[1]
  ith_k_hat <- ith_glm$k.hat

  # Estimate tail
  ith_theta_hat <- lomax.St(x = ith_thresh, alpha = ith_alpha_hat, k = ith_k_hat, log.scale=TRUE)
  allin1df$lomax_theta_hat[i] <- ith_theta_hat

  # Generalized Pareto
  ith_quant <- quantile(ith_samps$x_star, 0.5)
  allin1df$sampl_quant[i] <- ith_quant
  ith_evd <- fevd(ith_samps$x_star, threshold = ith_quant, type = "GP")
  allin1df$gp_theta_hat[i] <- pextRemes(ith_evd, ith_thresh, lower.tail = FALSE)

  # Nonparametric bootstrapping
  bth_df <- data.frame(n_B = 1:B)

  for(b in 1:B){
    # Sample with replacement from the given sample
    bth_samps <- data.frame(x_star = sample(ith_samps$x_star, ith_n, replace = TRUE))

    # Estimate params using glm Lomax
    bth_glm <- lomax.glm(formula=~1, my.dataf=bth_samps, response=bth_samps$x_star)
    bth_alpha_star <- bth_glm$alphas.hat[1]
    bth_k_star <- bth_glm$k.hat
    bth_df$lomax_theta_star[b] <- lomax.St(x = ith_thresh, alpha = bth_alpha_star, k = bth_k_star, log.scale=TRUE)

    # Generalized Pareto
    bth_quant <- quantile(bth_samps$x_star, 0.5)
    bth_evd <- fevd(bth_samps$x_star, threshold = bth_quant, type = "GP")
    bth_df$quant_star[b] <- bth_quant
    bth_df$gp_theta_star[b] <- pextRemes(bth_evd, ith_thresh, lower.tail = FALSE)

    bth_df$ith_row[b] <- paste(i)
    bth_df$ith_n[b] <- ith_n
    bth_df$ith_thresh <- ith_thresh

  }

  boots_df[[i]] <- bth_df

  allin1df$lomax_theta_stars[i] <- mean(bth_df$lomax_theta_star)
  allin1df$gp_theta_stars[i] <- mean(bth_df$gp_theta_star)

}

allin1df %>%
  right_join(., thetas[,c(1,3)]) %>%
  mutate(lomax_theta_bar = 2*lomax_theta_hat - lomax_theta_stars,
         gp_theta_bar = 2*gp_theta_hat - gp_theta_stars,
         lomax_ratio = lomax_theta_bar - log(theta),
         gp_ratio = log(gp_theta_bar/theta)) -> allin1df

save(allin1df, file = paste0("Ch3_samplesize/", dir_scenario, "/allin1df.RData"))

























