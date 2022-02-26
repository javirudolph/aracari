
# DATASIM -------------------------------------------------------
## Mixture distribution --------------------------------------

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


## Simulate data -------------------------------------------------------

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
tru_tail <- data.frame(values = truth_df$x_samps, y = 100) %>% arrange(desc(values)) %>% filter(values >=100)

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


# TEST VALUES --------------------------------------
## Thresholds ----------------
# Thresholds based on quantiles: because when we have a sample, we still don't know the truth
# But we can know the quantiles of our sample.
# summary(truth_df$x_samps)
# thresh_tests <- round(quantile(truth_df$x_samps, c(0.5, 0.75, 0.9, 0.99, 0.999, 0.9999), names = FALSE))

thresh_tests <- c(50, 100, 150, 250, 500, 750, 1000)
tibble(thresh_tests) %>%
  mutate(n_overthresh = map_dbl(1:length(thresh_tests),
                                function(y) length(which(truth_df$x_samps >= thresh_tests[y]))),
         theta = n_overthresh/tru_n) -> thetas
thetas

## Sample sizes --------------------
samp_n_tests <- c(80, 200, 500, 800, 1000, 1600)

all_test_combos <- data.frame(expand.grid(thresh_tests = thresh_tests, samp_n_tests = samp_n_tests))

# COMPARE TAIL METHODS ----------------------------------------------

## Estimate mles --------------------------------

mles_df <- all_test_combos

for(i in 1:nrow(mles_df)){
  ith_n <- mles_df$samp_n_tests[i]
  ith_thresh <- mles_df$thresh_tests[i]
  ith_samps <- data.frame(x_star = sample(truth_df$x_samps, ith_n))
  ith_tail <- length(which(ith_samps$x_star >= ith_thresh))
  mles_df$tail_p[i] <- ith_tail/ith_n
  
  # Fit Lomax
  optim_out <- optim(par=log(c(1.5, 1.5)), fn=nllike.simp, method="BFGS", Y=ith_samps$x_star)
  mles_star <- exp(optim_out$par)
  alpha_star <- mles_star[1]
  k_star <- mles_star[2]
  mles_df$alpha_star[i] <-alpha_star
  mles_df$k_star[i] <- k_star
  
  # Check with Lomax GLM
  glm_out <- lomax.glm(formula=~1, my.dataf=ith_samps, response=ith_samps$x_star)
  alpha.2 <- glm_out$alphas.hat[1]
  k.2     <- glm_out$k.hat
  mles_df$alpha_glm[i] <-alpha.2
  mles_df$k_glm[i] <- k.2
  
  # Check with GP fit
  ith_evd <- fevd(ith_samps$x_star, threshold = 0, type = "GP")
  gp_scale <- summary(ith_evd)$par[1]
  gp_shape <- summary(ith_evd)$par[2]
  mles_df$alpha_GP[i] <- gp_scale*gp_shape
  mles_df$k_GP[i] <- 1/gp_shape
  mles_df$scale[i] <- gp_scale
  mles_df$shape[i] <- gp_shape
  
  
  
  # Estimate the tail
  #Lomax simple
  log_St_star <- lomax.St(x = ith_thresh,alpha = alpha_star,k = k_star,log.scale=TRUE)
  mles_df$log_St_hat[i] <- log_St_star
  
  # Lomax glm
  log.st.star2 <- lomax.St(x=ith_thresh,alpha=alpha.2,k=k.2,log.scale=TRUE)
  mles_df$log_St_hat2[i] <- log.st.star2
  
  # GP tail
  mles_df$GP_theta_hat[i] <- pextRemes(ith_evd, ith_thresh, lower.tail = FALSE)
  
  
}

mles_df %>%
  right_join(., thetas[, c(1,3)]) -> mles_df



### Fig2 ---------------------------------------------
# Compare the approximation (theta_hat) of the tail to the truth (theta) for all three methods
# As the threshold increases, the variation in the estimation of the threshold increases
# Makes sense. It's harder to get precise estimates for very rare events

mles_df %>%
  ggplot(., aes(x = thresh_tests, y = log_St_hat - log(theta))) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept=  0), color = "red") +
  labs(x = "Threshold Value", y = "Log St - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() -> f2_lomax

mles_df %>%
  ggplot(., aes(x = thresh_tests, y = log_St_hat2 - log(theta))) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept=  0), color = "red") +
  labs(x = "Threshold Value", y = "Log St.glm - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() -> f2_glm

mles_df %>%
  ggplot(., aes(x = thresh_tests, y = GP_theta_hat/theta)) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept=  1), color = "red") +
  labs(x = "Threshold Value", y = "GP theta hat / theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() -> f2_evd

plot_grid(f2_lomax, f2_glm, f2_evd, nrow = 3)
ggsave(paste0("Ch3_samplesize/", dir_scenario, "/Figure2.png"), width = 8, height = 8)


### Fig3 ------------------
# Q: How does this change with sample size?
# I would expect that by increasing the sample size, we get closer to the true estimator

palette4_sampsize <- c("#95190C", "#610345", "#E3B505")

mles_df %>%
  ggplot(., aes(x = thresh_tests, y = log_St_hat - log(theta))) +
  geom_point(size = 3, alpha = 0.7, aes(color = samp_n_tests)) +
  geom_hline(aes(yintercept=  0), color = palette4_sampsize[1]) +
  scale_color_gradient(low = palette4_sampsize[2], high = palette4_sampsize[3]) +
  labs(x = "Threshold Value", y = "Log St - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none") -> f3_lomax

mles_df %>%
  ggplot(., aes(x = thresh_tests, y = log_St_hat2 - log(theta))) +
  geom_point(size = 3, alpha = 0.7, aes(color = samp_n_tests)) +
  geom_hline(aes(yintercept=  0), color = palette4_sampsize[1]) +
  scale_color_gradient(low = palette4_sampsize[2], high = palette4_sampsize[3]) +
  labs(x = "Threshold Value", y = "Log St.glm - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none") -> f3_glm

mles_df %>%
  ggplot(., aes(x = thresh_tests, y = GP_theta_hat/theta)) +
  geom_point(size = 3, alpha = 0.7, aes(color = samp_n_tests)) +
  geom_hline(aes(yintercept=  1), color = palette4_sampsize[1]) +
  scale_color_gradient(low = palette4_sampsize[2], high = palette4_sampsize[3]) +
  labs(x = "Threshold Value", y = "GP theta hat / theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_colorbar(barheight = 0.5, barwidth = 15, title = "Sample \n Size")) -> f3_evd

pleg <- get_legend(f3_evd)

plot_grid(f3_lomax, f3_glm, f3_evd + theme(legend.position = "none"), pleg, nrow = 4, rel_heights = c(1,1,1, 0.3))
ggsave(paste0("Ch3_samplesize/", dir_scenario, "/Figure3.png"), width = 8, height = 8)

## Repeat estimation n times ----------------------------------------------------

nreps_mles_df <- bxplt_data_fx(points_for_boxplots)

save(nreps_mles_df, file = paste0("Ch3_samplesize/", dir_scenario, "/nreps_mles_df.RData"))

## Fig4 bxplts ratio ----------------------------
# Boxplots to see what's going on with the estimation - compare estimate and truth

nreps_mles_df %>%
  mutate(ratio = log_St_hat - log(theta)) %>%
  ggplot(., aes(y = ratio, group = thresh_tests, color = samp_n_tests)) +
  facet_wrap(~samp_n_tests, nrow = 1) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(aes(yintercept=  0), color = palette4_sampsize[1]) +
  scale_color_gradient(low = palette4_sampsize[2], high = palette4_sampsize[3]) +
  labs(x = "Threshold Value", y = "Log St - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none") -> f4_lomax

nreps_mles_df %>%
  ggplot(., aes(y = log_St_hat2 - log(theta), group = thresh_tests, color = samp_n_tests)) +
  facet_wrap(~samp_n_tests, nrow = 1) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(aes(yintercept=  0), color = palette4_sampsize[1]) +
  scale_color_gradient(low = palette4_sampsize[2], high = palette4_sampsize[3]) +
  labs(x = "Threshold Value", y = "Log St.glm - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none") -> f4_glm

# nreps_mles_df %>%
#   ggplot(., aes(y = GP_theta_hat/theta, group = thresh_tests, color = samp_n_tests)) +
#   facet_wrap(~samp_n_tests, nrow = 1) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_hline(aes(yintercept=  1), color = palette4_sampsize[1]) +
#   scale_color_gradient(low = palette4_sampsize[2], high = palette4_sampsize[3]) +
#   labs(x = "Threshold Value", y = "GP theta hat/theta") +
#   scale_x_continuous(breaks = thresh_tests) +
#   scale_y_continuous(limits = c(-0.5, 3)) +
#   theme_bw() +
#   theme(legend.position = "bottom") +
#   guides(color = guide_colorbar(barheight = 0.5, barwidth = 15, title = "Sample \n Size"))  -> f4_evd

nreps_mles_df %>%
  filter(., GP_theta_hat != 0) %>% 
  ggplot(., aes(y = log(GP_theta_hat/theta), group = thresh_tests, color = samp_n_tests)) +
  facet_wrap(~samp_n_tests, nrow = 1) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(aes(yintercept=  0), color = palette4_sampsize[1]) +
  scale_color_gradient(low = palette4_sampsize[2], high = palette4_sampsize[3]) +
  labs(x = "Threshold Value", y = "Log GP - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_colorbar(barheight = 0.5, barwidth = 15, title = "Sample \n Size"))  -> f4_evd

pleg <- get_legend(f4_evd)

plot_grid(f4_lomax, f4_glm, f4_evd + theme(legend.position = "none"), pleg, nrow = 4, rel_heights = c(1,1,1, 0.3))
ggsave(paste0("Ch3_samplesize/", dir_scenario, "/Figure4.png"), width = 8, height = 8)

## Fig5 jitter ratio -------------------------------------
# So, the boxplot can hide important trends, so we explore the values themselves
nreps_mles_df %>%
  ggplot(., aes(y = log_St_hat - log(theta), color = samp_n_tests)) +
  facet_wrap(~samp_n_tests, nrow = 1) +
  geom_jitter(aes(x = thresh_tests), alpha = 0.5) +
  geom_hline(aes(yintercept=  0), color = palette4_sampsize[1]) +
  scale_color_gradient(low = palette4_sampsize[2], high = palette4_sampsize[3]) +
  labs(x = "Threshold Value", y = "Log St - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none") -> f5_lomax

nreps_mles_df %>%
  ggplot(., aes(y = log_St_hat2 - log(theta), group = thresh_tests, color = samp_n_tests)) +
  facet_wrap(~samp_n_tests, nrow = 1) +
  geom_jitter(aes(x = thresh_tests), alpha = 0.5) +
  geom_hline(aes(yintercept=  0), color = palette4_sampsize[1]) +
  scale_color_gradient(low = palette4_sampsize[2], high = palette4_sampsize[3]) +
  labs(x = "Threshold Value", y = "Log St.glm - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none") -> f5_glm

# nreps_mles_df %>%
#   ggplot(., aes(y = GP_theta_hat/theta, color = samp_n_tests)) +
#   facet_wrap(~samp_n_tests, nrow = 1) +
#   geom_jitter(aes(x = thresh_tests), alpha = 0.5) +
#   geom_hline(aes(yintercept=  1), color = palette4_sampsize[1]) +
#   scale_color_gradient(low = palette4_sampsize[2], high = palette4_sampsize[3]) +
#   labs(x = "Threshold Value", y = "GP theta hat /theta") +
#   scale_x_continuous(breaks = thresh_tests) +
#   theme_bw() +
#   theme(legend.position = "bottom") +
#   guides(color = guide_colorbar(barheight = 0.5, barwidth = 15, title = "Sample \n Size")) -> f5_evd

nreps_mles_df %>%
  filter(., GP_theta_hat != 0) %>% 
  ggplot(., aes(y = log(GP_theta_hat/theta), color = samp_n_tests)) +
  facet_wrap(~samp_n_tests, nrow = 1) +
  geom_jitter(aes(x = thresh_tests), alpha = 0.5) +
  geom_hline(aes(yintercept=  0), color = palette4_sampsize[1]) +
  scale_color_gradient(low = palette4_sampsize[2], high = palette4_sampsize[3]) +
  labs(x = "Threshold Value", y = "Log GP - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_colorbar(barheight = 0.5, barwidth = 15, title = "Sample \n Size")) -> f5_evd

pleg <- get_legend(f5_evd)

plot_grid(f5_lomax, f5_glm, f5_evd + theme(legend.position = "none"), pleg, nrow = 4, rel_heights = c(1,1,1, 0.3))
ggsave(paste0("Ch3_samplesize/", dir_scenario, "/Figure5.png"), width = 8, height = 8)

## Fig6 jitter raw values -----------------
# We see some strange trends, so we explore what happens when we consider the raw values and not the ratio

nreps_mles_df %>%
  ggplot(., aes(y = log_St_hat, color = samp_n_tests)) +
  facet_wrap(~samp_n_tests, nrow = 1) +
  geom_jitter(aes(x = thresh_tests), alpha = 0.5) +
  geom_point(data = thetas, aes(x = thresh_tests, y = log(theta)), color = "black", fill = "white", shape = 22, size = 1) +
  scale_color_gradient(low = palette4_sampsize[1], high = palette4_sampsize[3]) +
  labs(x = "Threshold Value", y = "Log St") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none") -> f6_lomax

nreps_mles_df %>%
  ggplot(., aes(y = log_St_hat2, group = thresh_tests, color = samp_n_tests)) +
  facet_wrap(~samp_n_tests, nrow = 1) +
  geom_jitter(aes(x = thresh_tests), alpha = 0.5) +
  geom_point(data = thetas, aes(x = thresh_tests, y = log(theta)), color = "black", fill = "white", shape = 22, size = 1) +
  scale_color_gradient(low = palette4_sampsize[1], high = palette4_sampsize[3]) +
  labs(x = "Threshold Value", y = "Log St.glm") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none") -> f6_glm


# nreps_mles_df %>%
#   ggplot(., aes(y = GP_theta_hat, color = samp_n_tests)) +
#   facet_wrap(~samp_n_tests, nrow = 1) +
#   geom_jitter(aes(x = thresh_tests), alpha = 0.5) +
#   geom_point(data = thetas, aes(x = thresh_tests, y = theta), color = "black", fill = "white", shape = 22, size = 1) +
#   scale_color_gradient(low = palette4_sampsize[1], high = palette4_sampsize[3]) +
#   labs(x = "Threshold Value", y = "GP theta hat") +
#   scale_x_continuous(breaks = thresh_tests) +
#   theme_bw() +
#   theme(legend.position = "bottom") +
#   guides(color = guide_colorbar(barheight = 0.5, barwidth = 15, title = "Sample \n Size"))  -> f6_evd

nreps_mles_df %>%
  filter(., GP_theta_hat != 0) %>% 
  ggplot(., aes(y = log(GP_theta_hat), color = samp_n_tests)) +
  facet_wrap(~samp_n_tests, nrow = 1) +
  geom_jitter(aes(x = thresh_tests), alpha = 0.5) +
  geom_point(data = thetas, aes(x = thresh_tests, y = log(theta)), color = "black", fill = "white", shape = 22, size = 1) +
  scale_color_gradient(low = palette4_sampsize[1], high = palette4_sampsize[3]) +
  labs(x = "Threshold Value", y = "GP theta hat") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_colorbar(barheight = 0.5, barwidth = 15, title = "Sample \n Size"))  -> f6_evd

pleg <- get_legend(f6_evd)

plot_grid(f6_lomax, f6_glm, f6_evd + theme(legend.position = "none"), pleg, nrow = 4, rel_heights = c(1,1,1, 0.3))
ggsave(paste0("Ch3_samplesize/", dir_scenario, "/Figure6.png"), width = 8, height = 8)

# MLES estimation error --------------------------
## Fig7 all param space ----------------------------

# This is an issue with parameter space, the estimation.
# Ok, so the survival function for the simple lomax starts to just give zero after a certain
# Range of sample sizes.
# Why?

ggplot(data = nreps_mles_df) +
  geom_point(aes(x = alpha_star, y = k_star, color = "Simple Lomax"), size = 3, alpha = 0.5) +
  geom_point(aes(x = alpha_glm, y = k_glm, color = "Lomax GLM"), size = 3, alpha = 0.5) +
  geom_point(aes(x = alpha_GP, y = k_GP, color = "GP estimate of Lomax"), size = 3, alpha = 0.5) +
  scale_y_log10() +
  scale_x_log10() +
  scale_color_discrete(name = "Model") +
  labs(x = "Log(Alpha_Star)", y = "Log(K_Star") +
  theme_bw() -> all_param_space
all_param_space
ggsave(paste0("Ch3_samplesize/", dir_scenario,"/Figure7.png"), width = 8, height = 8)

## Fig8 Param space by sampsize ----------------------------
sample_sizes <- samp_n_tests
PS <- list()
for(i in 1:length(sample_sizes)){
  ggplot(data = nreps_mles_df %>% filter(., samp_n_tests == sample_sizes[i])) +
    geom_point(aes(x = alpha_star, y = k_star, color = "Simple Lomax"), size = 3, alpha = 0.5) +
    geom_point(aes(x = alpha_glm, y = k_glm, color = "Lomax GLM"), size = 3, alpha = 0.5) +
    geom_point(aes(x = alpha_GP, y = k_GP, color = "GP estimate of Lomax"), size = 3, alpha = 0.5) +
    scale_y_log10() +
    scale_x_log10() +
    scale_color_discrete(name = "Model") +
    labs(x = "Log(Alpha_Star)", y = "Log(K_Star", title = paste("SampSize = ", sample_sizes[i])) +
    theme_bw() -> a
  PS[[i]] <- a + theme(legend.position = "none")
}


PS[[length(sample_sizes) + 1]] <- get_legend(a + theme(legend.position = "bottom"))

PS1 <- plot_grid(PS[[1]], PS[[2]], PS[[3]], PS[[4]], PS[[5]], PS[[6]], nrow = 2)
plot_grid(PS1, PS[[7]], nrow = 2, rel_heights = c(1, 0.1))
ggsave(paste0("Ch3_samplesize/", dir_scenario,"/Figure8.png"), width = 8, height = 8)



# GLM BIAS Correction ----------------------------------------------

## Fig9 Bias for glm model --------------------------------------

# Estimating bias as the sum of the difference between 

nreps_mles_df %>% 
  dplyr::select(., thresh_tests, samp_n_tests, log_St_hat2, theta) %>% 
  mutate(logtheta = log(theta),
         theta_error = log_St_hat2 - logtheta) %>% 
  group_by(thresh_tests, samp_n_tests) %>% 
  summarize(theta_bias = mean(theta_error)) %>% 
  ggplot(., aes(y = theta_bias, group = thresh_tests, color = samp_n_tests)) +
  geom_jitter(aes(x = thresh_tests), alpha = 0.9, size = 3) +
  geom_hline(aes(yintercept=  0), color = palette4_sampsize[2]) +
  scale_color_gradient(low = palette4_sampsize[1], high = palette4_sampsize[3]) +
  labs(x = "Threshold Value", y = "Bias Log St.glm") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_colorbar(barheight = 0.5, barwidth = 15, title = "Sample \n Size"))

ggsave(paste0("Ch3_samplesize/", dir_scenario,"/Figure9.png"), width = 6, height = 3.5)



## Nonparametric Bootstrap ------------------------

combo_tests <- data.frame(expand.grid(thresh_tests = thresh_tests, samp_n_tests = samp_n_tests))
nonpar_glm <- combo_tests

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
  nonpar_glm$theta_hat[i] <- ith_theta_hat
  
  # Nonparametric bootstrapping
  bth_df <- data.frame(n_B = 1:B)
  
  for(b in 1:B){
    # Sample with replacement from the given sample
    bth_samps <- data.frame(x_star = sample(ith_samps$x_star, ith_n, replace = TRUE))
    
    # Estimate params using glm Lomax
    bth_glm <- lomax.glm(formula=~1, my.dataf=bth_samps, response=bth_samps$x_star)
    bth_alpha_star <- bth_glm$alphas.hat[1]
    bth_k_star <- bth_glm$k.hat

    # Estimate the tail
    bth_df$log_theta_star[b] <- lomax.St(x = ith_thresh, alpha = bth_alpha_star, k = bth_k_star, log.scale=TRUE)
    bth_df$theta_star[b] <- lomax.St(x = ith_thresh, alpha = bth_alpha_star, k = bth_k_star, log.scale=FALSE)
    
  }
  
  nonpar_glm$log_theta_stars[i] <- mean(bth_df$log_theta_star)
  nonpar_glm$theta_stars[i] <- mean(bth_df$theta_star)
  
}

nonpar_glm %>%
  mutate(theta_bar = 2*theta_hat - log_theta_stars) %>% 
  right_join(., thetas[, c(1,3)]) %>% 
  mutate(theta = log(theta)) -> nonpar_glm


## Fig10 -------------------------------

nonpar_glm %>% 
  ggplot(., aes(y = theta_hat, color = samp_n_tests)) +
  geom_jitter(aes(x = thresh_tests), alpha = 0.7, size = 4, shape = 18) +
  geom_jitter(aes(x = thresh_tests, y = theta_bar), alpha = 0) +
  geom_point(aes(x = thresh_tests, y = theta), color = "black", size = 3) +
  scale_color_gradient(low = palette4_sampsize[1], high = palette4_sampsize[3]) +
  labs(x = "Threshold Value", y = "Theta hats") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) -> f10_top

nonpar_glm %>% 
  ggplot(., aes(y = theta_hat, color = samp_n_tests)) +
  # geom_jitter(aes(x = thresh_tests), alpha = 0.9, size = 6, shape = 5) +
  geom_jitter(aes(x = thresh_tests, y = theta_bar), size = 2, alpha = 0.8) +
  geom_point(aes(x = thresh_tests, y = theta), color = "black", size = 3) +
  scale_color_gradient(low = palette4_sampsize[1], high = palette4_sampsize[3]) +
  labs(x = "Threshold Value", y = "Theta bars") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) -> f10_bottom

f10_legend <- get_legend(f10_bottom + theme(legend.position = "bottom",
                                            axis.title.x = element_blank()) +
                           guides(color = guide_colorbar(barheight = 0.5, barwidth = 15, title = "Sample \n Size")))

plot_grid(f10_top, f10_bottom + theme(legend.position = "none"), f10_legend, nrow = 3, rel_heights = c(1,1,0.3))
ggsave(paste0("Ch3_samplesize/", dir_scenario,"/Figure10.png"), width = 8, height = 8)



### Multiple runs ---------------------------------------
nruns <- paste0("nrun_", 1:points_for_boxplots)
combo_tests <- data.frame(expand.grid(thresh_tests = thresh_tests, samp_n_tests = samp_n_tests, nruns = nruns))
nonpar_glm <- combo_tests

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
  nonpar_glm$theta_hat[i] <- ith_theta_hat
  
  # Nonparametric bootstrapping
  bth_df <- data.frame(n_B = 1:B)
  
  for(b in 1:B){
    # Sample with replacement from the given sample
    bth_samps <- data.frame(x_star = sample(ith_samps$x_star, ith_n, replace = TRUE))
    
    # Estimate params using glm Lomax
    bth_glm <- lomax.glm(formula=~1, my.dataf=bth_samps, response=bth_samps$x_star)
    bth_alpha_star <- bth_glm$alphas.hat[1]
    bth_k_star <- bth_glm$k.hat
    
    # Estimate the tail
    bth_df$log_theta_star[b] <- lomax.St(x = ith_thresh, alpha = bth_alpha_star, k = bth_k_star, log.scale=TRUE)
    bth_df$theta_star[b] <- lomax.St(x = ith_thresh, alpha = bth_alpha_star, k = bth_k_star, log.scale=FALSE)
    
  }
  
  nonpar_glm$log_theta_stars[i] <- mean(bth_df$log_theta_star)
  nonpar_glm$theta_stars[i] <- mean(bth_df$theta_star)
  
}

nonpar_glm %>%
  mutate(theta_bar = 2*theta_hat - log_theta_stars) %>% 
  right_join(., thetas[, c(1,3)]) %>% 
  mutate(theta = log(theta)) -> nonpar_glm

save(nonpar_glm, file = paste0("Ch3_samplesize/", dir_scenario, "/nonpar_glm.RData"))


### Fig11 ------

nonpar_glm %>% 
  mutate(hat = theta_hat - theta,
         bar = theta_bar - theta) %>% 
  dplyr::select(thresh_tests, samp_n_tests, hat, bar) %>% 
  pivot_longer(cols = c(hat, bar), names_to = "type") -> long_corr_df


long_corr_df %>% 
  group_by(thresh_tests, samp_n_tests, type) %>% 
  summarise(meanval = mean(value)) %>% 
  ungroup() %>% 
  right_join(., thetas[, c(1,3)]) %>% 
  mutate(thresh_tests = factor(thresh_tests)) -> means_df
  

long_corr_df %>% 
  mutate(thresh_tests = factor(thresh_tests)) %>% 
  ggplot(., aes(x = thresh_tests, y = value, fill = type)) +
  facet_wrap(~samp_n_tests, ncol = 2) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  scale_fill_manual(values = c("white", "grey")) +
  geom_hline(aes(yintercept=  0), color = palette4_sampsize[1]) +
  labs(x = "Threshold Value", y = "estimation - truth", title = "Nonparametric correction") +
  scale_y_continuous(limits = quantile(long_corr_df$value, c(0.1, 0.9))) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title = element_blank())

ggsave(paste0("Ch3_samplesize/", dir_scenario,"/Figure11.png"), width = 8, height = 8)




## Parametric Bootstrap ------------------------

combo_tests <- data.frame(expand.grid(thresh_tests = thresh_tests, samp_n_tests = samp_n_tests))
par_glm <- combo_tests

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
  par_glm$theta_hat[i] <- ith_theta_hat
  
  # Nonparametric bootstrapping
  bth_df <- data.frame(n_B = 1:B)
  
  for(b in 1:B){
    # Randrom draw from a Lomax with the estimated parameters
    bth_samps <- data.frame(x_star = rlomax(n = ith_n, alpha = ith_alpha_hat, k = ith_k_hat))
    
    # Estimate params using glm Lomax
    bth_glm <- lomax.glm(formula=~1, my.dataf=bth_samps, response=bth_samps$x_star)
    bth_alpha_star <- bth_glm$alphas.hat[1]
    bth_k_star <- bth_glm$k.hat
    
    # Estimate the tail
    bth_df$log_theta_star[b] <- lomax.St(x = ith_thresh, alpha = bth_alpha_star, k = bth_k_star, log.scale=TRUE)
    bth_df$theta_star[b] <- lomax.St(x = ith_thresh, alpha = bth_alpha_star, k = bth_k_star, log.scale=FALSE)
    
  }
  
  par_glm$log_theta_stars[i] <- mean(bth_df$log_theta_star)
  par_glm$theta_stars[i] <- mean(bth_df$theta_star)
  
}

par_glm %>%
  mutate(theta_bar = 2*theta_hat - log_theta_stars) %>% 
  right_join(., thetas[, c(1,3)]) %>% 
  mutate(theta = log(theta)) -> par_glm


## Fig12 -------------------------------

par_glm %>% 
  ggplot(., aes(y = theta_hat, color = samp_n_tests)) +
  geom_jitter(aes(x = thresh_tests), alpha = 0.7, size = 4, shape = 18) +
  geom_jitter(aes(x = thresh_tests, y = theta_bar), alpha = 0) +
  geom_point(aes(x = thresh_tests, y = theta), color = "black", size = 3) +
  scale_color_gradient(low = palette4_sampsize[1], high = palette4_sampsize[3]) +
  labs(x = "Threshold Value", y = "Theta hats") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) -> f12_top

par_glm %>% 
  ggplot(., aes(y = theta_hat, color = samp_n_tests)) +
  # geom_jitter(aes(x = thresh_tests), alpha = 0.9, size = 6, shape = 5) +
  geom_jitter(aes(x = thresh_tests, y = theta_bar), size = 2, alpha = 0.8) +
  geom_point(aes(x = thresh_tests, y = theta), color = "black", size = 3) +
  scale_color_gradient(low = palette4_sampsize[1], high = palette4_sampsize[3]) +
  labs(x = "Threshold Value", y = "Theta bars") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) -> f12_bottom

f12_legend <- get_legend(f12_bottom + theme(legend.position = "bottom",
                                            axis.title.x = element_blank()) +
                           guides(color = guide_colorbar(barheight = 0.5, barwidth = 15, title = "Sample \n Size")))

plot_grid(f12_top, f12_bottom + theme(legend.position = "none"), f12_legend, nrow = 3, rel_heights = c(1,1,0.3))
ggsave(paste0("Ch3_samplesize/", dir_scenario,"/Figure12.png"), width = 8, height = 8)



### Multiple runs ---------------------------------------
nruns <- paste0("nrun_", 1:points_for_boxplots)
combo_tests <- data.frame(expand.grid(thresh_tests = thresh_tests, samp_n_tests = samp_n_tests, nruns = nruns))
par_glm <- combo_tests

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
  par_glm$theta_hat[i] <- ith_theta_hat
  
  # Nonparametric bootstrapping
  bth_df <- data.frame(n_B = 1:B)
  
  for(b in 1:B){
    # Sample with replacement from the given sample
    bth_samps <- data.frame(x_star = sample(ith_samps$x_star, ith_n, replace = TRUE))
    
    # Estimate params using glm Lomax
    bth_glm <- lomax.glm(formula=~1, my.dataf=bth_samps, response=bth_samps$x_star)
    bth_alpha_star <- bth_glm$alphas.hat[1]
    bth_k_star <- bth_glm$k.hat
    
    # Estimate the tail
    bth_df$log_theta_star[b] <- lomax.St(x = ith_thresh, alpha = bth_alpha_star, k = bth_k_star, log.scale=TRUE)
    bth_df$theta_star[b] <- lomax.St(x = ith_thresh, alpha = bth_alpha_star, k = bth_k_star, log.scale=FALSE)
    
  }
  
  par_glm$log_theta_stars[i] <- mean(bth_df$log_theta_star)
  par_glm$theta_stars[i] <- mean(bth_df$theta_star)
  
}

par_glm %>%
  mutate(theta_bar = 2*theta_hat - log_theta_stars) %>% 
  right_join(., thetas[, c(1,3)]) %>% 
  mutate(theta = log(theta)) -> par_glm

save(par_glm, file = paste0("Ch3_samplesize/", dir_scenario, "/par_glm.RData"))


### Fig13 ------

par_glm %>% 
  mutate(hat = theta_hat - theta,
         bar = theta_bar - theta) %>% 
  dplyr::select(thresh_tests, samp_n_tests, hat, bar) %>% 
  pivot_longer(cols = c(hat, bar), names_to = "type") -> long_corr_df


long_corr_df %>% 
  group_by(thresh_tests, samp_n_tests, type) %>% 
  summarise(meanval = mean(value)) %>% 
  ungroup() %>% 
  right_join(., thetas[, c(1,3)]) %>% 
  mutate(thresh_tests = factor(thresh_tests)) -> means_df


long_corr_df %>% 
  mutate(thresh_tests = factor(thresh_tests)) %>% 
  ggplot(., aes(x = thresh_tests, y = value, fill = type)) +
  facet_wrap(~samp_n_tests, ncol = 2) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  scale_fill_manual(values = c("white", "grey")) +
  geom_hline(aes(yintercept=  0), color = palette4_sampsize[1]) +
  labs(x = "Threshold Value", y = "estimation - truth", title = "Parametric correction") +
  scale_y_continuous(limits = quantile(long_corr_df$value, c(0.1, 0.9))) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title = element_blank())

ggsave(paste0("Ch3_samplesize/", dir_scenario,"/Figure13.png"), width = 8, height = 8)


