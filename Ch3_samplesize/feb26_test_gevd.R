# Libraries -----------------------------------------------
set.seed(22227)

library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(purrr)
library(fitdistrplus)
library(extRemes)

source("Ch3_samplesize/PoncianoDesktopFeb26/Ch3_samplesize/Ch3_functions.R")

# ORIGINAL SCENARIO --------------------------------------
dir_scenario <- "orig_feb26_tests"

if(dir.exists(paste0("Ch3_samplesize/", dir_scenario)) == FALSE){
  dir.create(paste0("Ch3_samplesize/", dir_scenario))}

points_for_boxplots <- 10
B <- 10

desired_means <- c(28, 32, 40, 50)
desired_sds <- c(49.7, 39.9, 33.3, 31.1)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)

pis <- c( 0.1, 0.2, 0.3, 0.4)
tru_n <- 50000
w_samp_sizes <- pis*tru_n
purrrsampls <- tibble(gID = c(1:length(w_samp_sizes)), pars) %>%
  mutate(x_samps = purrr::map(gID, function(y) rlnorm(w_samp_sizes[y], pars$meanlog[y], pars$sdlog[y])))

purrrsampls %>%
  dplyr::select(., gID, x_samps) %>%
  unnest(cols = c(x_samps)) -> truth_df


# TEST VALUES --------------------------------------
#Thresholds

thresh_tests <- c(50, 100, 150, 250, 500, 750, 1000)
tibble(thresh_tests) %>%
  mutate(n_overthresh = map_dbl(1:length(thresh_tests),
                                function(y) length(which(truth_df$x_samps >= thresh_tests[y]))),
         theta = n_overthresh/tru_n) -> thetas


#Samlesizes
samp_n_tests <- c(80, 200, 500, 800, 1000, 1600)

all_test_combos <- data.frame(expand.grid(thresh_tests = thresh_tests, samp_n_tests = samp_n_tests))



# RERUN MLES with threshold -------------------------------------

# Load the nreps from the ponciano desktop folder
# Run nreps with a threshold for GP (nreps is 50)

out <- data.frame()

for(j in 1:50){
  samp_n_tests <- c(80, 200, 500, 800, 1000, 1600)
  mles_df <- data.frame(expand.grid(thresh_tests = thresh_tests, samp_n_tests = samp_n_tests))
  mles_df$nrep <- j

  for(i in 1:nrow(mles_df)){
    ith_n <- mles_df$samp_n_tests[i]
    ith_thresh <- mles_df$thresh_tests[i]
    ith_samps <- data.frame(x_star = sample(truth_df$x_samps, ith_n))
    ith_tail <- length(which(ith_samps$x_star >= ith_thresh))
    mles_df$tail_p[i] <- ith_tail/ith_n
    ith_quant <- round(quantile(ith_samps$x_star, 0.5))


    # Check with GP fit
    ith_evd <- fevd(ith_samps$x_star, threshold = ith_quant, type = "GP")
    gp_scale <- summary(ith_evd)$par[1]
    gp_shape <- summary(ith_evd)$par[2]

    # GP tail
    # mles_df$log_GP[i] <- log(pextRemes(ith_evd, ith_thresh, lower.tail = FALSE))
    mles_df$GP_t[i] <- pextRemes(ith_evd, ith_thresh, lower.tail = FALSE)
  }

  mles_df %>%
    right_join(., thetas[, c(1,3)]) -> mles_df
  out <- rbind.data.frame(out, mles_df)
}


# Compare them -------------------------------

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

out %>%
  filter(., GP_t != 0) %>%
  ggplot(., aes(y = log(GP_t/theta), group = thresh_tests, color = samp_n_tests)) +
  facet_wrap(~samp_n_tests, nrow = 1) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(aes(yintercept=  0), color = palette4_sampsize[1]) +
  scale_color_gradient(low = palette4_sampsize[2], high = palette4_sampsize[3]) +
  labs(x = "Threshold Value", y = "Log GP - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none")  -> f5_evd


plot_grid(f4_lomax, f4_glm, f4_evd + theme(legend.position = "none"), f5_evd, pleg, nrow = 5, rel_heights = c(1,1,1,1, 0.3))

# FIND THE THRESHOLD -------------------

ith_n <- samp_n_tests[5]

ith_samps <- data.frame(x_star = sample(truth_df$x_samps, ith_n))
mrlplot(ith_samps$x_star)

threshrange.plot(ith_samps$x_star, r = c(0, 500), type = "GP")

fit <- fevd(ith_samps$x_star, threshold = 50, type = "GP")
plot(fit)

# How do we use these tools from the pareto stuff with the lomax distribution
# How do we inform the model to correct the bias in the tail estimate?


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
