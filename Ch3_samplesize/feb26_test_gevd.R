# Libraries -----------------------------------------------
set.seed(22227)

library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(purrr)
library(fitdistrplus)
library(extRemes)

# Lnorm parms functions --------------------------------------

desired_mean_sd <- function(mu_x, sd_x){

  sigsq <- sd_x^2

  mu <- log(mu_x^2/(sqrt(mu_x^2+sigsq)))
  sigma_sq <- log(1+(sigsq/mu_x^2))
  sigma <- sqrt(sigma_sq)

  return(data.frame(meanlog = mu, sigma_sq = sigma_sq, sdlog = sigma))
}

# Now, for a Lnorm(meanlog, sdlog), get mean and var

lnorm_mean_var <- function(mean_log, sd_log){

  lnorm_mean <- exp(mean_log + ((sd_log^2)/2))
  lnorm_var  <- (exp(sd_log^2)-1)*exp(2*mean_log+sd_log^2)
  lnorm_sd   <- sqrt(lnorm_var)

  return(data.frame(lnorm_mean, lnorm_var, lnorm_sd))
}


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


## Nonparametric bootstrap using Generalized Pareto

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
