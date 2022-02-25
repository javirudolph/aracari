#################################
## CH3 script
#################################


# Libraries -----------------------------------------------
# set.seed(20220201)

library(aracari)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(purrr)
library(fitdistrplus)

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

save(truth_df, file = paste0("Ch3_samplesize/simdata/mixturesamples", scenario, ".RData"))


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
truth_hist

# Density curve for the mixture based on the 50k samples.
truth_df %>%
  ggplot(., aes(x = x_samps)) +
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

ggsave(paste0("Ch3_samplesize/Figures/Figure1", scenario,".png"), width = 6, height = 5)


## Thresholds ----------------
# Thresholds based on quantiles: because when we have a sample, we still don't know the truth
# But we can know the quantiles of our sample.
# summary(truth_df$x_samps)
thresh_tests <- round(quantile(truth_df$x_samps, c(0.5, 0.75, 0.9, 0.99, 0.999, 0.9999), names = FALSE))

# thresh_tests <- c(50, 75, 100, 150, 250, 500, 750, 1000)
tibble(thresh_tests) %>%
  mutate(n_overthresh = map_dbl(1:length(thresh_tests),
                                function(y) length(which(truth_df$x_samps >= thresh_tests[y]))),
         theta = n_overthresh/tru_n) -> thetas
thetas



# Lomax functions ---------------------------------
lomax.pdf <- function(x,alpha,k, log.scale=FALSE){

  if(log.scale==FALSE){out <- (k/(alpha+x))*(alpha/(alpha+x))^k
  }else{
    out <- log(k) + k*log(alpha) - (k+1)*log(alpha+x)
  }

  return(out)
}


lomax.cdf <- function(x,alpha,k){

  return(1-(alpha/(alpha+x))^k)

}

lomax.St <- function(x,alpha,k, log.scale=FALSE){

  if(log.scale==FALSE){out <- (alpha/(alpha+x))^k
  }else{
    out <- k*log(alpha) - k*log(alpha+x)
  }

  return(out)

}

# Simple negative log-likelihood
nllike.simp <- function(guess=c(1.5,1.5), Y=Y){

  parms         <- exp(guess)
  alpha         <- parms[1]
  k             <- parms[2]
  n             <- length(Y)
  lnft.yis     <- lomax.pdf(x=Y, alpha=alpha,k=k,log.scale=TRUE)
  lnL           <- sum(lnft.yis)
  nll           <- -lnL
  return(nll)
}


ft.nllike2 <- function(guess, designmat=designmat,Y=Y){

  # For the glm idea to work we want ln(E[Lomax]) = ln(alpha/(k-1)) = XB
  # That is the proper link function.  Then,
  # ln alpha = XB + ln(k-1), or
  # alpha = exp(XB + ln(k-1)) and it follows that
  # k = exp(ln(k-1)) +1.  To avoid problems with potential log(negative number)
  # we optimize 'k-1', not 'k'

  nbetasp1      <- length(guess)
  lnkm1         <- guess[1]
  Xbeta         <- designmat%*%guess[2:nbetasp1]
  ln.alphas     <- Xbeta + lnkm1
  alphas        <- exp(ln.alphas)
  k             <- exp(lnkm1)+1
  n             <- length(Y)

  #sumlogapy     <- sum(log(alphas+Y))
  #k.hat         <- n/(sumlogapy - sum(Xbeta))
  lnft.yis     <- lomax.pdf(x=Y, alpha=alphas,k=k,log.scale=TRUE)
  lnL           <- sum(lnft.yis)
  nll           <- -lnL

  return(nll)
}

lomax.glm <- function(formula=~1, my.dataf, response){

  Y           <- response
  n           <- length(Y)
  designmat   <- model.matrix(formula, data=my.dataf)
  nbetas      <- ncol(designmat)
  init.betas  <- c(4,rep(5,nbetas))

  opt.out <- optim(par=init.betas, fn=ft.nllike2, method = "Nelder-Mead",
                   designmat=designmat, Y=Y)

  mles          <- opt.out$par
  nll.hat       <- opt.out$value
  BIC.mod       <- 2*nll.hat + (nbetas+1)*log(length(Y))
  Xbeta.hat     <- designmat%*%mles[-1]
  lnkm1.hat     <- mles[1]
  k.hat         <- exp(lnkm1.hat)+1
  alphas.hat    <- exp(Xbeta.hat +lnkm1.hat)

  out.list <- list(opt.out = opt.out, designmat=designmat,Y=Y, mles=mles, nll.hat=nll.hat, BIC.mod = BIC.mod,
                   alphas.hat=alphas.hat, k.hat=k.hat,data=my.dataf)

  return(out.list)

}


# Estimate tails once -----------------------------------------
library(extRemes)

samp_n_tests <- c(80, 200, 500, 800, 1000, 1600)


mles_df <- data.frame(expand.grid(thresh_tests = thresh_tests, samp_n_tests = samp_n_tests))

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
  mles_df$log_GP[i] <- log(pextRemes(ith_evd, ith_thresh, lower.tail = FALSE))


}

mles_df %>%
  right_join(., thetas[, c(1,3)]) -> mles_df

### VIZ ---------------------------------------------

### ISSUE with parameter space

mles_df %>%
  ggplot() +
  geom_point(aes(x = alpha_star, y = k_star), color = "red", size = 3) +
  geom_point(aes(x = alpha_glm, y = k_glm), color = "blue", size = 3) +
  geom_point(aes(x = alpha_GP, y = k_GP), color = "green", size = 3) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw()


# Q: If we focus ONLY on the simple Lomax
# - How close to the truth is the estimate of the tai?

mles_df %>%
  ggplot(., aes(x = thresh_tests, y = log_St_hat - log(theta))) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept=  0), color = "red") +
  labs(x = "Threshold Value", y = "Log St - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() -> p1
p1
# As the threshold increases, the variation in the estimation of the threshold increases
# Makes sense. It's harder to get precise estimates for very rare events

# Q: How does it look like for all three methods?
mles_df %>%
  ggplot(., aes(x = thresh_tests, y = log_St_hat2 - log(theta))) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept=  0), color = "red") +
  labs(x = "Threshold Value", y = "Log St.glm - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() -> p2
mles_df %>%
  ggplot(., aes(x = thresh_tests, y = log_GP - log(theta))) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept=  0), color = "red") +
  labs(x = "Threshold Value", y = "Log GP - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() -> p3

plot_grid(p1, p2, p3, nrow = 3)
ggsave(paste0("Ch3_samplesize/Figures/Figure2", scenario,".png"), width = 8, height = 8)



# Q: How does this change with sample size?
# I would expect that by increasing the sample size, we get closer to the true estimator

mypal <- c("#95190C", "#610345", "#E3B505")

mles_df %>%
  ggplot(., aes(x = thresh_tests, y = log_St_hat - log(theta))) +
  geom_point(size = 3, alpha = 0.7, aes(color = samp_n_tests)) +
  geom_hline(aes(yintercept=  0), color = mypal[1]) +
  scale_color_gradient(low = mypal[2], high = mypal[3]) +
  labs(x = "Threshold Value", y = "Log St - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none") -> p1
p1
mles_df %>%
  ggplot(., aes(x = thresh_tests, y = log_St_hat2 - log(theta))) +
  geom_point(size = 3, alpha = 0.7, aes(color = samp_n_tests)) +
  geom_hline(aes(yintercept=  0), color = mypal[1]) +
  scale_color_gradient(low = mypal[2], high = mypal[3]) +
  labs(x = "Threshold Value", y = "Log St.glm - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none") -> p2
mles_df %>%
  ggplot(., aes(x = thresh_tests, y = log_GP - log(theta))) +
  geom_point(size = 3, alpha = 0.7, aes(color = samp_n_tests)) +
  geom_hline(aes(yintercept=  0), color = mypal[1]) +
  scale_color_gradient(low = mypal[2], high = mypal[3]) +
  labs(x = "Threshold Value", y = "Log GP - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_colorbar(barheight = 0.5, barwidth = 15, title = "Sample \n Size")) -> p3

pleg <- get_legend(p3)

plot_grid(p1, p2, p3 + theme(legend.position = "none"), pleg, nrow = 4, rel_heights = c(1,1,1, 0.3))
ggsave(paste0("Ch3_samplesize/Figures/Figure3", scenario,".png"), width = 8, height = 8)

### conclusion ------------------------------------------

# So, the GP tends to go way off
# The less biased if you will is the regular lomax and we see clearly that sample size matters.


# Make some boxplots ----------------------------------------------------

bxplt_data_fx <- function(nreps){
  out <- data.frame()

  for(j in 1:nreps){
    samp_n_tests <- c(80, 200, 500, 800, 1000, 1600)
    mles_df <- data.frame(expand.grid(thresh_tests = thresh_tests, samp_n_tests = samp_n_tests))
    mles_df$nrep <- j

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
      mles_df$alpha_star[i] <- alpha_star
      mles_df$k_star[i] <- k_star
      # Estimate the tail
      log_St_star <- lomax.St(x = ith_thresh,alpha = alpha_star,k = k_star,log.scale=TRUE)
      mles_df$log_St_hat[i] <- log_St_star

      # Check with GP fit
      ith_evd <- fevd(ith_samps$x_star, threshold = 0, type = "GP")
      evd_mles <- summary(ith_evd)$par
      mles_df$scale[i] <- evd_mles[1]
      mles_df$shape[i] <- evd_mles[2]
      mles_df$log_GP[i] <- log(pextRemes(ith_evd, ith_thresh, lower.tail = FALSE))

      # Check with Lomax GLM
      glm_out <- lomax.glm(formula=~1, my.dataf=ith_samps, response=ith_samps$x_star)
      alpha.2 <- glm_out$alphas.hat[1]
      k.2     <- glm_out$k.hat
      log.st.star2 <- lomax.St(x=ith_thresh,alpha=alpha.2,k=k.2,log.scale=TRUE)
      mles_df$log_St_hat2[i] <- log.st.star2

    }

    mles_df %>%
      right_join(., thetas[, c(1,3)]) -> mles_df
    out <- rbind.data.frame(out, mles_df)
  }

  return(out)

}

nreps_mles_df <- bxplt_data_fx(30)

## bxplt ratio ----------------------------
nreps_mles_df %>%
  #mutate(thresh_tests = factor(thresh_tests)) %>%
  ggplot(., aes(y = log_St_hat - log(theta), group = thresh_tests, color = samp_n_tests)) +
  facet_wrap(~samp_n_tests, nrow = 1) +
  geom_boxplot() +
  geom_hline(aes(yintercept=  0), color = mypal[1]) +
  scale_color_gradient(low = mypal[2], high = mypal[3]) +
  labs(x = "Threshold Value", y = "Log St - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none") -> p1

nreps_mles_df %>%
  #mutate(thresh_tests = factor(thresh_tests)) %>%
  ggplot(., aes(y = log_St_hat2 - log(theta), group = thresh_tests, color = samp_n_tests)) +
  facet_wrap(~samp_n_tests, nrow = 1) +
  geom_boxplot() +
  geom_hline(aes(yintercept=  0), color = mypal[1]) +
  scale_color_gradient(low = mypal[2], high = mypal[3]) +
  labs(x = "Threshold Value", y = "Log St.glm - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none") -> p2

nreps_mles_df %>%
  #mutate(thresh_tests = factor(thresh_tests)) %>%
  ggplot(., aes(y = log_GP - log(theta), group = thresh_tests, color = samp_n_tests)) +
  facet_wrap(~samp_n_tests, nrow = 1) +
  geom_boxplot() +
  geom_hline(aes(yintercept=  0), color = mypal[1]) +
  scale_color_gradient(low = mypal[2], high = mypal[3]) +
  labs(x = "Threshold Value", y = "Log GP - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_colorbar(barheight = 0.5, barwidth = 15, title = "Sample \n Size"))  -> p3

pleg <- get_legend(p3)

plot_grid(p1, p2, p3 + theme(legend.position = "none"), pleg, nrow = 4, rel_heights = c(1,1,1, 0.3))
ggsave(paste0("Ch3_samplesize/Figures/Figure4", scenario,".png"), width = 8, height = 8)

## jitter ratio -------------------------------------
nreps_mles_df %>%
  #mutate(thresh_tests = factor(thresh_tests)) %>%
  ggplot(., aes(y = log_St_hat - log(theta), color = samp_n_tests)) +
  facet_wrap(~samp_n_tests, nrow = 1) +
  geom_jitter(aes(x = thresh_tests)) +
  geom_hline(aes(yintercept=  0), color = mypal[1]) +
  scale_color_gradient(low = mypal[2], high = mypal[3]) +
  labs(x = "Threshold Value", y = "Log St - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none") -> p1

nreps_mles_df %>%
  #mutate(thresh_tests = factor(thresh_tests)) %>%
  ggplot(., aes(y = log_St_hat2 - log(theta), group = thresh_tests, color = samp_n_tests)) +
  facet_wrap(~samp_n_tests, nrow = 1) +
  geom_jitter(aes(x = thresh_tests)) +
  geom_hline(aes(yintercept=  0), color = mypal[1]) +
  scale_color_gradient(low = mypal[2], high = mypal[3]) +
  labs(x = "Threshold Value", y = "Log St.glm - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none") -> p2

nreps_mles_df %>%
  #mutate(thresh_tests = factor(thresh_tests)) %>%
  ggplot(., aes(y = log_GP - log(theta), color = samp_n_tests)) +
  facet_wrap(~samp_n_tests, nrow = 1) +
  geom_jitter(aes(x = thresh_tests)) +
  geom_hline(aes(yintercept=  0), color = mypal[1]) +
  scale_color_gradient(low = mypal[2], high = mypal[3]) +
  labs(x = "Threshold Value", y = "Log GP - Log theta") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_colorbar(barheight = 0.5, barwidth = 15, title = "Sample \n Size"))  -> p3

pleg <- get_legend(p3)

plot_grid(p1, p2, p3 + theme(legend.position = "none"), pleg, nrow = 4, rel_heights = c(1,1,1, 0.3))
ggsave(paste0("Ch3_samplesize/Figures/Figure5", scenario,".png"), width = 8, height = 8)

## jitter raw values -----------------
nreps_mles_df %>%
  #mutate(thresh_tests = factor(thresh_tests)) %>%
  ggplot(., aes(y = log_St_hat, color = samp_n_tests)) +
  facet_wrap(~samp_n_tests, nrow = 1) +
  geom_jitter(aes(x = thresh_tests)) +
  geom_hline(aes(yintercept=  0), color = mypal[1]) +
  scale_color_gradient(low = mypal[2], high = mypal[3]) +
  labs(x = "Threshold Value", y = "Log St") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none") -> p1

nreps_mles_df %>%
  #mutate(thresh_tests = factor(thresh_tests)) %>%
  ggplot(., aes(y = log_St_hat2, group = thresh_tests, color = samp_n_tests)) +
  facet_wrap(~samp_n_tests, nrow = 1) +
  geom_jitter(aes(x = thresh_tests)) +
  geom_hline(aes(yintercept=  0), color = mypal[1]) +
  scale_color_gradient(low = mypal[2], high = mypal[3]) +
  labs(x = "Threshold Value", y = "Log St.glm") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "none") -> p2

nreps_mles_df %>%
  #mutate(thresh_tests = factor(thresh_tests)) %>%
  ggplot(., aes(y = log_GP, color = samp_n_tests)) +
  facet_wrap(~samp_n_tests, nrow = 1) +
  geom_jitter(aes(x = thresh_tests)) +
  geom_hline(aes(yintercept=  0), color = mypal[1]) +
  scale_color_gradient(low = mypal[2], high = mypal[3]) +
  labs(x = "Threshold Value", y = "Log GP") +
  scale_x_continuous(breaks = thresh_tests) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_colorbar(barheight = 0.5, barwidth = 15, title = "Sample \n Size"))  -> p3

pleg <- get_legend(p3)

plot_grid(p1, p2, p3 + theme(legend.position = "none"), pleg, nrow = 4, rel_heights = c(1,1,1, 0.3))
ggsave(paste0("Ch3_samplesize/Figures/Figure6", scenario,".png"), width = 8, height = 8)












