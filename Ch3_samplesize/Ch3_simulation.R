###
## DATA GENERATION CH3 ------------------------------------
###

# Libraries -----------------------------------------------
#set.seed(20220201)

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
## Mixture components --------------------------------------
# Since movement lengths are only positive, we use a lognormal distribution to describe them


# scenario <- "_original"
# desired_means <- c(28, 32, 40, 50)
# desired_sds <- c(49.7, 39.9, 33.3, 31.1)
#
#
# pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)

# Just checking that it makes sense now
lnorm_mean_var(pars$meanlog, pars$sdlog)
pars


mus <- pars$meanlog
sigsqs <- pars$sdlog
dens_cols <- c("#264653", "#2a9d8f", "#f4a261", "#e76f51")

## Mixture proportions ----------------------------------------
# We assume four categories of individuals with increasing movement.
pis <- c( 0.1, 0.2, 0.3, 0.4)
sum(pis)
hist(pis)

### Plot the mixture components -------------------------------

# Using I function I wrote that I use to make a lot of these density curves for a plot
# can be found in the functions script.
lnorm_densities <- lnorm_densities_fx(pars$meanlog, pars$sdlog, dens_cols)

# But choosing to use the weighted curves instead
weighted_densities <- purrr::map(1:4, function(y) stat_function(fun = dlnorm,
                                                                args = list(meanlog = pars$meanlog[y], sdlog = pars$sdlog[y]),
                                                                color = dens_cols[y], size=pis[y]*5, alpha = 0.8))

ggplot() +
  #lnorm_densities +
  weighted_densities +
  theme_minimal() +
  labs(y = "Density") +
  lims(x = c(0, 150)) -> densities_plot
densities_plot


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


top_row <- plot_grid(densities_plot, means_sd_plot, rel_widths = c(2,1))
top_row

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

save(simplsamps, file = paste0("Ch3_samplesize/simdata/mixturesamples", scenario, ".RData"))

# Extract the tail (the highest values) to visualize later in the histogram since they are too few to show in the bins
data.tail <- data.frame(values = simplsamps$data, y = 100) %>% arrange(desc(values)) %>% filter(values >=100)

## Visualization and save -----------------------------------------------------

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

plot_grid(top_row, bottom_row,nrow = 2)

ggsave(paste0("Ch3_samplesize/Figures/Figure1", scenario,".png"))


###
## Probability of tails
###

thresh.vals <- c(100, 250,500,750,1000)

weighted.cdfs <- NULL
for(y in 1:4){
  a <- plnorm(thresh.vals, meanlog = pars$meanlog[y], sdlog = pars$sdlog[y], lower.tail = FALSE) * pis[y]
  weighted.cdfs <- rbind(weighted.cdfs, a)
}
w.cdfs <- colSums(weighted.cdfs)


tibble(thresh.vals) %>%
  mutate(w.cdfs = w.cdfs,
         w.cdfs = round(w.cdfs, 5),
         samp.n = map_dbl(1:length(thresh.vals), function(y) length(which(simplsamps$data >= thresh.vals[y]))),
         samp.p = signif(samp.n/samp.size, 3)) -> tru.cdfs



###
## FITTING EVD ----------------------------------
##

library(extRemes)

samp.sizes <- c(80, 200, 500, 800, 1000, 1600)

fevd.mles <- data.frame(samp.size = rep(samp.sizes, length(thresh.vals)), threshold = thresh.vals, scale = 0, shape = 0, nllh = 0)

for(i in 1:nrow(fevd.mles)){
  ith.n        <- fevd.mles$samp.size[i]
  ith.samples  <- data.frame(x = sample(simplsamps$data, ith.n))

  # So, setting the threshold to 0 so we can compare to Lomax.

  ith.fit      <- fevd(ith.samples$x, threshold = 0, type = "GP")

  mles          <- summary(ith.fit)$par
  nll.hat       <- summary(ith.fit)$nllh
  BIC.mod       <- summary(ith.fit)$BIC

  fevd.mles$scale[i] <- mles[1]
  fevd.mles$shape[i] <- mles[2]
  fevd.mles$nllh[i]  <- nll.hat

  ith.thresh   <- fevd.mles$threshold[i]
  fevd.mles$gp.tail[i] <- pextRemes(ith.fit, ith.thresh, lower.tail = FALSE)
}

# Plot parameter space ---------------------------
fevd.mles %>%
  mutate(samp.size = factor(samp.size),
         threshold = factor(threshold)) %>%
  ggplot(., aes(x = shape, y = scale, color = threshold
                #, shape = samp.size
                )) +
  facet_wrap(~samp.size) +
  geom_point() +
  theme_bw()

# Plot tail estimates -----------------------------

fevd.mles$tru.tail <- tru.cdfs$w.cdfs

fevd.mles %>%
  mutate(threshold = factor(threshold),
         samp.size = factor(samp.size),
         tail.ratio = gp.tail/tru.tail) %>%
  ggplot(., aes(x = samp.size, y = tail.ratio, color = threshold)) +
  geom_point() +
  labs(y = "GP tail.ratio") +
  theme_bw() +
  theme(legend.position = "none") -> a
a

# This is the conversion to a Lomax using the GP parameters.
# Then calculate the tail using the lomax function
# And getting the ratio
# We get the same plot as above


fevd.mles %>%
  mutate(k = 1/shape,
         alpha = scale*k,
         lomax.tail = lomax.st(threshold, alpha, k),
         tail.ratio = lomax.tail/tru.tail,
         threshold = factor(threshold),
         samp.size = factor(samp.size)) %>%
  ggplot(., aes(x = samp.size, y = tail.ratio, color = threshold)) +
  geom_point() +
  labs(y = "est Lomax tail.ratio") +
  theme_bw() -> b
b

plot_grid(a,b, rel_widths = c(0.8, 1))
ggsave(paste0("Ch3_samplesize/Figures/tail_ratio", scenario, ".png"))



## MC Samples -----------------------------------------------

nreps <- 100
gp.mles.reps <- data.frame(NULL)

for(j in 1:nreps){

  ith.mle.df <- data.frame(samp.size = rep(samp.sizes, length(thresh.vals)), threshold = thresh.vals,
                           scale = 0, shape = 0, nllh = 0, gp.tail = 0, rep = 0)

  for(i in 1:nrow(ith.mle.df)){
    ith.n        <- ith.mle.df$samp.size[i]
    ith.samples  <- data.frame(x = sample(simplsamps$data, ith.n))

    # So, setting the threshold to 0 so we can compare to Lomax.

    ith.fit      <- fevd(ith.samples$x, threshold = 0, type = "GP")

    mles          <- summary(ith.fit)$par
    nll.hat       <- summary(ith.fit)$nllh
    BIC.mod       <- summary(ith.fit)$BIC

    ith.mle.df$scale[i] <- mles[1]
    ith.mle.df$shape[i] <- mles[2]
    ith.mle.df$nllh[i]  <- nll.hat
    ith.mle.df$rep[i] <- paste0("rep", j)

    ith.thresh   <- ith.mle.df$threshold[i]
    ith.mle.df$gp.tail[i] <- pextRemes(ith.fit, ith.thresh, lower.tail = FALSE)

  }

  gp.mles.reps <- rbind.data.frame(gp.mles.reps, ith.mle.df)

}



gp.mles.reps$tru.tail <- tru.cdfs$w.cdfs


gp.mles.reps %>%
  mutate(k = 1/shape,
         alpha = scale*k,
         lomax.tail = lomax.st(threshold, alpha, k),
         lm.tail.ratio = lomax.tail/tru.tail,
         gp.tail.ratio = gp.tail/tru.tail,
         threshold = factor(threshold),
         samp.size = factor(samp.size)) -> gp.mles.reps


save(gp.mles.reps, file = paste0("Ch3_samplesize/simdata/GPmles", scenario, ".RData"))


# Boxplots -------------------

gp.mles.reps %>%
  ggplot(., aes(x = samp.size, y = gp.tail.ratio, color = threshold)) +
  geom_boxplot() +
  labs(y = "GP tail.ratio") +
  lims(y = c(-0.1, 5)) +
  theme_bw()

ggsave(paste0("Ch3_samplesize/Figures/GPtail_boxplt", scenario, ".png"))

# Mean and Standard Error plots -------------

head(gp.mles.reps)

gp.mles.reps %>%
  mutate(sqrd_diff = (gp.tail-tru.tail)^2) %>%
  group_by(samp.size, threshold) %>%
  summarise(mean.ratio = mean(gp.tail.ratio),
            ste.gp.ratio = sd(gp.tail.ratio)/sqrt(length(gp.tail.ratio)),
            mse = (1/length(gp.tail.ratio))*sum(sqrd_diff)) -> tail.ratio.means


tail.ratio.means %>%
  mutate(lo = mean.ratio - ste.gp.ratio,
         hi = mean.ratio + ste.gp.ratio) %>%
  ggplot(., aes(color = threshold,
                y = mean.ratio,
                # y = mse,
                x = samp.size)) +
  geom_jitter(size = 2, width = 0.2, alpha = 0.8) +
  #scale_y_log10(name = "Log(mean.ratio)") +
  theme_bw() -> a
a


tail.ratio.means %>%
  mutate(lo = mean.ratio - ste.gp.ratio,
         hi = mean.ratio + ste.gp.ratio) %>%
  ggplot(., aes(x = threshold,
                y = mean.ratio,
                # y = mse,
                color = samp.size)) +
  # facet_wrap(~samp.size) +
  # geom_point() +
  # geom_errorbar(aes(ymin = lo, ymax = hi), width = 0) +
  geom_jitter(size = 2, width = 0.2, alpha = 0.8) +
  #scale_y_log10(name = "Log(mean.ratio)") +
  theme_bw() -> b
b

top <- plot_grid(a,b, ncol=2)

bottom <- plot_grid(a + scale_y_log10(name = "Log(mean.ratio)"), b + scale_y_log10(name = "Log(mean.ratio)"))

plot_grid(top, bottom, nrow =2)

ggsave(paste0("Ch3_samplesize/Figures/GPtail_mean", scenario, ".png"))





# UNBIASED ESTIMATOR ------------------------------------------------

## Testing ----------------------------------------------------------

ith.n <- 1000 # the sample size
ith.df <- data.frame(x = sample(simplsamps$data, ith.n))

# Estimate the shape and scale parameters
ith.evd <- fevd(x = ith.df$x, type = "GP", threshold = 0)
ith.scale <- ith.evd$results$par[1]
ith.shape <- ith.evd$results$par[2]


# Using the same threshold values as above
# These are the true thetas for the different thresholds
ith.thetas <- pevd(thresh.vals, loc = 0, scale = ith.scale, shape= ith.shape, lower.tail = FALSE)


# with those parameters, simulate B samples
# Then for each sample get the mles

B <- 100
jth_df <- data.frame(jscale = 0, jshape = 0)

for(j in 1:B){


  jth.draw <- revd(n = ith.n, loc = 0, scale = ith.scale, shape = ith.shape)
  jth.fit <- fevd(x = jth.draw, type = "GP", threshold = 0)
  jth.scale <- jth.fit$results$par[1]
  jth.shape <- jth.fit$results$par[2]

  jth_df[j,] <- c(jth.scale, jth.shape)

}

# Using the new mles, estimate the thetas for each threshold
# Calculate the bias and the unbiased estimator for each threshold

bias_df <- NULL

for(k in 1:5){
  kth_thresh <- thresh.vals[k]

  est_theta <- jth_df %>%
    mutate(kth_theta = map_dbl(1:B, function(y) pevd(kth_thresh, loc = 0,
                                                 scale = jscale[y], shape = jshape[y], lower.tail = FALSE)))
  est_theta$thresh <- kth_thresh
  est_theta$tru_theta <- ith.thetas[k]

  bias_df <- rbind.data.frame(bias_df, est_theta)
}

# Ok, so this works. Make a function to do over sample sizes

# Bayas FXs ----------------------------

bias_fx <- function(B, theta_i, tru.theta){
  (1/B)*sum(theta_i - tru.theta)
}

bayas_fx <- function(samplesize = samplesize, B=100){
  ith.n <- samplesize # the sample size
  ith.df <- data.frame(x = sample(simplsamps$data, ith.n))

  # Estimate the shape and scale parameters
  ith.evd <- fevd(x = ith.df$x, type = "GP", threshold = 0)
  ith.scale <- ith.evd$results$par[1]
  ith.shape <- ith.evd$results$par[2]


  # Using the same threshold values as above
  # These are the true thetas for the different thresholds
  ith.thetas <- pevd(thresh.vals, loc = 0, scale = ith.scale, shape= ith.shape, lower.tail = FALSE)


  # with those parameters, simulate B samples
  # Then for each sample get the mles

  jth_df <- data.frame(jscale = 0, jshape = 0)

  for(j in 1:B){


    jth.draw <- revd(n = ith.n, loc = 0, scale = ith.scale, shape = ith.shape)
    jth.fit <- fevd(x = jth.draw, type = "GP", threshold = 0)
    jth.scale <- jth.fit$results$par[1]
    jth.shape <- jth.fit$results$par[2]

    jth_df[j,] <- c(jth.scale, jth.shape)

  }

  # Using the new mles, estimate the thetas for each threshold
  # Calculate the bias and the unbiased estimator for each threshold

  bias_df <- NULL

  for(k in 1:5){
    kth_thresh <- thresh.vals[k]

    est_theta <- jth_df %>%
      mutate(kth_theta = map_dbl(1:B, function(y) pevd(kth_thresh, loc = 0,
                                                       scale = jscale[y], shape = jshape[y], lower.tail = FALSE)))
    est_theta$thresh <- kth_thresh
    est_theta$tru_theta <- ith.thetas[k]

    bias_df <- rbind.data.frame(bias_df, est_theta)
  }

  return(bias_df)
}

###
# Loop over sample sizes ------------------
###

samp.sizes

B <- 1000
allthebayas <- NULL

for(i in 1:length(samp.sizes)){
  ith_bayas <- bayas_fx(samplesize = samp.sizes[i], B = B)
  ith_bayas$sampsize <- samp.sizes[i]
  allthebayas <- rbind.data.frame(allthebayas, ith_bayas)
}

save(allthebayas, file = paste0("Ch3_samplesize/simdata/Bias", scenario, ".RData"))

allthebayas %>%
  mutate(sampsize = factor(sampsize),
         thresh = factor(thresh)) %>%
  group_by(sampsize, thresh) %>%
  summarise(bias_hat = (1/B)*sum(kth_theta-tru_theta),
            unbiased1 = tru_theta - bias_hat,
            unbiased2 = (2*tru_theta) - mean(kth_theta)) -> summ_bias

summ_bias %>%
  ggplot(., aes(x= thresh, y = bias_hat, color = sampsize)) +
  geom_point() +
  #scale_y_log10(name = "Log(bias_hat") +
  #scale_y_sqrt(name = "sqrt(bias_hat") +
  theme_bw()



