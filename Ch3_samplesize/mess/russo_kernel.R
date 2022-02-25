
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

# 5) Dispersal Kernel from Russo et.al. 2006
scenario <- "_russo"

mode1 <- data.frame(dist = "norm", mean = 5.5, sd = 2.15) #non-dispersed
mode2 <- data.frame(dist = "weibull", shape = 3.26, scale = 221.97) # p=shape, s=scale to sleeping sites
mode3 <- data.frame(dist = "lnorm", mean = 4.50, sd = 1.71) # in transit <500m
mode4 <- data.frame(dist = "norm", mean = 901.5, sd = 121.23) # in transit >500

prms.df <- data.frame(dist = c("norm", "weibull", "lnorm", "norm"),
                      param1 = c(5.5, 3.26, 4.5, 901.5),
                      param2 = c(2.15, 221.97, 1.71, 121.23),
                      mean = c(5.5))

# assume 92% of seeds dispersed happened by spider monkeys. Of those, half were to sleeping sites and half in transit.
pis <- c(0.2, 0.4, 0.3, 0.1)
pis <- c(0.25, 0.25, 0.25, 0.25)


# MIXTURE MODEL --------------------------------------------------
## Mixture components --------------------------------------

dens_cols <- c("#264653", "#2a9d8f", "#f4a261", "#e76f51")


### Plot the mixture components -------------------------------


# But choosing to use the weighted curves instead
ggplot() +
  stat_function(fun = dnorm,
                args = list(mean = 5.5, sd = 2.15),
                color = dens_cols[1], size=pis[1]*5, alpha = 0.8)+
  stat_function(fun = dweibull,
                args = list(shape = 3.26, scale = 221.97),
                color = dens_cols[2], size=pis[2]*5, alpha = 0.8)+
  stat_function(fun = dlnorm,
                args = list(meanlog = 4.5, sdlog = 1.71),
                color = dens_cols[3], size=pis[3]*5, alpha = 0.8)+
  stat_function(fun = dnorm,
                args = list(mean = 901.5, sd = 121.23),
                color = dens_cols[4], size=pis[4]*5, alpha = 0.8)+
  theme_minimal() +
  labs(y = "Density") +
  #scale_y_sqrt() +
  lims(x = c(0, 1500)) -> densities_plot
densities_plot

top_row <- plot_grid(densities_plot)
top_row

# SAMPLING --------------------------------------------------------------

## Draw from mixture ----------------------------------------------------
# Using purrr to get samples instead of a for loop. More efficient.

# We are getting 50k samples from the mixture
samp.n <- 50000


# These category samples are based on the weights
cat.samp.sizes <- samp.n*pis

xmode1 <- data.frame(gID = 1, x.samps = rnorm(cat.samp.sizes[1], mean = prms.df$param1[1], sd = prms.df$param2[1]))
xmode2 <- data.frame(gID = 2, x.samps = rweibull(cat.samp.sizes[2], shape = prms.df$param1[2], scale = prms.df$param2[2]))
xmode3 <- data.frame(gID = 3, x.samps = rlnorm(cat.samp.sizes[3], meanlog = prms.df$param1[3], sdlog = prms.df$param2[3]))
xmode4 <- data.frame(gID = 4, x.samps = rnorm(cat.samp.sizes[4], mean = prms.df$param1[4], sd = prms.df$param2[4]))

simplsamps <- rbind.data.frame(xmode1, xmode2, xmode3, xmode4)

save(simplsamps, file = paste0("Ch3_samplesize/simdata/mixturesamples", scenario, ".RData"))

# Extract the tail (the highest values) to visualize later in the histogram since they are too few to show in the bins
data.tail <- data.frame(values = simplsamps$x.samps, y = 100) %>% arrange(desc(values)) %>% filter(values >=100)

## Visualization and save -----------------------------------------------------

# Histogram plus the points in the tail.
simplsamps %>%
  ggplot(., aes(x = x.samps)) +
  geom_histogram(bins = 100) +
  geom_point(data = data.tail[1:50,], aes(x = values, y = y), color = "black", alpha = 0.5) +
  labs(y = "Frequency", x = "Distance") +
  #lims(x = c(0, 150)) +
  theme_minimal() -> sampleshist
sampleshist

# Density curve for the mixture based on the 50k samples.

simplsamps %>%
  ggplot(., aes(x = x.samps)) +
  geom_density() +
  labs(y = "Density", x = "Distance") +
  lims(x = c(0, 1500)) +
  # scale_y_sqrt(name = "Sqrt(density)", limits = c(0, 0.005)) +
  theme_minimal() -> samplesdens
samplesdens

# Visualize both in one frame
plot_grid(sampleshist, samplesdens)


# Generate Figure 1 that describes the mixture components, and the samples

bottom_row <- plot_grid(sampleshist, samplesdens)

plot_grid(top_row, bottom_row,nrow = 2)

ggsave(paste0("Ch3_samplesize/Figures/Figure1", scenario,".png"))


###
## Probability of tails
###

# Explore the quantiles
quantile(simplsamps$x.samps, 0.9)
summary(simplsamps$x.samps)

# This threshold is between the 75-100% quantiles.
# thresh.vals <- seq(50, 500, length.out = 5)

# new thresholds
thresh.vals <- c(50, 75, 100, 150, 250, 500, 1000)


tibble(thresh.vals) %>%
  mutate(n.overthresh = map_dbl(1:length(thresh.vals),
                                function(y) length(which(simplsamps$x.samps >= thresh.vals[y]))),
         theta = n.overthresh/samp.n) -> tru.cdfs
tru.cdfs


###
## FITTING EVD ----------------------------------
##

library(extRemes)

samp.sizes <- c(80, 200, 500, 800, 1000, 1600)


fevd.mles <- data.frame(expand.grid(threshold = thresh.vals, samp.size = samp.sizes))
for(i in 1:nrow(fevd.mles)){
  ith.n        <- fevd.mles$samp.size[i]
  ith.samples  <- data.frame(x = sample(simplsamps$x.samps, ith.n))
  # So, setting the threshold to 0 so we can compare to Lomax.
  ith.fit      <- fevd(ith.samples$x, threshold = 0, type = "GP")
  mles          <- summary(ith.fit)$par
  nll.hat       <- summary(ith.fit)$nllh
  BIC.mod       <- summary(ith.fit)$BIC
  fevd.mles$scale[i] <- mles[1]
  fevd.mles$shape[i] <- mles[2]
  fevd.mles$nllh[i]  <- nll.hat
  ith.thresh   <- fevd.mles$threshold[i]
  ith.tail <- length(which(ith.samples$x >= ith.thresh))
  fevd.mles$samps.tail[i] <- ith.tail/ith.n
  fevd.mles$gp.tail[i] <- pextRemes(ith.fit, ith.thresh, lower.tail = FALSE)

  ith.lomax <- lomax.glm(formula = ~1, ith.samples, ith.samples$x)
  ith.alpha <- ith.lomax$alphas.hat[1]
  ith.k     <- ith.lomax$k.hat
  fevd.mles$alpha[i] <- ith.alpha
  fevd.mles$k[i] <- ith.k
  fevd.mles$lomax.tail[i] <- lomax.st(ith.thresh, alpha = ith.alpha, k = ith.k)

  fevd.mles$calc.k[i] <- 1/mles[2]
  fevd.mles$calc.alpha[i] <- mles[1]/mles[2]
  fevd.mles$calc.tail[i] <- lomax.st(ith.thresh, alpha = mles[1]/mles[2], k = 1/mles[2])
}

### Plot parameter space ---------------------------

fevd.mles %>%
  mutate(samp.size = factor(samp.size),
         threshold = factor(threshold)) %>%
  ggplot(., aes(x = shape, y = scale, color = threshold)) +
  facet_wrap(~samp.size) +
  geom_point(size = 2) +
  labs(title = "GP fits") +
  theme_bw() -> gp.param.space
gp.param.space

fevd.mles %>%
  mutate(samp.size = factor(samp.size),
         threshold = factor(threshold)) %>%
  ggplot(., aes(x = alpha, y = k, color = threshold)) +
  facet_wrap(~samp.size) +
  geom_point(size = 2) +
  labs(title = "Lomax fits") +
  theme_bw() -> lomax.param.space
lomax.param.space

param.space.legend <- get_legend(lomax.param.space)

plot_grid(gp.param.space + theme(legend.position = "none"), lomax.param.space + theme(legend.position = "none"), param.space.legend, rel_widths = c(1,1,0.2), ncol = 3)

ggsave(paste0("Ch3_samplesize/Figures/Param_space", scenario,".png"))


# Plot tail estimates -----------------------------

fevd.mles$w.cdfs <- tru.cdfs$w.cdfs
fevd.mles$theta <- tru.cdfs$theta

fevd.mles %>%
  mutate(threshold = factor(threshold),
         samp.size = factor(samp.size),
         tail.ratio = gp.tail/theta) %>%
  ggplot(., aes(x = samp.size, y = tail.ratio, color = threshold)) +
  geom_point() +
  labs(y = "GP tail.ratio") +
  theme_bw() +
  theme(legend.position = "none") -> a
a


fevd.mles %>%
  mutate(threshold = factor(threshold),
         samp.size = factor(samp.size),
         tail.ratio = lomax.tail/theta) %>%
  ggplot(., aes(x = samp.size, y = tail.ratio, color = threshold)) +
  geom_point() +
  labs(y = "Lomax tail.ratio") +
  theme_bw() -> b
b

plot_grid(a,b, rel_widths = c(0.7, 1))
ggsave(paste0("Ch3_samplesize/Figures/tail_ratio", scenario, ".png"))



## MC Samples -----------------------------------------------

#CHANGE nreps ================
nreps <- 1000

gp.mles.reps <- data.frame(NULL)

for(j in 1:nreps){

  ith.mles.df <- data.frame(expand.grid(threshold = thresh.vals, samp.size = samp.sizes))
  for(i in 1:nrow(ith.mles.df)){
    ith.n        <- ith.mles.df$samp.size[i]
    ith.samples  <- data.frame(x = sample(simplsamps$x.samps, ith.n))
    # So, setting the threshold to 0 so we can compare to Lomax.
    ith.fit      <- fevd(ith.samples$x, threshold = 0, type = "GP")
    mles          <- summary(ith.fit)$par
    nll.hat       <- summary(ith.fit)$nllh
    BIC.mod       <- summary(ith.fit)$BIC
    ith.mles.df$scale[i] <- mles[1]
    ith.mles.df$shape[i] <- mles[2]
    ith.mles.df$nllh[i]  <- nll.hat
    ith.thresh   <- ith.mles.df$threshold[i]
    ith.tail <- length(which(ith.samples$x >= ith.thresh))
    ith.mles.df$samps.tail[i] <- ith.tail/ith.n
    ith.mles.df$gp.tail[i] <- pextRemes(ith.fit, ith.thresh, lower.tail = FALSE)

    ith.lomax <- lomax.glm(formula = ~1, ith.samples, ith.samples$x)
    ith.alpha <- ith.lomax$alphas.hat[1]
    ith.k     <- ith.lomax$k.hat
    ith.mles.df$alpha[i] <- ith.alpha
    ith.mles.df$k[i] <- ith.k
    ith.mles.df$lomax.tail[i] <- lomax.st(ith.thresh, alpha = ith.alpha, k = ith.k)

    fevd.mles$calc.k[i] <- 1/mles[2]
    fevd.mles$calc.alpha[i] <- mles[1]/mles[2]
    fevd.mles$calc.tail[i] <- lomax.st(ith.thresh, alpha = mles[1]/mles[2], k = 1/mles[2])


    ith.mles.df$rep[i] <- paste0("rep", j)

  }

  gp.mles.reps <- rbind.data.frame(gp.mles.reps, ith.mles.df)

}

gp.mles.reps$w.cdfs <- tru.cdfs$w.cdfs
gp.mles.reps$theta <- tru.cdfs$theta


save(gp.mles.reps, file = paste0("Ch3_samplesize/simdata/GPmles", scenario, ".RData"))


# Ratio Boxplots -------------------

gp.mles.reps %>%
  mutate(lomax.tail.ratio = lomax.tail/theta,
         gp.tail.ratio = gp.tail/theta,
         threshold = factor(threshold),
         samp.size = factor(samp.size)) %>%
  ggplot(., aes(x = samp.size, y = gp.tail.ratio, color = threshold)) +
  geom_boxplot() +
  # geom_point()
  labs(y = "GP - i/theta") +
  #lims(y = c(-0.1, 5)) +
  theme_bw()->box1
box1
#
#
# gp.mles.reps %>%
#   mutate(gp.tail.ratio = gp.tail/samps.tail,
#          threshold = factor(threshold),
#          samp.size = factor(samp.size)) %>%
#   ggplot(., aes(x = samp.size, y = gp.tail.ratio, color = threshold)) +
#   geom_boxplot() +
#   labs(y = "GP - i/samps.tail") +
#   #lims(y = c(-0.1, 5)) +
#   theme_bw()->b
# b
#
# c <- get_legend(b)
#
# plot_grid(a + theme(legend.position = "none"),
#           b + theme(legend.position = "none"),
#           c, rel_widths = c(1,1,0.2), ncol = 3)


ggsave(paste0("Ch3_samplesize/Figures/GPtail_boxplt", scenario, ".png"), width = 10, height = 4)

# Mean and Standard Error plots -------------

head(gp.mles.reps)

gp.mles.reps %>%
  mutate(gp.tail.ratio = gp.tail/theta,
         threshold = factor(threshold),
         samp.size = factor(samp.size)) %>%
  #filter(gp.tail != 0) %>%
  mutate(sqrd_diff = (gp.tail-theta)^2) %>%
  group_by(samp.size, threshold) %>%
  summarise(mean.ratio = mean(gp.tail.ratio),
            ste.gp.ratio = sd(gp.tail.ratio)/sqrt(length(gp.tail.ratio)),
            mse = (1/length(gp.tail.ratio))*sum(sqrd_diff)) -> tail.ratio.means
tail.ratio.means

tail.ratio.means %>%
  ggplot(., aes(color = threshold,
                y = mse,
                x = samp.size)) +
  geom_jitter(size = 2, width = 0.2, alpha = 0.8) +
  scale_y_log10(name = "Log(mse)") +
  labs(subtitle = "Is MSE smaller or is it driven by zeroes?") +
  theme_bw() -> a
a


tail.ratio.means %>%
  mutate(lo = mean.ratio - ste.gp.ratio,
         hi = mean.ratio + ste.gp.ratio) %>%
  ggplot(., aes(x = threshold,
                y = mean.ratio)) +
  facet_wrap(~samp.size) +
  geom_point() +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0) +
  # scale_y_log10(name = "Log(mean.ratio)") +
  theme_bw() -> b
b

plot_grid(a,b, nrow=2, rel_heights = c(1,1.5))


ggsave(paste0("Ch3_samplesize/Figures/GPtail_mean", scenario, ".png"))



# UNBIASED ESTIMATOR ------------------------------------------------


# Bayas FXs ----------------------------

bias_fx <- function(B, theta_i, tru.theta){
  (1/B)*sum(theta_i - tru.theta)
}

bayas_fx <- function(samplesize = 100, B=10, thresh.vals = thresh.vals, data.vec = simplsamps$x.samps){
  ith.n <- samplesize # the sample size
  ith.df <- data.frame(x = sample(data.vec, ith.n))

  ith.tail <- map_dbl(1:length(thresh.vals), function(y) length(which(ith.df$x >= thresh.vals[y]))/ith.n)

  # Estimate the shape and scale parameters
  ith.evd <- fevd(x = ith.df$x, type = "GP", threshold = 0)
  ith.scale <- ith.evd$results$par[1]
  ith.shape <- ith.evd$results$par[2]


  # Using the same threshold values as above
  # These are the true thetas for the different thresholds
  ith.theta.hats <- pevd(thresh.vals, loc = 0, scale = ith.scale, shape= ith.shape, lower.tail = FALSE)


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

  for(k in 1:length(thresh.vals)){
    kth_thresh <- thresh.vals[k]

    k.frame <- jth_df %>%
      mutate(jth_theta = map_dbl(1:B, function(y) pevd(kth_thresh, loc = 0,
                                                       scale = jscale[y], shape = jshape[y], lower.tail = FALSE)))
    k.frame$thresh <- kth_thresh
    k.frame$ith.theta.hats <- ith.theta.hats[k]
    k.frame$ith.tail <- ith.tail[k]
    k.frame$ith.shape <- ith.shape
    k.frame$ith.scale <- ith.scale

    bias_df <- rbind.data.frame(bias_df, k.frame)
  }

  return(bias_df)
}

###
# MBAYAS Loop ------------------
###

# samp.sizes
#
# B <- 10
# allthebayas <- NULL
#
# for(i in 1:length(samp.sizes)){
#   ith_bayas <- bayas_fx(samplesize = samp.sizes[i], B = B)
#   ith_bayas$sampsize <- samp.sizes[i]
#   allthebayas <- rbind.data.frame(allthebayas, ith_bayas)
# }

mbaya <- NULL
nruns <- 100
B <- 1000

for(m in 1:nruns){

  allthebayas <- NULL

  for(i in 1:length(samp.sizes)){
    ith_bayas <- bayas_fx(samplesize = samp.sizes[i],
                          B = B, thresh.vals = thresh.vals, data.vec = simplsamps$x.samps)
    ith_bayas$sampsize <- samp.sizes[i]
    allthebayas <- rbind.data.frame(allthebayas, ith_bayas)
  }

  allthebayas$run <- paste0("run_", m)
  mbaya <- rbind.data.frame(mbaya, allthebayas)

}


# Ok, so the TRUE TRUE tail is the one coming from the original simulation using the mixtures
tru.cdfs

tru.cdfs %>%
  dplyr::select(thresh.vals, theta) %>%
  rename(thresh = thresh.vals) %>%
  right_join(., mbaya, by = "thresh") -> mbaya

save(mbaya, file = paste0("Ch3_samplesize/simdata/Bias_df", scenario, ".RData"))

# Need to check these, multiple rows are the same, the fx summarise isn't working well.
mbaya %>%
  mutate(sampsize = factor(sampsize),
         thresh = factor(thresh),
         theta.diff = jth_theta-ith.theta.hats) %>%
  group_by(run, sampsize, thresh) %>%
  summarise(theta = mean(theta),
            theta_hat = mean(ith.theta.hats),
            bias_hat = (1/B)*sum(theta.diff),
            theta_bar = (2*theta_hat) - mean(jth_theta)) -> summ_bias

summ_bias %>%
  ggplot(., aes(x= thresh, y = bias_hat, color = sampsize)) +
  geom_point() +
  #scale_y_log10(name = "Log(bias_hat") +
  #scale_y_sqrt(name = "sqrt(bias_hat") +
  theme_bw() +
  theme(legend.position = "none") -> baya1

summ_bias %>%
  ggplot(., aes(x= thresh, y = bias_hat, color = sampsize)) +
  geom_point() +
  scale_y_log10(name = "Log") +
  theme_bw() +
  theme(legend.position = "none") -> baya2

summ_bias %>%
  ggplot(., aes(x= thresh, y = bias_hat, color = sampsize)) +
  geom_point() +
  scale_y_sqrt(name = "sqrt") +
  theme_bw() -> baya3

top_baya <- plot_grid(baya1, baya2, ncol = 2)
plot_grid(top_baya, baya3, nrow = 2)

ggsave(paste0("Ch3_samplesize/Figures/Bias_hat", scenario, ".png"))

# True theta vs theta bar

summ_bias %>%
  ggplot(., aes(x = theta, y = theta_bar, color = thresh)) +
  geom_point() +
  theme_bw() -> a

summ_bias %>%
  ggplot(., aes(x = theta, y = theta_bar, color = sampsize)) +
  geom_point() +
  theme_bw() -> b

summ_bias %>%
  ggplot(., aes(x = sampsize, y = theta_bar/theta, color = thresh)) +
  geom_boxplot() +
  theme_bw() -> c

toprow <- plot_grid(a, b, nrow = 1)

plot_grid(toprow, c, nrow=2)

ggsave(paste0("Ch3_samplesize/Figures/theta_bar", scenario, ".png"))


