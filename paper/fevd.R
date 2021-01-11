
# Fit extreme value distributions to the dispersal kernel data

load("paper/dispersal_kernels.RData")

library(dplyr)
library(ggplot2)
library(cowplot)
library(extRemes)

#################################################################################################################
# Just testing

null_sample <- sample(null_dispersal$dispersal, 10000)
hist(null_sample)
#hist(log(null_sample))
#qqnorm(log(null_sample))
qqnorm(null_sample)
# Qq plots don't show straight lines, not consistent with a normal distribution of the data,

a <- threshrange.plot(null_sample, r = c(300, 700), nint = 200)
b <- mrlplot(null_sample, xlim = c(300, 700), nint = 100)

c <- decluster(null_sample, threshold = 500)
plot(c)

fit_D <- fevd(null_sample, threshold = 500, type = "GP")
#fit_D$results
summary(fit_D)
plot(fit_D)
t <- plot(fit_D, type = "qq")
pextRemes(fit_D, c(500, 1000, 1500, 2000, 3000, 3500, 4000), lower.tail = FALSE)

scale <- fit_D$results$par[1]
shape <- fit_D$results$par[2]
dens <- revd(10000, scale = scale, shape = shape, type = "GP")
hist(dens)
d <- density(dens)
plot(d)

ggplot(data = data.frame(x = c(0.1, 6)), aes(x)) +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = 1, shape = 0)) +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = 1, shape = 0.5))
  ylab("") +
  scale_y_continuous(breaks = NULL)


lines <- purrr::map(seq(0, 1, 0.1), function(y) stat_function(fun = devd, args = list(type = "GP", scale = 1, shape = y), color = "grey"))

ggplot(data = data.frame(x = c(0.1, 6)), aes(x)) +
  ylab("") +
  scale_y_continuous(breaks = NULL) -> p

p + lines +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = 1, shape = 0), linetype = "dashed", size = 2) +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = 1, shape = 0.5))

################################################################################################################################

#*****************************************************************************************
# NULL MODEL
orig_thres <- 500
# 1. Determine threshold

null_threshplot <- threshrange.plot(null_dispersal$dispersal, type = "GP", nint = 200)
null_mrl <- mrlplot(null_dispersal$dispersal, nint = 200)

null_thresh <- 500
# 2. Fit Generalized Pareto Distribution

null_fit <- fevd(null_dispersal$dispersal, threshold = null_thresh, type = "GP")
null_qq <- plot(null_fit, type = "qq")
null_summ <- summary(null_fit)
null_scale <- null_summ$par[1]
null_shape <- null_summ$par[2]
null_scale_se <- null_summ$se.theta[1]
null_shape_se <- null_summ$se.theta[2]

# 3. Plot the distribution
#   a. Plot density with estimated parameters
#   b. Plot density with the 95% confidence intervals
#   c. Plot bootstrap? estimates: sample 10000 from the data, 1000 times and get those parameters. Make those density plots

ggplot(data = data.frame(x = c(100, 2000)), aes(x)) +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = null_scale, shape = null_shape), size = 1) -> null_base_plot


# Don't know if this is called bootstrapping

null_boot <- NULL
null_boot_probs <- NULL
for(i in 1:1000){
  samp <- sample(null_dispersal$dispersal, 1000)
  fitD <- fevd(samp, threshold = null_thresh, type = "GP")
  scale <- fitD$results$par[1]
  shape <- fitD$results$par[2]
  out <- data.frame(boot = paste0("boot_", i), scale = scale, shape = shape)
  null_boot <- rbind.data.frame(null_boot, out)
  dist_range <- seq(null_thresh, 4000, by = 100)
  out2 <- data.frame(distance = dist_range,
                           prob = pextRemes(null_fit, dist_range, lower.tail = FALSE))
  null_boot_probs <- rbind.data.frame(null_boot_probs, out2)
}

null_lines <- purrr::map(seq(1:1000), function(y) stat_function(fun = devd, args = list(type = "GP", scale = null_boot$scale[y], shape = null_boot$shape[i]), color = "grey"))

# CI
hist(null_boot$scale)
hist(null_boot$shape)

# Calculate confidence interval according to t-distribution
confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}

null_ci_scale <- confidence_interval(null_boot$scale, 0.95)
null_ci_shape <- confidence_interval(null_boot$shape, 0.95)

null_base_plot + null_lines +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = null_scale, shape = null_shape), size = 1) +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = null_ci_scale$lower, shape = null_ci_shape$lower),
                linetype = "dashed", size = 1) +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = null_ci_scale$upper, shape = null_ci_shape$upper),
                linetype = "dashed", size = 1)

ggsave("null_dens.png")

# 4. Calculate probability of getting those LDD events
#   a. make a figure with probability on the y axis, and distance in the x.

dist_range <- seq(null_thresh, 4000, by = 100)
null_probs <- data.frame(distance = dist_range,
                         prob = pextRemes(null_fit, dist_range, lower.tail = FALSE))









