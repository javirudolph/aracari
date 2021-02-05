load("paper/fevd.RData")

library(dplyr)
library(ggplot2)
library(cowplot)
library(extRemes)
library(Rmisc)

#*****************************************************************************************

# Find source for this:

# the influence of the outlier is much smaller when using L-moments estimation.

# 3. Plot the distribution
#   a. Plot density with estimated parameters
#   b. Plot density with the 95% confidence intervals
#   c. Plot bootstrap? estimates: sample 10000 from the data, 1000 times and get those parameters. Make those density plot


# Don't know if this is called bootstrapping


null_boot <- NULL
null_boot_probs <- NULL
for(i in 1:1000){
  samp <- sample(null_dispersal$dispersal, 500)
  fitD <- fevd(samp, threshold = null_thresh, type = "GP")
  scale <- fitD$results$par[1]
  shape <- fitD$results$par[2]
  out <- data.frame(boot = paste0("boot_", i), scale = scale, shape = shape)
  null_boot <- rbind.data.frame(null_boot, out)

  dist_range <- seq(null_thresh, 5000, by = 100)
  out2 <- data.frame(boot = paste0("boot_", i), distance = dist_range,
                     prob = pextRemes(null_fit, dist_range, lower.tail = FALSE))
  null_boot_probs <- rbind.data.frame(null_boot_probs, out2)
}

null_lines <- purrr::map(seq(1:1000), function(y)
  stat_function(fun = devd, args = list(type = "GP",
                                        scale = null_boot$scale[y], shape = null_boot$shape[y]), color = "grey"))

# ggplot(data = data.frame(x = c(null_thresh, 1000)), aes(x)) +
#   null_lines +
#   stat_function(aes(linetype = "CI"), fun = devd, args = list(scale = CI(null_boot$scale)[1], shape = CI(null_boot$shape)[1], type = "GP")) +
#   stat_function(aes(linetype = "Model"), fun = devd, args = list(scale = null_scale, shape = null_shape, type = "GP")) +
#   stat_function(aes(linetype = "CI"), fun = devd, args = list(scale = CI(null_boot$scale)[3], shape = CI(null_boot$shape)[3], type = "GP")) +
#   ylab("Density") + xlab("Distance (m)") +
#   theme_bw()  +
#   theme(legend.position = c(0.7, 0.7),
#         legend.title = element_blank()) -> null_boot_plot

ggplot(data = data.frame(x = c(null_thresh, 1000)), aes(x)) +
  null_lines +
  stat_function(aes(linetype = "CI"), fun = devd, args = list(scale = null_ci[1,1], shape = null_ci[2,1], type = "GP")) +
  stat_function(aes(linetype = "Model"), fun = devd, args = list(scale = null_ci[1,2], shape = null_ci[2,2], type = "GP")) +
  stat_function(aes(linetype = "CI"), fun = devd, args = list(scale = null_ci[1,3], shape = null_ci[2,3], type = "GP")) +
  ylab("Density") + xlab("Distance (m)") +
  theme_bw()  +
  theme(legend.position = c(0.7, 0.7),
        legend.title = element_blank()) -> null_boot_plot

indiv_boot <- NULL
indiv_boot_probs <- NULL
for(i in 1:1000){
  samp <- sample(indiv_dispersal$dispersal, 500)
  fitD <- fevd(samp, threshold = indiv_thresh, type = "GP")
  scale <- fitD$results$par[1]
  shape <- fitD$results$par[2]
  out <- data.frame(boot = paste0("boot_", i), scale = scale, shape = shape)
  indiv_boot <- rbind.data.frame(indiv_boot, out)

  dist_range <- seq(indiv_thresh, 5000, by = 100)
  out2 <- data.frame(boot = paste0("boot_", i), distance = dist_range,
                     prob = pextRemes(indiv_fit, dist_range, lower.tail = FALSE))
  indiv_boot_probs <- rbind.data.frame(indiv_boot_probs, out2)
}

indiv_lines <- purrr::map(seq(1:1000), function(y)
  stat_function(fun = devd, args = list(type = "GP",
                                        scale = indiv_boot$scale[y], shape = indiv_boot$shape[y]), color = "grey"))

ggplot(data = data.frame(x = c(indiv_thresh, 1000)), aes(x)) +
  indiv_lines +
  stat_function(aes(linetype = "CI"), fun = devd, args = list(scale = indiv_ci[1,1], shape = indiv_ci[2,1], type = "GP")) +
  stat_function(aes(linetype = "Model"), fun = devd, args = list(scale = indiv_ci[1,2], shape = indiv_ci[2,2], type = "GP")) +
  stat_function(aes(linetype = "CI"), fun = devd, args = list(scale = indiv_ci[1,3], shape = indiv_ci[2,3], type = "GP")) +
  ylab("Density") + xlab("Distance (m)") +
  theme_bw()  +
  theme(legend.position = c(0.7, 0.7),
        legend.title = element_blank()) -> indiv_boot_plot

fam_boot <- NULL
fam_boot_probs <- NULL
for(i in 1:1000){
  samp <- sample(fam_dispersal$dispersal, 500)
  fitD <- fevd(samp, threshold = fam_thresh, type = "GP")
  scale <- fitD$results$par[1]
  shape <- fitD$results$par[2]
  out <- data.frame(boot = paste0("boot_", i), scale = scale, shape = shape)
  fam_boot <- rbind.data.frame(fam_boot, out)

  dist_range <- seq(fam_thresh, 5000, by = 100)
  out2 <- data.frame(boot = paste0("boot_", i), distance = dist_range,
                     prob = pextRemes(fam_fit, dist_range, lower.tail = FALSE))
  fam_boot_probs <- rbind.data.frame(fam_boot_probs, out2)
}

fam_lines <- purrr::map(seq(1:1000), function(y)
  stat_function(fun = devd, args = list(type = "GP",
                                        scale = fam_boot$scale[y], shape = fam_boot$shape[y]), color = "grey"))

ggplot(data = data.frame(x = c(fam_thresh, 1000)), aes(x)) +
  fam_lines +
  stat_function(aes(linetype = "CI"), fun = devd, args = list(scale = fam_ci[1,1], shape = fam_ci[2,1], type = "GP")) +
  stat_function(aes(linetype = "Model"), fun = devd, args = list(scale = fam_ci[1,2], shape = fam_ci[2,2], type = "GP")) +
  stat_function(aes(linetype = "CI"), fun = devd, args = list(scale = fam_ci[1,3], shape = fam_ci[2,3], type = "GP")) +
  ylab("Density") + xlab("Distance (m)") +
  theme_bw()  +
  theme(legend.position = c(0.7, 0.7),
        legend.title = element_blank()) -> fam_boot_plot

plot_grid(null_boot_plot,
          indiv_boot_plot,
          fam_boot_plot,
          nrow = 3)

# The confidence intervals are for the bootstrapped data, not for the MLE estimate
# Figure out how to do that.



ggplot(data = data.frame(x = c(160, 1000)), aes(x)) +
  stat_function(aes(linetype = "Null"), fun = devd, args = list(scale = null_ci[1,2], shape = null_ci[2,2], type = "GP")) +
  stat_function(aes(linetype = "Individual"), fun = devd, args = list(scale = indiv_ci[1,2], shape = indiv_ci[2,2], type = "GP")) +
  stat_function(aes(linetype = "Family"), fun = devd, args = list(scale = fam_ci[1,2], shape = fam_ci[2,2], type = "GP")) +
  ylab("Density") + xlab("Distance (m)") +
  theme_bw()  +
  theme(legend.position = c(0.7, 0.7),
        legend.title = element_blank())


# 4. Calculate probability of getting those LDD events
#   a. make a figure with probability on the y axis, and distance in the x.

dist_range <- c(250, 500, 1000, 1250, 1500, 1750, 2000)
null_probs <- pextRemes(null_fit, dist_range, lower.tail = FALSE)
indiv_probs <- pextRemes(indiv_fit, dist_range, lower.tail = FALSE)
fam_probs <- pextRemes(fam_fit, dist_range, lower.tail = FALSE)

evd_probs <- as.data.frame(rbind(null_probs, indiv_probs, fam_probs))
names(evd_probs) <- dist_range
row.names(evd_probs) <- c("Null", "Individual", "Family")

save.image("paper/resampling_evd_densplots.RData")
