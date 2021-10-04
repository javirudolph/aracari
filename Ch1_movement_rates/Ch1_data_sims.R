# Script to run simulations using parameters from aracari movement rates
# LIBRARIES --------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(fitdistrplus)
library(Hmisc)

set.seed(98)

# LOAD the data ------------------------------------------------------------------------

data(ptpl)

null_moverate <- data.frame(Bird_ID = unique(ptpl$Bird_ID), movrate = mean(ptpl$mpm))

indiv_moverate <- ptpl %>%
  group_by(Bird_ID, fam_g) %>%
  summarise(movrate = mean(mpm))

fam_moverate <- ptpl %>%
  group_by(fam_g) %>%
  summarise(movrate = mean(mpm))

ids <- ptpl %>% distinct(., Bird_ID, fam_g)

# Visualize the movement rates on top of a density distribution, so that we can "assume" these movement rates are sampled from that distribution. The case of theory is with a lognormal distribution.

indiv_moverate %>%
  arrange(., movrate) %>%
  ggplot(., aes(x = movrate)) +
  geom_histogram()


# Average movement rate for this 'population'
mean(indiv_moverate$movrate)
mean(ptpl$mpm)
sum.indivs <- summary(indiv_moverate$movrate)


# Color palette ----------------------------------------
my.cols1 <- c("#23262f","#717492","#b3a82a","#c94f21","#980012","#0d907a","#b9bec3")

# Vis Movement Rates ---------------------------------------------
# Fit a lognormal to the distribution of movement rates
# Keep in mind these are all males.

logfit <- fitdist(indiv_moverate$movrate, distr = 'lnorm')
# normfit <- fitdist(indiv_moverate$movrate, distr = "norm")
# Mean using the fit
exp(logfit$estimate[1])

# Visualize distribution and data points
vlines <- summary(indiv_moverate$movrate)[c(1,3,6)]

indiv_moverate %>%
  ggplot(., aes(x = movrate, y = 0.002)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_function(fun = dlnorm, args = list(meanlog = logfit$estimate[1], sdlog = logfit$estimate[2]),
                color = "black", alpha = 0.8, size = 1) +
  labs(x = "Movement Rate (meters/minute)", y = "Density") +
  geom_vline(xintercept = vlines, color = my.cols1[4], size = 1, lty = 2) +
  scale_x_continuous(limits = c(0, 70), breaks = c(0, vlines[1], 20, vlines[2], 40, vlines[3], 60),
                     labels = c(0, "Min.", 20, "Median.", 40, "Max.", 60)) +
  #scale_x_discrete(labels = c(0, vlines[1], 20, vlines[2], 40, vlines[3], 60)) +
  theme_bw()

ggsave2(filename = "Ch1_movement_rates/Figures/Lnorm_fit2data.png", width = 6, height = 4, units = "in")

# Complete pooling ---------------------------
# We've fitted this lognormal distribution to our movement rates of 12 individuals (12 movement rates). The average movement rate for this sample of the population (the 12 individuals) is the average of the movement rate for each individual:
mean(indiv_moverate$movrate)

# 27.94211 meters/minute

# In a complete pooling approach, we take this average movement rate and use it as the parameter for our exponential distribution, from which we draw the movement distance at each time step.

## Compare exponential movement distances to full dataset ----------------------------------

movrate_cp <- mean(indiv_moverate$movrate)

# test <- data.frame(mpm = rexp(500, rate = 1/movrate_cp))
#
# test %>%
#   ggplot(., aes(x = mpm)) +
#   geom_histogram(aes(y = ..density..)) +
#   stat_function(fun = dexp, args = list(rate = 1/movrate_cp), color = "red")

ptpl %>%
  ggplot(., aes(x = mpm)) +
  geom_histogram(aes(y = ..density..)) +
  stat_function(fun = dexp, args = list(rate = 1/movrate_cp), color = my.cols1[4], size = 1) +
  stat_function(fun = dexp, args = list(rate = 1/exp(logfit$estimate[1])), color = my.cols1[3], size = 1) +
  labs(title = "Hiistogram of movement distances per minute from data set",
       x = "Movement distance per minute",
       y = "Density",
       caption = "Orange line is exponential density using average movement rate for population \n
       Yellow line shows exponential density using mean estimate from lognormal fit to movement rates") +
  theme_bw()








## Simulate movement rates -------------------------------------------------------------
n.individuals <- 20
m_1 <- sort(round(rlnorm(n.individuals, meanlog = logfit$estimate[1], sdlog = logfit$estimate[2]),3))

logfitdraw <- data.frame(indv = 1:20, movrate = m_1)

mycols <- viridis::magma(20)

indiv_moverate %>%
  ggplot(., aes(x = movrate, y = 0)) +
  geom_point(size = 5) +
  stat_function(fun = dlnorm, args = list(meanlog = logfit$estimate[1], sdlog = logfit$estimate[2]), color = "black", alpha = 0.8) +
  geom_point(data = logfitdraw, aes(x = movrate, color = factor(indv), y = 0.01), size = 5) +
  # scale_color_viridis_d(option="magma") +
  scale_color_viridis_d() +
  xlim(0, 60) +
  labs(title = "Mov rate sampling lnorm") +
  theme_bw() -> plot.logfitdraw
plot.logfitdraw

### Fist population--------------------------------------------------------------------

expcurves <- purrr::map(seq(1:n.individuals), function(y)
  stat_function(fun = dexp, args = list(rate = 1/logfitdraw$movrate[y]), color = mycols[y], alpha = 0.8))

logfitdraw %>%
  ggplot() +
  expcurves +
  xlim(-1, 125) +
  labs(title = "MD dist exp", x = 'movement distance') +
  theme_bw() -> expdraws
expdraws

cowplot::plot_grid(plot.logfitdraw, expdraws)

