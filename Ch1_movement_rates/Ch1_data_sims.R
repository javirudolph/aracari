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

# Complete pooling ---------------------------------------------
# Fit a lognormal to the distribution of movement rates
# Keep in mind these are all males.

logfit <- fitdist(indiv_moverate$movrate, distr = 'lnorm')
# normfit <- fitdist(indiv_moverate$movrate, distr = "norm")
# Mean using the fit
exp(logfit$estimate[1])

# Visualize distribution and data points

indiv_moverate %>%
  ggplot(., aes(x = movrate, y = 0, color = fam_g)) +
  geom_point() +
  stat_function(fun = dlnorm, args = list(meanlog = logfit$estimate[1], sdlog = logfit$estimate[2]), color = "black", alpha = 0.8) +
  xlim(0, 50)

# Partial pooling ---------------------------------------------
# We don't have enough individuals to fit individual lognormal distributions to each family group


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
  scale_color_viridis_d(option="magma") +
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

