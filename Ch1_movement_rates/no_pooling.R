# Script to run simulations using parameters from aracari movement rates
# LIBRARIES --------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(fitdistrplus)
library(Hmisc)

set.seed(98)

# Source functions -----------------------------------------------------------------
source("Ch1_movement_rates/Ch1_functions.R")

# LOAD the data ------------------------------------------------------------------------

data(ptpl)

null_moverate <- data.frame(Bird_ID = unique(ptpl$Bird_ID), movrate = mean(ptpl$mpm))

indiv_moverate <- ptpl %>%
  group_by(Bird_ID, fam_g) %>%
  summarise(movrate = mean(mpm))

fam_moverate <- ptpl %>%
  group_by(fam_g) %>%
  summarise(movrate = mean(mpm),
            sd = sd(mpm))

ids <- ptpl %>% distinct(., Bird_ID, fam_g)

# Average movement rate for this 'population'
mean(indiv_moverate$movrate)
mean(ptpl$mpm)

# Color palette ----------------------------------------
my.cols1 <- c("#23262f","#717492","#b3a82a","#c94f21","#980012","#0d907a","#b9bec3")
# No pooling --------------------------------------------
# This one looks at individual variation.
# We will sample 20 individuals, 3 times from the lognormal distribution that describes movement rates for the population
logfit <- fitdist(indiv_moverate$movrate, distr = 'lnorm')

n.individuals <- 20
m_1 <- sort(round(rlnorm(n.individuals, meanlog = logfit$estimate[1], sdlog = logfit$estimate[2]),3))
m_2 <- sort(round(rlnorm(n.individuals, meanlog = logfit$estimate[1], sdlog = logfit$estimate[2]),3))
m_3 <- sort(round(rlnorm(n.individuals, meanlog = logfit$estimate[1], sdlog = logfit$estimate[2]),3))

logfitdraw <- data.frame(indv = 1:20, m_1 = m_1, m_2 = m_2, m_3 = m_3)

indiv_moverate %>%
  ggplot(., aes(x = movrate, y = 0)) +
  geom_point(size = 5) +
  stat_function(fun = dlnorm, args = list(meanlog = logfit$estimate[1], sdlog = logfit$estimate[2]), color = "black", alpha = 0.8) +
  geom_point(data = logfitdraw, aes(x = m_1, color = factor(indv), y = 0.005), size = 3) +
  geom_point(data = logfitdraw, aes(x = m_2, color = factor(indv), y = 0.006), size = 3) +
  geom_point(data = logfitdraw, aes(x = m_3, color = factor(indv), y = 0.007), size = 3) +
  # scale_color_viridis_d(option="magma") +
  scale_color_viridis_d() +
  xlim(0, 60) +
  labs(title = "Sampling individuals from lognormal",
       x = "Movement rate",
       y = "Density") +
  theme_bw() +
  theme(legend.position = "none") -> plot.logfitdraw

plot.logfitdraw

mycols <- viridisLite::viridis(20)

expcurves_m1 <- purrr::map(seq(1:n.individuals), function(y)
  stat_function(fun = dexp, args = list(rate = 1/logfitdraw$m_1[y]), color = mycols[y], alpha = 0.8))

expcurves_m2 <- purrr::map(seq(1:n.individuals), function(y)
  stat_function(fun = dexp, args = list(rate = 1/logfitdraw$m_2[y]), color = mycols[y], alpha = 0.8))

expcurves_m3 <- purrr::map(seq(1:n.individuals), function(y)
  stat_function(fun = dexp, args = list(rate = 1/logfitdraw$m_3[y]), color = mycols[y], alpha = 0.8))

logfitdraw %>%
  ggplot() +
  expcurves_m1 +
  expcurves_m2 +
  expcurves_m3 +
  xlim(-1, 125) +
  labs(title = "Movement distance Exponential Draws", x = 'Movement Distance', y = "Density") +
  theme_bw() -> expdraws
expdraws

cowplot::plot_grid(plot.logfitdraw, expdraws)


ggsave2(filename = "Ch1_movement_rates/Figures/Sampling_indivs_lnorm.png")

## Generate data -------------------------------------------------
m.data <- data.frame(m_1, m_2, m_3)
kruns <- 1000
nseeds <- 5

np.df <- NULL
np.summ.df <- NULL

for(m in 1:3){
  m.0 <- m.data[m]
  for(j in 1:n.individuals){
    for(k in 1:kruns){
      a <- sim_seeds(m.prms = m.0[j,], nseeds = nseeds) %>%
        mutate(indiv = as.factor(j),
               run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
               popu = as.factor(m))

      b <- summ_seeds(a) %>%
        mutate(indiv = as.factor(j),
               run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
               popu = as.factor(m))

      np.df <- rbind.data.frame(np.df, a)
      np.summ.df <- rbind.data.frame(np.summ.df, b)
    }
  }
}




save(np.df, np.summ.df, file = "Ch1_movement_rates/sims_backup/datagen_np.RData")


### Visualize -----
# np.summ.df %>%
#   group_by(., popu) %>%
#   sample_n(., 100) %>%
#   ggplot(., aes(y = dsprsn, x = factor(popu), color = factor(popu))) +
#   geom_boxplot() +
#   geom_point(color = "grey", alpha = 0.5) +
#   labs(title = "Dispersion") +
#   theme_bw() +
#   theme(legend.position = "none") -> np_p1
# # p1
#
# np.summ.df %>%
#   group_by(., popu) %>%
#   sample_n(., 100) %>%
#   ggplot(., aes(y = av.disp, x = factor(popu), color = factor(popu))) +
#   geom_boxplot() +
#   geom_point(color = "grey", alpha = 0.5) +
#   labs(title = "Average dispersal per run") +
#   theme_bw() +
#   theme(legend.position = "none") -> np_p2
# # p2
#
# plot_grid(np_p1, np_p2)

### Kernel -----
n.boots <- 100
samp.size <- 50
weib.boot.np <- NULL

for(j in 1:n.boots){
  s.df <- np.df %>% drop_na(s.id) %>%
    group_by(popu) %>%
    sample_n(., samp.size) %>%
    mutate(disp = round(disp, digits = 2))

  # s.df %>%
  #   ggplot(., aes(x = disp, fill = popu)) +
  #   geom_histogram()

  # s.df <- df %>% drop_na(s.id)

  weib.fits <- NULL

  for(i in 1:3){
    dat <- s.df %>%
      filter(popu == i)

    g <- fitdist(dat$disp, distr = "weibull", method = 'mle', lower = c(0,0))
    #g <- fitdistr(dat$disp, densfun = "weibull", lower = c(0,0))
    prms.weib <- data.frame(est.shape = as.numeric(g$estimate[1]),
                            est.scale = as.numeric(g$estimate[2]),
                            loglik = g$loglik,
                            popu = i)
    weib.fits <- rbind.data.frame(weib.fits, prms.weib)
    # plot(g)
  }

  weib.boot.np <- rbind.data.frame(weib.boot.np, weib.fits %>% mutate(boot = j))
}


save(weib.boot.np, file = "Ch1_movement_rates/sims_backup/weib_np.RData")

weib.boot.np %>%
  group_by(popu) %>%
  #sample_n(., 100) %>%
  ggplot(., aes(x = factor(popu), y = est.shape, color = factor(popu))) +
  geom_violin() +
  # scale_color_manual(values = c("black", mycols)) +
  # geom_boxplot(width = 0.01) +
  # geom_point(color = "grey", alpha = 0.5) +
  # stat_summary(fun.data=mean_sdl, mult=1,
  #              geom="pointrange", color="black") +
  geom_jitter(position = position_jitter(0.1)) +
  labs(title = "Shape") +
  theme_bw() -> np_weib_p1
# p1

weib.boot.np %>%
  group_by(popu) %>%
  #sample_n(., 100) %>%
  ggplot(., aes(x = factor(popu), y = est.scale, color = factor(popu))) +
  geom_violin() +
  # scale_color_manual(values = c("black", mycols)) +
  #geom_boxplot(width = 0.1) +
  # geom_point() +
  # stat_summary(fun.data=mean_sdl, mult=1,
  #              geom="pointrange", color="black") +
  geom_jitter(position = position_jitter(0.1)) +
  labs(title = "Scale") +
  theme_bw() -> np_weib_p2

plot_grid(np_weib_p1, np_weib_p2)
