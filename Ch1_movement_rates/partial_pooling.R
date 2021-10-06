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

fam_moverate

ptpl %>%
  distinct(fam_g, Bird_ID)

ptpl %>%
  dplyr::select(fam_g, Bird_ID) %>%
  distinct(Bird_ID, .keep_all = TRUE) %>%
  group_by(fam_g) %>%
  summarise(n = n()) %>%
  right_join(., fam_moverate) %>%
  mutate(logpar = log(movrate),
         fitted_sd = logfit$estimate[2]) -> fam_moverate
fam_moverate

# Using the average movement rate per family group as the meanlog for a lognormal distribution, and the sdlog from the fitted distribution to the whole dataset since we don't have enough individuals per fmaily group to fit a new lognormal.


# My question now is whether or not I would be considering individual variation. I'm saying no, because that's like douible variaion of no pooling? Like, we have a lognormal for each family group, and then get individuals from each family group. That is one thing.
# The other option is just using the average movement rates per family group to go into the exponential distribution for movement distances. This is what I would consider partial pooling.
# I guess the other is also a type of partial pooling, as you are only allowing the individual variation to occur within the bound of the family group expectations.

## Visualize mov rates ---------------
pp_1 <- sort(round(fam_moverate$movrate,3))

### Compare lognorm fits CP and PP -----------

logfit_fam <- fitdist(pp_1, distr = 'lnorm')
# normfit <- fitdist(indiv_moverate$movrate, distr = "norm")

# Mean using the fit
exp(logfit_fam$estimate[1])

# Visualize distribution and data points
vlines_fam <- summary(fam_moverate$movrate)[c(1,3,6)]

indiv_moverate %>%
  ggplot(., aes(x = movrate, y = 0.002)) +
  # stat_function(fun = dlnorm, args = list(meanlog = logfit$estimate[1], sdlog = logfit$estimate[2]),
  #               color = "black", alpha = 0.8, size = 1) +
  # geom_point(size = 3, alpha = 0.8) +
  # geom_vline(xintercept = vlines, color = "black", size = 1, lty = 2) +
  # scale_x_continuous(limits = c(0, 70), breaks = c(0, vlines[1], 20, vlines[2], 40, vlines[3], 60),
  #                    labels = c(0, "Min.", 20, "Median.", 40, "Max.", 60)) +
  stat_function(fun = dlnorm, args = list(meanlog = logfit_fam$estimate[1], sdlog = logfit_fam$estimate[2]),
                color = my.cols1[4], alpha = 0.8, size = 1) +
  geom_point(data = fam_moverate, aes(x = movrate, y = 0.004),size = 3, alpha = 0.8, color = my.cols1[4]) +
  geom_vline(xintercept = vlines_fam, color = my.cols1[4], size = 1, lty = 2) +
  scale_x_continuous(limits = c(0, 70), breaks = c(0, vlines_fam[1], 20, vlines_fam[2], vlines_fam[3], 60),
                     labels = c(0, "Min.", 20, "Median.", "Max.", 60)) +
  labs(x = "Movement Rate (meters/minute)", y = "Density",
       title = "Lognormal fit to family group movement rates") +
  theme_bw()

ggsave2(filename = "Ch1_movement_rates/Figures/logfit_pp.png")

## PP w/fam moverates ---------------------------------------
# Meaning, we use the average movement rate for each family group as the parameter for an exponential from which we draw movement distances for the simulation.

fam_moverate %>%
  arrange(., movrate) -> fam_moverate

fam_moverate %>%
  ggplot(., aes(x = movrate, y = 0, color = movrate)) +
  geom_point(size = 3) +
  stat_function(fun = dlnorm, args = list(meanlog = logfit_fam$estimate[1], sdlog = logfit_fam$estimate[2]),
                color = "black", alpha = 0.8, size = 1) +
  xlim(0, 60) +
  labs(title = "Lognormal fitted to family group",
       x = "Movement rate",
       y = "Density") +
  scale_color_viridis_c() +
  theme_bw() +
  theme(legend.position = "none") -> pp_logfit

mycols2 <- viridisLite::viridis(7)
expcurves_pp1 <- purrr::map(seq(1:7), function(y)
  stat_function(fun = dexp, args = list(rate = 1/fam_moverate$movrate[y]), color = mycols2[y], alpha = 0.8))

ggplot() +
  expcurves_pp1 +
  xlim(-1, 125) +
  labs(title = "Movement distance exponential draws", x = 'Movement Distance', y = "Density") +
  theme_bw() -> pp_expdraws


cowplot::plot_grid(pp_logfit, pp_expdraws)
ggsave2(filename = "Ch1_movement_rates/Figures/pp_logfit_plus_exp.png")

### PP Generate data ------------------------------------------------------

kruns <- 1500
nseeds <- 5

pp.df <- NULL
pp.summ.df <- NULL

for(j in 1:7){
  for(k in 1:kruns){
    a <- sim_seeds(m.prms = fam_moverate$movrate[j], nseeds = nseeds) %>%
      mutate(fam_g = fam_moverate$fam_g[j],
             run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
             model = paste0("pp"))

    b <- summ_seeds(a) %>%
      mutate(fam_g = fam_moverate$fam_g[j],
             run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
             model = paste0("pp"))

    pp.df <- rbind.data.frame(pp.df, a)
    pp.summ.df <- rbind.data.frame(pp.summ.df, b)
  }
}


save(pp.df, pp.summ.df, file = "Ch1_movement_rates/sims_backup/datagen_pp.RData")


# sample with replacement groups of 100 seeds and fit weibull , here are the parameters.
n.boots <- 1000
samp.size <- 100
weib.boot.pp <- NULL

for(j in 1:n.boots){
  s.df <- pp.df %>% drop_na(s.id) %>%
    sample_n(., samp.size) %>%
    mutate(disp = round(disp, digits = 2))

  g <- fitdist(s.df$disp, distr = "weibull", method = 'mle', lower = c(0,0))


  prms.weib <- data.frame(est.shape = as.numeric(g$estimate[1]),
                          est.scale = as.numeric(g$estimate[2]),
                          loglik = g$loglik,
                          model = "pp")
  weib.boot.pp <- rbind.data.frame(weib.boot.pp, prms.weib %>% mutate(boot = j))
}

save(weib.boot.pp, file = "Ch1_movement_rates/sims_backup/weib_pp.RData")

# Visualize -----

samp2plot <- 50
weib.boot.pp %>%
  sample_n(., samp2plot) %>%
  ggplot(., aes(x = model, y = est.shape)) +
  geom_violin() +
  # scale_color_manual(values = c("black", mycols)) +
  # geom_boxplot(width = 0.01) +
  # geom_point(color = "grey", alpha = 0.5) +
  # stat_summary(fun.data=mean_sdl, mult=1,
  #              geom="pointrange", color="black") +
  geom_jitter(position = position_jitter(0.1)) +
  labs(title = "Shape") +
  theme_bw() -> p1
# p1

weib.boot.pp %>%
  sample_n(., samp2plot) %>%
  ggplot(., aes(x = model, y = est.scale)) +
  geom_violin() +
  # scale_color_manual(values = c("black", mycols)) +
  #geom_boxplot(width = 0.1) +
  # geom_point() +
  # stat_summary(fun.data=mean_sdl, mult=1,
  #              geom="pointrange", color="black") +
  geom_jitter(position = position_jitter(0.1)) +
  labs(title = "Scale") +
  theme_bw() -> p2

plot_grid(p1, p2)
