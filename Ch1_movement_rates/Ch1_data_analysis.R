# Analyze data from the Ch1_data_generation.R script
# Includes fitting weibull distributions to seed dispersal distances to estimate seed dispersal kernel
# LIBRARIES --------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(fitdistrplus)
library(Hmisc)

set.seed(98)

# Bring data --------------------------------------------
load("Ch1_movement_rates/sims_backup/datagen_cp.RData")
load("Ch1_movement_rates/sims_backup/datagen_pp.RData")
#load("Ch1_movement_rates/sims_backup/datagen_ppi.RData")
load("Ch1_movement_rates/sims_backup/datagen_np.RData")

# Weibull Seed Dispersal Kernels --------------------------------------------------------
## CP Kernel ------------------------------------------------------------
# sample with replacement groups of 100 seeds and fit weibull
n.boots <- 1000
samp.size <- 100
weib.boot.cp <- NULL

for(j in 1:n.boots){
  s.df <- cp.df %>% drop_na(s.id) %>%
    filter(disp != 0) %>%
    sample_n(., samp.size) %>%
    mutate(disp = round(disp, digits = 2))

  g <- fitdist(s.df$disp, distr = "weibull", method = 'mle', lower = c(0,0))


  prms.weib <- data.frame(est.shape = as.numeric(g$estimate[1]),
                          est.scale = as.numeric(g$estimate[2]),
                          loglik = g$loglik,
                          model = "cp")
  weib.boot.cp <- rbind.data.frame(weib.boot.cp, prms.weib %>% mutate(boot = j))
}

save(weib.boot.cp, file = "Ch1_movement_rates/sims_backup/weib_cp.RData")

## PP Kernel ---------------------------------------------------------------------
n.boots <- 1000
samp.size <- 100
weib.boot.pp <- NULL

for(j in 1:n.boots){
  s.df <- pp.df %>% drop_na(s.id) %>%
    filter(disp != 0) %>%
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

## PPi Kernel -------------------------------------------------------------------
#
# n.boots <- 1000
# samp.size <- 100
# weib.boot.ppi <- NULL
#
# for(j in 1:n.boots){
#   s.df <- ppi.df %>% drop_na(s.id) %>%
#     group_by(model) %>%
#     sample_n(., samp.size) %>%
#     mutate(disp = round(disp, digits = 2))
#
#   # s.df %>%
#   #   ggplot(., aes(x = disp, fill = popu)) +
#   #   geom_histogram()
#
#   # s.df <- df %>% drop_na(s.id)
#
#   weib.fits <- NULL
#
#   for(i in 1:7){
#     dat <- s.df %>%
#       filter(model == paste0("ppi_", i)) %>%
#       filter(disp != 0)
#
#     g <- fitdist(dat$disp, distr = "weibull", method = 'mle', lower = c(0,0))
#     #g <- fitdistr(dat$disp, densfun = "weibull", lower = c(0,0))
#     prms.weib <- data.frame(est.shape = as.numeric(g$estimate[1]),
#                             est.scale = as.numeric(g$estimate[2]),
#                             loglik = g$loglik,
#                             popu = i,
#                             model = "ppi")
#     weib.fits <- rbind.data.frame(weib.fits, prms.weib)
#     # plot(g)
#   }
#
#   weib.boot.ppi <- rbind.data.frame(weib.boot.ppi, weib.fits %>% mutate(boot = j))
# }
#
#
# save(weib.boot.ppi, file = "Ch1_movement_rates/sims_backup/weib_ppi.RData")


## NP Kernel ------------------------------------------------------------------
n.boots <- 1000
samp.size <- 100
weib.boot.np <- NULL

for(j in 1:n.boots){
  s.df <- np.df %>% drop_na(s.id) %>%
    #group_by(popu) %>%
    sample_n(., samp.size) %>%
    mutate(disp = round(disp, digits = 2))

  # s.df %>%
  #   ggplot(., aes(x = disp, fill = popu)) +
  #   geom_histogram()

  # s.df <- df %>% drop_na(s.id)

  weib.fits <- NULL

  for(i in 1:1){
    dat <- s.df %>%
      #filter(popu == i) %>%
      filter(disp != 0)

    g <- fitdist(dat$disp, distr = "weibull", method = 'mle', lower = c(0,0))
    #g <- fitdistr(dat$disp, densfun = "weibull", lower = c(0,0))
    prms.weib <- data.frame(est.shape = as.numeric(g$estimate[1]),
                            est.scale = as.numeric(g$estimate[2]),
                            loglik = g$loglik,
                            popu = i,
                            model = paste0("np_", i))
    weib.fits <- rbind.data.frame(weib.fits, prms.weib)
    # plot(g)
  }

  weib.boot.np <- rbind.data.frame(weib.boot.np, weib.fits %>% mutate(boot = j))
}


save(weib.boot.np, file = "Ch1_movement_rates/sims_backup/weib_np.RData")


# Long distance dispersal percentage ------------------

calc_ldd <- function(df, ...){
  t.seeds <- df %>%
    drop_na(., s.id) %>%
    count()

  n.ldd <- df %>%
    drop_na(., s.id) %>%
    filter(., disp >= 500) %>%
    count()

  value <- n.ldd/t.seeds

  return(value)
}

calc_ldd(cp.df)
calc_ldd(pp.df)
calc_ldd(np.df)



