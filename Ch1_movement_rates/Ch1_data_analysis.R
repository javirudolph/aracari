# Analyze data from the Ch1_data_generation.R script
# Includes fitting weibull distributions to seed dispersal distances to estimate seed dispersal kernel
# LIBRARIES --------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(fitdistrplus)
library(Hmisc)
library(moments)

set.seed(98)

# Bring data --------------------------------------------
load("Ch1_movement_rates/sims_backup/datagen_cp.RData")
load("Ch1_movement_rates/sims_backup/datagen_pp.RData")
#load("Ch1_movement_rates/sims_backup/datagen_ppi.RData")
load("Ch1_movement_rates/sims_backup/datagen_np.RData")

load("Ch1_movement_rates/sims_backup/datagen_cpr.RData")
load("Ch1_movement_rates/sims_backup/datagen_ppr.RData")
#load("Ch1_movement_rates/sims_backup/datagen_ppi.RData")
load("Ch1_movement_rates/sims_backup/datagen_npr.RData")

# Weibull Seed Dispersal Kernels --------------------------------------------------------
## CP Kernel ------------------------------------------------------------
# sample with replacement groups of 100 seeds and fit weibull
n.boots <- 1000
samp.size <- 100
weib.boot.cp <- NULL

for(j in 1:n.boots){
  s.df <- cp.df %>% drop_na(s.id) %>%
    sample_n(., samp.size) %>%
    mutate(disp = round(disp, digits = 2))%>%
    dplyr::filter(., disp > 0)

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

pp.df %>% drop_na(s.id) %>%
  dplyr::filter(disp == 0) -> pp.df.zeros

for(j in 1:n.boots){
  s.df <- pp.df %>% drop_na(s.id) %>%
    sample_n(., samp.size) %>%
    mutate(disp = round(disp, digits = 2)) %>%
    dplyr::filter(., disp > 0)

  g <- fitdist(s.df$disp, distr = "weibull", method = 'mle', lower = c(0,0))
  print(j)

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

calc_ldd <- function(df, thresh = 500, ...){
  t.seeds <- df %>%
    drop_na(., s.id) %>%
    count()

  n.ldd <- df %>%
    drop_na(., s.id) %>%
    filter(., disp >= thresh) %>%
    count()

  value <- n.ldd/t.seeds

  return(value)
}

cp.seed <- cp.df %>%
  drop_na(., s.id)

pp.seed <- pp.df %>%
  drop_na(., s.id)

np.seed <- np.df %>%
  drop_na(., s.id)



calc_ldd(cp.seed)
calc_ldd(pp.seed)
calc_ldd(np.seed)


disp_table <- data.frame(
  Model = c("CP", "PP", "NP"),
  Mean.dispersal_sd = c(paste0(signif(mean(cp.seed$disp), 4),
                               " (", signif(sd(cp.seed$disp), 2), ")"),
                        paste0(signif(mean(pp.seed$disp), 4),
                               " (", signif(sd(pp.seed$disp), 2), ")"),
                        paste0(signif(mean(np.seed$disp), 4),
                               " (", signif(sd(np.seed$disp), 2), ")")),
  # Seed.dispersion_sd = c(paste0(signif(mean(null_dispersion$seed_dispersion), 4),
  #                               " (", signif(sd(null_dispersion$seed_dispersion), 2), ")"),
  #                        paste0(signif(mean(indiv_dispersion$seed_dispersion), 4),
  #                               " (", signif(sd(indiv_dispersion$seed_dispersion), 2), ")"),
  #                        paste0(signif(mean(fam_dispersion$seed_dispersion), 4),
  #                               " (", signif(sd(fam_dispersion$seed_dispersion), 2), ")")),
  kurtosis = c(signif(kurtosis(cp.seed$disp), 3),
               signif(kurtosis(pp.seed$disp), 3),
               signif(kurtosis(np.seed$disp), 3)),
  Max_dispersal = c(signif(max(cp.seed$disp), 4),
                    signif(max(pp.seed$disp), 4),
                    signif(max(np.seed$disp), 4)),
  percentile_90 = c(signif(quantile(cp.seed$disp, probs = 0.9), 4),
                    signif(quantile(pp.seed$disp, probs = 0.9), 4),
                    signif(quantile(np.seed$disp, probs = 0.9), 4)),
  # Max_dispersal = c(paste0(signif(null_ldd$max_mean, 3), " (", signif(null_ldd$max_sd, 1), ")"),
  #                   paste0(signif(indiv_ldd$max_mean, 3), " (", signif(indiv_ldd$max_sd, 1), ")"),
  #                   paste0(signif(fam_ldd$max_mean, 3), " (", signif(fam_ldd$max_sd, 1), ")")),
  LDD_500 = c(paste0(signif(calc_ldd(cp.seed)*100, 4), "%"),
          paste0(signif(calc_ldd(pp.seed)*100, 4), "%"),
          paste0(signif(calc_ldd(np.seed)*100, 4), "%")),
  LDD_150 = c(paste0(signif(calc_ldd(cp.seed, thresh = 100)*100, 4), "%"),
            paste0(signif(calc_ldd(pp.seed, thresh = 100)*100, 4), "%"),
            paste0(signif(calc_ldd(np.seed, thresh = 100)*100, 4), "%")))

write.csv(disp_table, "Ch1_movement_rates/Figures/table2.csv")

# weibull_table <- data.frame(
#   Weibull_Shape = c(paste0(signif(null_weibull$estimate[1], 4), " (", signif(null_weibull$sd[1], 2), ")"),
#                     paste0(signif(indiv_weibull$estimate[1], 4), " (", signif(indiv_weibull$sd[1], 2), ")"),
#                     paste0(signif(fam_weibull$estimate[1], 4), " (", signif(fam_weibull$sd[1], 2), ")")),
#   Weibull_Scale = c(paste0(signif(null_weibull$estimate[2], 4), " (", signif(null_weibull$sd[2], 2), ")"),
#                     paste0(signif(indiv_weibull$estimate[2], 4), " (", signif(indiv_weibull$sd[2], 2), ")"),
#                     paste0(signif(fam_weibull$estimate[2], 4), " (", signif(fam_weibull$sd[2], 2), ")"))
# )


# EVD -------------------------------------------------------------

library(extRemes)

nint <- 100
r <- c(0, 700)
nsamps <- 1000
#*****************************************************************************************
# CP MODEL

# 1. Determine threshold
cp.df %>%
  drop_na(., s.id) %>%
  sample_n(., nsamps) -> cp.evd


cp_threshplot <- threshrange.plot(cp.evd$disp, r = r, type = "GP", nint = nint)
as.data.frame(cp_threshplot) %>%
  mutate(u.i = seq(r[1], r[2], length.out = nint)) -> cp_threshplot


cp_mrl <- mrlplot(cp.evd$disp, nint = nint)
r_mrl <- range(cp.evd$disp, finite=TRUE)
as.data.frame(cp_mrl) %>%
  mutate(u.i_mrl = seq(r_mrl[1], r_mrl[2], length.out = nint),
         slope = `Mean Excess` - lag(`Mean Excess`)) -> cp_mrl

# PP MODEL
# 1. Determine threshold

pp.df %>%
  drop_na(., s.id) %>%
  sample_n(., nsamps) -> pp.evd

pp_threshplot <- threshrange.plot(pp.evd$disp, r = r, type = "GP", nint = nint)
as.data.frame(pp_threshplot) %>%
  mutate(u.i = seq(r[1], r[2], length.out = nint)) -> pp_threshplot

pp_mrl <- mrlplot(pp.evd$disp, nint = nint)
r_mrl <- range(pp.evd$disp, finite=TRUE)
as.data.frame(pp_mrl) %>%
  mutate(u.i_mrl = seq(r_mrl[1], r_mrl[2], length.out = nint),
         slope = `Mean Excess` - lag(`Mean Excess`)) -> pp_mrl



# NP MODEL
# 1. Determine threshold

np.df %>%
  drop_na(., s.id) %>%
  sample_n(., nsamps) -> np.evd

np_threshplot <- threshrange.plot(np.evd$disp, r = r, type = "GP", nint = nint)
as.data.frame(np_threshplot) %>%
  mutate(u.i = seq(r[1], r[2], length.out = nint)) -> np_threshplot


np_mrl <- mrlplot(np.evd$disp, nint = nint)
r_mrl <- range(np.evd$disp, finite=TRUE)
as.data.frame(np_mrl) %>%
  mutate(u.i_mrl = seq(r_mrl[1], r_mrl[2], length.out = nint),
         slope = `Mean Excess` - lag(`Mean Excess`)) -> np_mrl



# Plots

threshplot_fx <- function(thresh_data, th, title = NULL){

  thresh_data %>%
    ggplot(., aes(y = t.scale, x = u.i)) +
    geom_point(shape = 1, size = 2) +
    # geom_line(linetype = "dashed") +
    geom_linerange(aes(x = u.i, ymin = low.t.scale, ymax = up.t.scale)) +
    labs(title = title, x = "Threshold", y = "Scale") +
    theme_bw() +
    geom_vline(xintercept = th, color = "red", linetype = "dashed") +
    scale_x_continuous(limits = c(0, 700), n.breaks = 7 ) +
    # scale_x_continuous(name = NULL, labels = NULL, breaks = NULL) +
    # theme(axis.title.x = element_blank(),
    #       axis.ticks.x = element_blank(),
    #       axis.text.x = element_blank()) +
    # coord_cartesian(xlim = c(200, 550)) +
    NULL -> thresh_scale

  thresh_data %>%
    ggplot(., aes(y = shape, x = u.i)) +
    geom_point(shape = 1, size = 2) +
    geom_linerange(aes(x = u.i, ymin = low.shape, ymax = up.shape)) +
    labs(x = "Threshold", y = "Shape") +
    theme_bw() +
    geom_vline(xintercept = th, color = "red", linetype = "dashed") +
    scale_x_continuous(limits = c(0, 700), n.breaks = 7 ) +
    # coord_cartesian(xlim = c(200, 550)) +
    NULL -> thresh_shape

  plot_grid(thresh_scale, thresh_shape, nrow = 2)

}

mrl_plot <- function(mrl_data, title = NULL){

  mrl_data %>%
    ggplot(., aes(x = u.i_mrl, y = `Mean Excess`)) +
    geom_line() +
    geom_line(aes(y = `95% lower`), linetype = "dashed", color = "black") +
    geom_line(aes(y = `95% upper`), linetype = "dashed", color = "black") +
    labs(title = title, x = "Threshold values") +
    scale_x_continuous(n.breaks = ceiling(round(max(mrl_data$u.i_mrl))/100)) +
    theme_bw()
}

aligned_plots <- function(mrl_data, thresh_data, threshold, title = NULL){
  title <- ggdraw() +
    draw_label(
      title,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )


  plot_grid(mrl_plot(mrl_data) +
              #coord_equal() +
              geom_vline(xintercept = threshold, color = "red", linetype = "dashed"),
            threshplot_fx(thresh_data, threshold),
            nrow = 2,
            rel_heights = c(2,3),
            labels = "auto") -> actual_plots

  plot_grid(title,
            actual_plots,
            ncol = 1,
            rel_heights = c(0.1,1))
}

# CP
cp_thresh <- round(cp_mrl$u.i_mrl[15])
mrl_plot(cp_mrl, title = "CP mean excess plot") +
  coord_equal() +
  geom_vline(xintercept = cp_thresh, color = "red", linetype = "dashed")
threshplot_fx(cp_threshplot, cp_thresh, title = "CP threshold plots")

# PP
pp_thresh <- round(pp_mrl$u.i_mrl[15])
mrl_plot(pp_mrl) +
  coord_equal()+
  geom_vline(xintercept = pp_thresh, color = "red", linetype = "dotted")
threshplot_fx(pp_threshplot, pp_thresh, title = "PP threshold plots")


# NP
np_thresh <- round(np_mrl$u.i_mrl[15])
mrl_plot(np_mrl) +
  coord_equal()+
  geom_vline(xintercept = np_thresh, color = "red", linetype = "dotted")
threshplot_fx(np_threshplot, np_thresh, title = "NP threshold plots")


aligned_plots(cp_mrl, cp_threshplot, cp_thresh, title = "CP model")
aligned_plots(pp_mrl, pp_threshplot, pp_thresh, title = "PP model")
aligned_plots(np_mrl, np_threshplot, np_thresh, title = "NP model")

# REAL data ----------------------

# Weibull Seed Dispersal Kernels --------------------------------------------------------
## CPr Kernel ------------------------------------------------------------
# sample with replacement groups of 100 seeds and fit weibull
n.boots <- 1000
samp.size <- 100
weib.boot.cpr <- NULL

for(j in 1:n.boots){
  s.df <- cpr.df %>% drop_na(s.id) %>%
    sample_n(., samp.size) %>%
    mutate(disp = round(disp, digits = 2))%>%
    dplyr::filter(., disp > 0)

  g <- fitdist(s.df$disp, distr = "weibull", method = 'mle', lower = c(0,0))


  prms.weib <- data.frame(est.shape = as.numeric(g$estimate[1]),
                          est.scale = as.numeric(g$estimate[2]),
                          loglik = g$loglik,
                          model = "cpr")
  weib.boot.cpr <- rbind.data.frame(weib.boot.cpr, prms.weib %>% mutate(boot = j))
}

save(weib.boot.cpr, file = "Ch1_movement_rates/sims_backup/weib_cpr.RData")

## PPr Kernel ---------------------------------------------------------------------
n.boots <- 1000
samp.size <- 100
weib.boot.ppr <- NULL

ppr.df %>% drop_na(s.id) %>%
  dplyr::filter(disp == 0) -> ppr.df.zeros

for(j in 1:n.boots){
  s.df <- ppr.df %>% drop_na(s.id) %>%
    sample_n(., samp.size) %>%
    mutate(disp = round(disp, digits = 2)) %>%
    dplyr::filter(., disp > 0)

  g <- fitdist(s.df$disp, distr = "weibull", method = 'mle', lower = c(0,0))
  print(j)

  prms.weib <- data.frame(est.shape = as.numeric(g$estimate[1]),
                          est.scale = as.numeric(g$estimate[2]),
                          loglik = g$loglik,
                          model = "ppr")
  weib.boot.ppr <- rbind.data.frame(weib.boot.ppr, prms.weib %>% mutate(boot = j))
}

save(weib.boot.ppr, file = "Ch1_movement_rates/sims_backup/weib_ppr.RData")


## NPr Kernel ------------------------------------------------------------------
n.boots <- 1000
samp.size <- 100
weib.boot.npr <- NULL

for(j in 1:n.boots){
  s.df <- npr.df %>% drop_na(s.id) %>%
    #group_by(popu) %>%
    sample_n(., samp.size) %>%
    mutate(disp = round(disp, digits = 2))%>%
    dplyr::filter(., disp > 0)

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
                            model = "npr")
    weib.fits <- rbind.data.frame(weib.fits, prms.weib)
    # plot(g)
  }

  weib.boot.npr <- rbind.data.frame(weib.boot.npr, weib.fits %>% mutate(boot = j))
}


save(weib.boot.npr, file = "Ch1_movement_rates/sims_backup/weib_npr.RData")


# Long distance dispersal percentage ------------------

calc_ldd <- function(df, thresh = 500, ...){
  t.seeds <- df %>%
    drop_na(., s.id) %>%
    count()

  n.ldd <- df %>%
    drop_na(., s.id) %>%
    filter(., disp >= thresh) %>%
    count()

  value <- n.ldd/t.seeds

  return(value)
}

cpr.seed <- cpr.df %>%
  drop_na(., s.id)

ppr.seed <- ppr.df %>%
  drop_na(., s.id)

npr.seed <- npr.df %>%
  drop_na(., s.id)



calc_ldd(cpr.seed)
calc_ldd(ppr.seed)
calc_ldd(npr.seed)


disp_table_r <- data.frame(
  Model = c("CPr", "PPr", "NPr"),
  Mean.dispersal_sd = c(paste0(signif(mean(cpr.seed$disp), 4),
                               " (", signif(sd(cpr.seed$disp), 2), ")"),
                        paste0(signif(mean(ppr.seed$disp), 4),
                               " (", signif(sd(ppr.seed$disp), 2), ")"),
                        paste0(signif(mean(npr.seed$disp), 4),
                               " (", signif(sd(npr.seed$disp), 2), ")")),
  # Seed.dispersion_sd = c(paste0(signif(mean(null_dispersion$seed_dispersion), 4),
  #                               " (", signif(sd(null_dispersion$seed_dispersion), 2), ")"),
  #                        paste0(signif(mean(indiv_dispersion$seed_dispersion), 4),
  #                               " (", signif(sd(indiv_dispersion$seed_dispersion), 2), ")"),
  #                        paste0(signif(mean(fam_dispersion$seed_dispersion), 4),
  #                               " (", signif(sd(fam_dispersion$seed_dispersion), 2), ")")),
  kurtosis = c(signif(kurtosis(cpr.seed$disp), 3),
               signif(kurtosis(ppr.seed$disp), 3),
               signif(kurtosis(npr.seed$disp), 3)),
  Max_dispersal = c(signif(max(cpr.seed$disp), 4),
                    signif(max(ppr.seed$disp), 4),
                    signif(max(npr.seed$disp), 4)),
  percentile_90 = c(signif(quantile(cpr.seed$disp, probs = 0.9), 4),
                    signif(quantile(ppr.seed$disp, probs = 0.9), 4),
                    signif(quantile(npr.seed$disp, probs = 0.9), 4)),
  # Max_dispersal = c(paste0(signif(null_ldd$max_mean, 3), " (", signif(null_ldd$max_sd, 1), ")"),
  #                   paste0(signif(indiv_ldd$max_mean, 3), " (", signif(indiv_ldd$max_sd, 1), ")"),
  #                   paste0(signif(fam_ldd$max_mean, 3), " (", signif(fam_ldd$max_sd, 1), ")")),
  LDD_500 = c(paste0(signif(calc_ldd(cpr.seed)*100, 4), "%"),
              paste0(signif(calc_ldd(ppr.seed)*100, 4), "%"),
              paste0(signif(calc_ldd(npr.seed)*100, 4), "%")),
  LDD_150 = c(paste0(signif(calc_ldd(cpr.seed, thresh = 100)*100, 4), "%"),
              paste0(signif(calc_ldd(ppr.seed, thresh = 100)*100, 4), "%"),
              paste0(signif(calc_ldd(npr.seed, thresh = 100)*100, 4), "%")))

write.csv(disp_table_r, "Ch1_movement_rates/Figures/table2r.csv")

# weibull_table <- data.frame(
#   Weibull_Shape = c(paste0(signif(null_weibull$estimate[1], 4), " (", signif(null_weibull$sd[1], 2), ")"),
#                     paste0(signif(indiv_weibull$estimate[1], 4), " (", signif(indiv_weibull$sd[1], 2), ")"),
#                     paste0(signif(fam_weibull$estimate[1], 4), " (", signif(fam_weibull$sd[1], 2), ")")),
#   Weibull_Scale = c(paste0(signif(null_weibull$estimate[2], 4), " (", signif(null_weibull$sd[2], 2), ")"),
#                     paste0(signif(indiv_weibull$estimate[2], 4), " (", signif(indiv_weibull$sd[2], 2), ")"),
#                     paste0(signif(fam_weibull$estimate[2], 4), " (", signif(fam_weibull$sd[2], 2), ")"))
# )


# EVD -------------------------------------------------------------

library(extRemes)

nint <- 100
r <- c(0, 700)
nsamps <- 1000
#*****************************************************************************************
# CPr MODEL

# 1. Determine threshold
cpr.df %>%
  drop_na(., s.id) %>%
  sample_n(., nsamps) -> cpr.evd


cpr_threshplot <- threshrange.plot(cpr.evd$disp, r = r, type = "GP", nint = nint)
as.data.frame(cpr_threshplot) %>%
  mutate(u.i = seq(r[1], r[2], length.out = nint)) -> cpr_threshplot


cpr_mrl <- mrlplot(cpr.evd$disp, nint = nint)
r_mrl <- range(cpr.evd$disp, finite=TRUE)
as.data.frame(cpr_mrl) %>%
  mutate(u.i_mrl = seq(r_mrl[1], r_mrl[2], length.out = nint),
         slope = `Mean Excess` - lag(`Mean Excess`)) -> cpr_mrl

# PPr MODEL
# 1. Determine threshold

ppr.df %>%
  drop_na(., s.id) %>%
  sample_n(., nsamps) -> ppr.evd

ppr_threshplot <- threshrange.plot(ppr.evd$disp, r = r, type = "GP", nint = nint)
as.data.frame(ppr_threshplot) %>%
  mutate(u.i = seq(r[1], r[2], length.out = nint)) -> ppr_threshplot

ppr_mrl <- mrlplot(ppr.evd$disp, nint = nint)
r_mrl <- range(ppr.evd$disp, finite=TRUE)
as.data.frame(ppr_mrl) %>%
  mutate(u.i_mrl = seq(r_mrl[1], r_mrl[2], length.out = nint),
         slope = `Mean Excess` - lag(`Mean Excess`)) -> ppr_mrl



# NPr MODEL
# 1. Determine threshold

npr.df %>%
  drop_na(., s.id) %>%
  sample_n(., nsamps) -> npr.evd

npr_threshplot <- threshrange.plot(npr.evd$disp, r = r, type = "GP", nint = nint)
as.data.frame(npr_threshplot) %>%
  mutate(u.i = seq(r[1], r[2], length.out = nint)) -> npr_threshplot


npr_mrl <- mrlplot(npr.evd$disp, nint = nint)
r_mrl <- range(npr.evd$disp, finite=TRUE)
as.data.frame(npr_mrl) %>%
  mutate(u.i_mrl = seq(r_mrl[1], r_mrl[2], length.out = nint),
         slope = `Mean Excess` - lag(`Mean Excess`)) -> npr_mrl



# Plots

threshplot_fx <- function(thresh_data, th, title = NULL){

  thresh_data %>%
    ggplot(., aes(y = t.scale, x = u.i)) +
    geom_point(shape = 1, size = 2) +
    # geom_line(linetype = "dashed") +
    geom_linerange(aes(x = u.i, ymin = low.t.scale, ymax = up.t.scale)) +
    labs(title = title, x = "Threshold", y = "Scale") +
    theme_bw() +
    geom_vline(xintercept = th, color = "red", linetype = "dashed") +
    scale_x_continuous(limits = c(0, 700), n.breaks = 7 ) +
    # scale_x_continuous(name = NULL, labels = NULL, breaks = NULL) +
    # theme(axis.title.x = element_blank(),
    #       axis.ticks.x = element_blank(),
    #       axis.text.x = element_blank()) +
    # coord_cartesian(xlim = c(200, 550)) +
    NULL -> thresh_scale

  thresh_data %>%
    ggplot(., aes(y = shape, x = u.i)) +
    geom_point(shape = 1, size = 2) +
    geom_linerange(aes(x = u.i, ymin = low.shape, ymax = up.shape)) +
    labs(x = "Threshold", y = "Shape") +
    theme_bw() +
    geom_vline(xintercept = th, color = "red", linetype = "dashed") +
    scale_x_continuous(limits = c(0, 700), n.breaks = 7 ) +
    # coord_cartesian(xlim = c(200, 550)) +
    NULL -> thresh_shape

  plot_grid(thresh_scale, thresh_shape, nrow = 2)

}

mrl_plot <- function(mrl_data, title = NULL){

  mrl_data %>%
    ggplot(., aes(x = u.i_mrl, y = `Mean Excess`)) +
    geom_line() +
    geom_line(aes(y = `95% lower`), linetype = "dashed", color = "black") +
    geom_line(aes(y = `95% upper`), linetype = "dashed", color = "black") +
    labs(title = title, x = "Threshold values") +
    scale_x_continuous(n.breaks = ceiling(round(max(mrl_data$u.i_mrl))/100)) +
    theme_bw()
}

aligned_plots <- function(mrl_data, thresh_data, threshold, title = NULL){
  title <- ggdraw() +
    draw_label(
      title,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )


  plot_grid(mrl_plot(mrl_data) +
              #coord_equal() +
              geom_vline(xintercept = threshold, color = "red", linetype = "dashed"),
            threshplot_fx(thresh_data, threshold),
            nrow = 2,
            rel_heights = c(2,3),
            labels = "auto") -> actual_plots

  plot_grid(title,
            actual_plots,
            ncol = 1,
            rel_heights = c(0.1,1))
}

# CPr
cpr_thresh <- round(cpr_mrl$u.i_mrl[15])
mrl_plot(cpr_mrl, title = "CPr mean excess plot") +
  coord_equal() +
  geom_vline(xintercept = cpr_thresh, color = "red", linetype = "dashed")
threshplot_fx(cpr_threshplot, cpr_thresh, title = "CPr threshold plots")

# PPr
ppr_thresh <- round(ppr_mrl$u.i_mrl[15])
mrl_plot(ppr_mrl) +
  coord_equal()+
  geom_vline(xintercept = ppr_thresh, color = "red", linetype = "dotted")
threshplot_fx(ppr_threshplot, ppr_thresh, title = "PPr threshold plots")


# NPr
npr_thresh <- round(npr_mrl$u.i_mrl[15])
mrl_plot(npr_mrl) +
  coord_equal()+
  geom_vline(xintercept = npr_thresh, color = "red", linetype = "dotted")
threshplot_fx(npr_threshplot, npr_thresh, title = "NPr threshold plots")


aligned_plots(cpr_mrl, cpr_threshplot, cp_thresh, title = "CP model")
aligned_plots(ppr_mrl, ppr_threshplot, pp_thresh, title = "PP model")
aligned_plots(npr_mrl, npr_threshplot, np_thresh, title = "NP model")
