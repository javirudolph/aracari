# Find thresholds to fit evd

load("paper/dispersal_kernels.RData")

library(dplyr)
library(ggplot2)
library(cowplot)
library(extRemes)

nint <- 100
r <- c(0, 700)
#*****************************************************************************************
# NULL MODEL

# 1. Determine threshold
null_threshplot <- threshrange.plot(null_dispersal$dispersal, r = r, type = "GP", nint = nint)
as.data.frame(null_threshplot) %>%
  mutate(u.i = seq(r[1], r[2], length.out = nint)) -> null_threshplot


null_mrl <- mrlplot(null_dispersal$dispersal, nint = nint)
r_mrl <- range(null_dispersal$dispersal, finite=TRUE)
as.data.frame(null_mrl) %>%
  mutate(u.i_mrl = seq(r_mrl[1], r_mrl[2], length.out = nint),
         slope = `Mean Excess` - lag(`Mean Excess`)) -> null_mrl

# null_thresh <- orig_thresh

# INDIVIDUAL MODEL
# 1. Determine threshold

indiv_threshplot <- threshrange.plot(indiv_dispersal$dispersal, r = r, type = "GP", nint = nint)
as.data.frame(indiv_threshplot) %>%
  mutate(u.i = seq(r[1], r[2], length.out = nint)) -> indiv_threshplot

indiv_mrl <- mrlplot(indiv_dispersal$dispersal, nint = nint)
r_mrl <- range(indiv_dispersal$dispersal, finite=TRUE)
as.data.frame(indiv_mrl) %>%
  mutate(u.i_mrl = seq(r_mrl[1], r_mrl[2], length.out = nint),
         slope = `Mean Excess` - lag(`Mean Excess`)) -> indiv_mrl

# indiv_thresh <- orig_thresh


# FAMILY MODEL
# 1. Determine threshold

fam_threshplot <- threshrange.plot(fam_dispersal$dispersal, r = r, type = "GP", nint = nint)
as.data.frame(fam_threshplot) %>%
  mutate(u.i = seq(r[1], r[2], length.out = nint)) -> fam_threshplot


fam_mrl <- mrlplot(fam_dispersal$dispersal, nint = nint)
r_mrl <- range(fam_dispersal$dispersal, finite=TRUE)
as.data.frame(fam_mrl) %>%
  mutate(u.i_mrl = seq(r_mrl[1], r_mrl[2], length.out = nint),
         slope = `Mean Excess` - lag(`Mean Excess`)) -> fam_mrl

# fam_thresh <- orig_thresh


# Plots

threshplot_fx <- function(thresh_data){

  thresh_data %>%
    ggplot(., aes(y = t.scale, x = u.i)) +
    geom_point(shape = 1, size = 2) +
    # geom_line(linetype = "dashed") +
    geom_linerange(aes(x = u.i, ymin = low.t.scale, ymax = up.t.scale)) +
    labs(x = "Threshold", y = "Reparameterized \n scale") +
    theme_bw() +
    geom_vline(xintercept = th, color = "red", linetype = "dashed") +
    # coord_cartesian(xlim = c(200, 550)) +
    NULL -> thresh_scale

  thresh_data %>%
    ggplot(., aes(y = shape, x = u.i)) +
    geom_point(shape = 1, size = 2) +
    geom_linerange(aes(x = u.i, ymin = low.shape, ymax = up.shape)) +
    labs(x = "Threshold", y = "Reparameterized \n shape") +
    theme_bw() +
    geom_vline(xintercept = th, color = "red", linetype = "dashed") +
    # coord_cartesian(xlim = c(200, 550)) +
    NULL -> thresh_shape

  plot_grid(thresh_scale, thresh_shape, nrow = 2)

}

mrl_plot <- function(mrl_data){
  mrl_data %>%
    ggplot(., aes(x = u.i_mrl, y = `Mean Excess`)) +
    geom_line() +
    geom_line(aes(y = `95% lower`), linetype = "dashed", color = "grey") +
    geom_line(aes(y = `95% upper`), linetype = "dashed", color = "grey") +
    theme_bw() +
    #coord_cartesian(xlim = c(200, 700)) +
    NULL
}

