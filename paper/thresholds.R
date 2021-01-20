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

# Null
null_thresh <- round(null_mrl$u.i_mrl[20])
mrl_plot(null_mrl, title = "Null model mean excess plot") +
  coord_equal() +
  geom_vline(xintercept = null_thresh, color = "red", linetype = "dashed")
threshplot_fx(null_threshplot, null_thresh, title = "Null model threshold plots")

# Indiv
indiv_thresh <- round(indiv_mrl$u.i_mrl[16])
mrl_plot(indiv_mrl) +
  coord_equal()+
  geom_vline(xintercept = indiv_thresh, color = "red", linetype = "dotted")
threshplot_fx(indiv_threshplot, indiv_thresh, title = "Individual model threshold plots")


# Fam
fam_thresh <- round(fam_mrl$u.i_mrl[17])
mrl_plot(fam_mrl) +
  coord_equal()+
  geom_vline(xintercept = fam_thresh, color = "red", linetype = "dotted")
threshplot_fx(fam_threshplot, fam_thresh, title = "Family model threshold plots")


aligned_plots(null_mrl, null_threshplot, null_thresh, title = "Null model")
aligned_plots(indiv_mrl, indiv_threshplot, indiv_thresh, title = "Individual model")
aligned_plots(fam_mrl, fam_threshplot, fam_thresh, title = "Family model")

save.image("paper/thresholds.RData")
