# Script to run simulations using parameters from aracari movement rates
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
load("Ch1_movement_rates/sims_backup/datagen_np.RData")

load("Ch1_movement_rates/sims_backup/weib_cp.RData")
load("Ch1_movement_rates/sims_backup/weib_pp.RData")
load("Ch1_movement_rates/sims_backup/weib_np.RData")

# Dispersion and Dispersal --------------------------------------------------------------------

cp.summ.df %>%
  dplyr::select(av.disp, dsprsn, popu) %>%
  mutate(model = "cp") %>%
  bind_rows(., pp.summ.df %>%
              dplyr::select(av.disp, dsprsn, model)) %>%
  bind_rows(., np.summ.df %>%
              dplyr::select(av.disp, dsprsn, popu) %>%
              mutate(model = paste0("np_", popu))) -> disp_df

disp_df %>%
  group_by(., model) %>%
  sample_n(., 100) %>%
  ggplot(., aes(y = dsprsn, x = model)) +
  geom_boxplot() +
  geom_point(color = "grey", alpha = 0.5) +
  labs(title = "Dispersion", y = "Dispersion in meters", x = "Model") +
  theme_bw() +
  theme(legend.position = "none") -> p_dispersion

disp_df %>%
  group_by(., model) %>%
  sample_n(., 100) %>%
  ggplot(., aes(y = av.disp, x = model)) +
  geom_boxplot() +
  geom_point(color = "grey", alpha = 0.5) +
  labs(title = "Average dispersal", y = "Dispersal in meters", x = "Model") +
  theme_bw() +
  theme(legend.position = "none") -> p_dispersal

plot_grid(p_dispersal, p_dispersion)

ggsave2(filename = "Ch1_movement_rates/Figures/disp_plot.png")



# Weibull fit ------------------------------------------------------------

weib.boot.cp %>%
  mutate(model = "cp") %>%
  dplyr::select(., est.shape, est.scale, boot, model) %>%
  bind_rows(., weib.boot.pp %>%
              dplyr::select(., est.shape, est.scale, boot, model)) %>%
  bind_rows(., weib.boot.np %>%
              mutate(model = paste0("np_", popu)) %>%
              dplyr::select(., est.shape, est.scale, boot, model)) -> joint_weib

joint_weib %>%
  group_by(model) %>%
  sample_n(., 100) %>%
  ggplot(., aes(x = model, y = est.shape)) +
  geom_violin() +
  # scale_color_manual(values = c("black", mycols)) +
  # geom_boxplot(width = 0.01) +
  # geom_point(color = "grey", alpha = 0.5) +
  # stat_summary(fun.data=mean_sdl, mult=1,
  #              geom="pointrange", color="black") +
  geom_jitter(position = position_jitter(0.1)) +
  labs(title = "Shape") +
  theme_bw() -> weib_shape

joint_weib %>%
  group_by(model) %>%
  sample_n(., 100) %>%
  ggplot(., aes(x = model, y = est.scale)) +
  geom_violin() +
  # scale_color_manual(values = c("black", mycols)) +
  # geom_boxplot(width = 0.01) +
  # geom_point(color = "grey", alpha = 0.5) +
  # stat_summary(fun.data=mean_sdl, mult=1,
  #              geom="pointrange", color="black") +
  geom_jitter(position = position_jitter(0.1)) +
  labs(title = "Scale") +
  theme_bw() -> weib_scale

plot_grid(weib_shape, weib_scale)

ggsave2(filename = "Ch1_movement_rates/Figures/weib_plot.png")
