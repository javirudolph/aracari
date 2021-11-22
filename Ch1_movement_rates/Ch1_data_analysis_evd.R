# Script for evd analysis

library(dplyr)
library(tidyr)
library(extRemes)
library(ggplot2)
library(cowplot)

# Bring data --------------------------------------------
load("Ch1_movement_rates/sims_backup/datagen_cp.RData")
load("Ch1_movement_rates/sims_backup/datagen_pp.RData")
load("Ch1_movement_rates/sims_backup/datagen_np.RData")


set.seed(563)

# Functions -----------------------------------------------
threshplot_fx <- function(thresh_data, th, r, title = NULL){

  thresh_data %>%
    ggplot(., aes(y = t.scale, x = u.i)) +
    geom_point(shape = 1, size = 2) +
    # geom_line(linetype = "dashed") +
    geom_linerange(aes(x = u.i, ymin = low.t.scale, ymax = up.t.scale)) +
    labs(title = title, x = "Threshold", y = "Reparameterized \n scale") +
    theme_bw() +
    geom_vline(xintercept = th, color = "red", linetype = "dashed") +
    scale_x_continuous(limits = r, n.breaks = 7 ) +
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
    labs(x = "Threshold", y = "Reparameterized \n shape") +
    theme_bw() +
    geom_vline(xintercept = th, color = "red", linetype = "dashed") +
    scale_x_continuous(limits = r, n.breaks = 7 ) +
    # coord_cartesian(xlim = c(200, 550)) +
    NULL -> thresh_shape

  plot_grid(thresh_scale, thresh_shape, nrow = 2)

}

# subset data ---------------------------------

nsamps <- 10000
itera <- 1000


cp.df %>%
  drop_na(., s.id) %>%
  sample_n(., nsamps) -> cp.evd

pp.df %>%
  drop_na(., s.id) %>%
  sample_n(., nsamps) -> pp.evd

np.df %>%
  drop_na(., s.id) %>%
  sample_n(., nsamps) -> np.evd


cp.evd %>%
  mutate(model = toupper(model)) %>%
  bind_rows(pp.evd %>% mutate(model = toupper(model))) %>%
  bind_rows(np.evd %>% mutate(model = toupper(model))) %>%
  group_by(model) %>%
  summarise(q75 = quantile(disp, 0.75),
            q90 = quantile(disp, 0.9),
            q99 = quantile(disp, 0.99),
            max = max(disp))

# Find threshold -------------------------------------
## CP ----------------------------

r <- c(0, 600)

nint <- 100
# r <- quantile(cp.evd$disp, c(0.75, 0.99))
u.i <- matrix(seq(r[1],r[2],,nint), ncol=1)

cp.mrl <- mrlplot(cp.evd$disp)
lm.cp <- lm(cp.mrl[,2]~1)
tval.cp <- confint(lm.cp, level = 0.9)[2]

thresh.cp <- threshrange.plot(cp.evd$disp, r = r, type = "GP", nint = nint) %>%
  as.data.frame() %>%
  mutate(u.i = u.i)

plot.cp <- threshplot_fx(thresh_data = thresh.cp, th = tval.cp, r = r)


## PP ----------------------------

nint <- 100
# r <- quantile(pp.evd$disp, c(0.75, 0.99))
u.i <- matrix(seq(r[1],r[2],,nint), ncol=1)

pp.mrl <- mrlplot(pp.evd$disp)
lm.pp <- lm(pp.mrl[,2]~1)
tval.pp <- confint(lm.pp, level = 0.9)[2]

thresh.pp <- threshrange.plot(pp.evd$disp, r = r, type = "GP", nint = nint) %>%
  as.data.frame() %>%
  mutate(u.i = u.i)

plot.pp <- threshplot_fx(thresh_data = thresh.pp, th = tval.pp, r = r)

## NP ----------------------------

nint <- 100
# r <- quantile(pp.evd$disp, c(0.75, 0.99))
u.i <- matrix(seq(r[1],r[2],,nint), ncol=1)

np.mrl <- mrlplot(np.evd$disp)
lm.np <- lm(np.mrl[,2]~1)
tval.np <- confint(lm.np, level = 0.9)[2]

thresh.np <- threshrange.plot(np.evd$disp, r = r, type = "GP", nint = nint) %>%
  as.data.frame() %>%
  mutate(u.i = u.i)

plot.np <- threshplot_fx(thresh_data = thresh.np, th = tval.pp, r = r)
