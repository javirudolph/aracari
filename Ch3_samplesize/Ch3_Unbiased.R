library(dplyr)
library(ggplot2)



load("Ch3_samplesize/simdata/Bias_df_original.RData")


mbaya %>%
  filter(., run == "run_1") %>%
  dplyr::select(., sampsize, thresh, tru_theta) %>%
  distinct() -> tru_thetas_only

mbaya %>%
  filter(., run == "run_1") %>%
  dplyr::select(., sampsize, thresh, kth_theta) %>%
  group_by(sampsize, thresh) %>%
  summarise(kth_t_mean = mean(kth_theta)) -> mean_theta_is

left_join(tru_thetas_only, mean_theta_is) %>%
  mutate(., unbiased_est = 2*tru_theta - kth_t_mean) -> unbiased_est_df


unbiased_est_df %>%
  filter(., unbiased_est > 0) %>%
  mutate(#sampsize = factor(sampsize),
         thresh = factor(thresh)) %>%
  ggplot() +
  geom_point(aes(x = thresh, y = tru_theta, color = sampsize), size = 3, alpha = 0.5) +
  geom_point(aes(x = thresh, y = unbiased_est, color = sampsize), size = 1) +
  # scale_y_sqrt() +
  scale_y_log10() +
  theme_bw()


library(mgcv)
mod_lm <- gam(unbiased_est ~ sampsize + thresh, data=unbiased_est_df)
summary(mod_lm)
