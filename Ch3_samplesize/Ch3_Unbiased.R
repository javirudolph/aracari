library(dplyr)
library(ggplot2)
library(cowplot)


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




testfit <- unbiased_est_df[, c(1,2,5)]
fits <- predict(mod_lm, newdata=testfit, type='response', se=T)

predicts <- data.frame(testfit, fits) %>%
  mutate(lower = fit - 1.96*se.fit,
         upper = fit + 1.96*se.fit)

ggplot(aes(x=unbiased_est,y=fit), data=predicts) +
  geom_ribbon(aes(ymin = lower, ymax=upper), fill='gray90') +
  geom_line(color='#00aaff') +
  theme_bw()



#################################


# This works and doesn't have repeats.
mbaya %>%
  filter(., run == "run_1") %>%
  mutate(sampsize = factor(sampsize),
         thresh = factor(thresh)) %>%
  group_by(run, sampsize, thresh) %>%
  summarize(bias_hat = (1/B)*sum(kth_theta-tru_theta),
            unbiased1 = tru_theta - bias_hat,
            unbiased2 = (2*tru_theta) - mean(kth_theta)) %>%
  distinct() -> summ_bias

# Visualize the bias in different y scales

summ_bias %>%
  ggplot(., aes(x= thresh, y = bias_hat, color = sampsize)) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none") -> baya1

summ_bias %>%
  ggplot(., aes(x= thresh, y = bias_hat, color = sampsize)) +
  geom_point() +
  scale_y_log10(name = "Log") +
  theme_bw() +
  theme(legend.position = "none") -> baya2

summ_bias %>%
  ggplot(., aes(x= thresh, y = bias_hat, color = sampsize)) +
  geom_point() +
  scale_y_sqrt(name = "sqrt") +
  theme_bw() -> baya3

top_baya <- plot_grid(baya1, baya2, ncol = 2)
plot_grid(top_baya, baya3, nrow = 2)


# Visualize the unbiased estimator

summ_bias %>%
  ggplot(., aes(x= thresh, y = unbiased2, color = sampsize)) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none") -> baya1

summ_bias %>%
  ggplot(., aes(x= thresh, y = unbiased2, color = sampsize)) +
  geom_point() +
  scale_y_log10(name = "Log") +
  theme_bw() +
  theme(legend.position = "none") -> baya2

summ_bias %>%
  ggplot(., aes(x= thresh, y = unbiased2, color = sampsize)) +
  geom_point() +
  scale_y_sqrt(name = "sqrt") +
  theme_bw() -> baya3

top_baya <- plot_grid(baya1, baya2, ncol = 2)
plot_grid(top_baya, baya3, nrow = 2)


# I have a question regarding the zeroes.
# What do I do when we have zeroes? Am I choosing thresholds that are too big? Should I make the threshold choice more continuous?



