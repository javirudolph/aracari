
load("paper/simple_exp_runs.RData")


library(fitdistrplus)
library(ggplot2)
library(dplyr)


# Figure for movement rates at the three levels

rates <- data.frame(level = "null", movrate = null_moverate$movrate[1]) %>%
  bind_rows(., indiv_moverate %>%
              rename(level = Bird_ID)) %>%
  bind_rows(., fam_moverate %>%
              rename(level = fam_g)) %>%
  mutate(model = ifelse(level == "null", "Null",
                        ifelse(level %in% paste0("f", 1:7), "Family", "Individual")))

n <- 1000
dens_rates <- NULL
for(i in 1:length(rates$level)) {
  dens_exp <- rexp(n, rate = 1/rates$movrate[i])

  out <- data.frame(level = rates$level[i], model = rates$model[i], value = dens_exp)
  dens_rates <- rbind.data.frame(dens_rates, out)
}

dens_rates %>%
  mutate(model = factor(model, levels = c("Null", "Individual", "Family"))) %>%
  ggplot(., aes(x = value, group = level, color = model)) +
  facet_wrap(~model) +
  geom_density() +
  scale_color_manual(values = c("red", "grey", "black"))
