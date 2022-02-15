# Test script

load("Ch2_distributions/Orig_data_KH/tidy_data.RData")

# This is the code included in the raw-data folder, but I don't want to restrict stuff yet.

# ptpl %>%
#   mutate(T_minutes = dt/60,
#          Bird_ID = id,
#          mpm = dist/T_minutes,
#          R2n = lead(R2n)) %>%
#   filter(T_minutes != 0 & T_minutes %in% c(15, 30, 45, 60, 75, 90)) %>%
#   group_by(Bird_ID) %>%
#   add_tally() %>%
#   filter(n >= 30) -> ptpl

# For home range, Kimberly used the cutoff of 38 point observations.

library(dplyr)
library(ggplot2)
library(tidyr)

ptpl <- ptpl %>%
  mutate(T_minutes = dt/60,
         Bird_ID = id,
         mpm = dist/T_minutes,
         R2n = lead(R2n))

id.order <- ptpl %>% dplyr::select(Bird_ID, group) %>% distinct()

ptpl %>% arrange(., group) %>%
  mutate(Bird_ID = factor(Bird_ID, levels = id.order$Bird_ID)) -> ptpl

obs.individual <- ptpl %>%
  group_by(Bird_ID) %>%
  summarise(N_Obs = n()) %>%
  mutate(status = ifelse(N_Obs<30, "exclude",
                         ifelse(N_Obs>40, "include", "edge")))

obs.individual$Bird_ID <- factor(obs.individual$Bird_ID, levels = unique(obs.individual$Bird_ID))

ggplot(obs.individual, aes(x=Bird_ID, y = N_Obs, fill = status)) +
  geom_col() +
  theme_bw() +
  scale_fill_grey(start = 0, end = 0.8) +
  geom_hline(yintercept = 40, color = "red") +
  labs(title = "Number of observations by individual",
       caption = "We use 40 observations as the minimum to consider individuals for our analysis. \n Individual Bird_ID 29 has 38 observations and is right in the edge of the threshold established")


ptpl %>%
  filter(T_minutes != 0) %>%
  group_by(Bird_ID, T_minutes) %>%
  count() %>%
  ggplot(aes(x = as.factor(Bird_ID), y = n, fill = as.factor(T_minutes))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "Individual bird tag",
       y = "Number of intervals recorded") +
  scale_fill_viridis_d(option = "magma", direction = -1) +
  guides(fill = guide_legend(title = "Tracking interval\n in minutes"))

ptpl %>%
  drop_na(mpm) %>%
  filter(mpm != 0) %>%
  group_by(Bird_ID) %>%
  add_tally() %>%
  filter(n >= 30) -> nonzerompm


nonzerompm %>%
  group_by(Bird_ID, T_minutes) %>%
  count() %>%
  ggplot(aes(x = as.factor(Bird_ID), y = n, fill = as.factor(T_minutes))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "Individual bird tag",
       y = "Number of intervals recorded") +
  scale_fill_viridis_d(option = "magma", direction = -1) +
  guides(fill = guide_legend(title = "Tracking interval\n in minutes"))
