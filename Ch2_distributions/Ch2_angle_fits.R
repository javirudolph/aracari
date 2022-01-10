# 20211207 Data analysis

# Libraries -----------------------------------------------
set.seed(20211207)

library(aracari)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(fitdistrplus)
library(CircStats)

# Data ----------------------------------------------------
# The script for this data is in the data-raw folder
load("Ch2_distributions/Orig_data_KH/tidy_data.RData")


ptpl %>%
  mutate(individual = paste0("B", id))%>%
  dplyr::select(individual, sgroup, rel.angle, abs.angle) %>%
  drop_na(individual, rel.angle) -> angles_df


# Wrapped Cauchy fit

?dwrpcauchy


# Estimate for complete pooling

mu0 <- circ.mean(angles_df$rel.angle)
rho0 <- est.rho(angles_df$rel.angle)
angle.params <- as.data.frame(wrpcauchy.ml(angles_df$rel.angle, mu0 = mu0, rho0 = rho0)) %>%
  mutate(model = "CP",
         group = NA)


# Partial pooling

sgroups <- paste0("G", 1:7)



for(i in 1:7){

  sub_data <- angles_df %>%
    drop_na(rel.angle) %>%
    filter(., sgroup == sgroups[i])

  mu0 <- circ.mean(sub_data$rel.angle)
  rho0 <- est.rho(sub_data$rel.angle)
  a <- as.data.frame(wrpcauchy.ml(sub_data$rel.angle, mu0 = mu0, rho0 = rho0)) %>%
    mutate(model = "PP",
           group = sgroups[i])

  angle.params <- rbind.data.frame(angle.params, a)

}

# No pooling

individuals <- unique(angles_df$individual)



for(i in 1:15){

  sub_data <- angles_df %>%
    drop_na(rel.angle) %>%
    filter(., individual == individuals[i])

  mu0 <- circ.mean(sub_data$rel.angle)
  rho0 <- est.rho(sub_data$rel.angle)
  a <- as.data.frame(wrpcauchy.ml(sub_data$rel.angle, mu0 = mu0, rho0 = rho0)) %>%
    mutate(model = "NP",
           group = individuals[i])

  angle.params <- rbind.data.frame(angle.params, a)

}



