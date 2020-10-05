

# Start getting seed dispersal simulations

library(tidyverse)
library(fitdistrplus)

directory <- "appendixB/sim10k"


## Parameter setup

nseeds <- 5
nruns <- 10000

source("appendixB/functions.R")

# Start with the population level parameters
# I need to set a seed
set.seed(192)

popind <- readRDS("appendixA/popindprms.RDS")

popind %>%
  mutate(dist_upper = dist,
         dist = str_to_lower(dist),
         dist = str_remove(dist, "onential"),
         dist = str_replace(dist, "lognormal", "lnorm"),
         param = str_to_lower(param)) -> popind

# Population level
popSims <- vector("list", nruns*12)
for(i in 1:length(popSims)){
  popSims[[i]] <- main_simulation(nseeds = nseeds, params = popind %>%
                                    dplyr::filter(data == "POP"))
}

saveRDS(popSims, paste0(directory, "populationSims.RDS"))

#popSims <- readRDS("appendixA/populationSims.RDS")

# These is the data for each seed that got dispersed
# Since there are 5 seeds in each run, with 4 models, each run gets 20 distances
popSeed_data <- get_seed_disp_info(popSims) %>%
  mutate(data = "POP")

ggplot(popSeed_data, aes(x = model, y = dispersal)) +
  geom_boxplot()

popSeed_data %>%
  #dplyr::filter(., dispersal > 500) %>%
  ggplot(., aes(x = dispersal, color = model)) +
  facet_grid(~model) +
  geom_line(stat = "density", lwd = 1) +
  #geom_line(aes(x = dist, group = id), stat = "density", alpha = 0.4) +
  geom_vline(xintercept = 0, color = "grey") +
  #geom_hline(yintercept = 0, color = "grey") +
  #ylim(0, 0.05) +
  labs(title = "Population") +
  #scale_x_log10() +
  theme_classic()


# Individual level simulation ---------------------------------------------

indSim_fx <- function(paramList, individual){
  sims <- vector("list", nruns)
  for(i in 1:length(sims)){
    sims[[i]] <- main_simulation(nseeds = nseeds, params = paramList %>%
                                   dplyr::filter(., data == individual))
  }
  return(sims)
}

IDs <- unique(popind$data)[2:13]

indSeed_data <- NULL
indSims <- vector("list", 12)
for(i in 1:12){
  a <- indSim_fx(paramList = popind, individual = IDs[i])

  b <- get_seed_disp_info(a) %>%
    mutate(data = IDs[i])

  indSeed_data <- rbind(indSeed_data, b)

  indSims[[i]] <- a
}

indSims <- setNames(indSims, IDs)

saveRDS(indSims, paste0(directory, "individualSims.RDS"))

#indSims <- readRDS("appendixA/individualSims.RDS")


seed_dispersal_popind <- rbind.data.frame(popSeed_data, indSeed_data)


# Mixed distribution model for seed dispersal

# These are the best models with BIC
popind %>%
  #filter(., data != "POP") %>%
  filter(., deltaBIC == 0) %>%
  dplyr::select(dist, data) %>%
  distinct() %>%
  mutate(best = paste0(dist,data)) %>%
  ungroup() %>%
  dplyr::pull(., best)-> best_bic

# Find the best with AIC

popind %>%
  #filter(., data != "POP") %>%
  filter(., deltaAIC == 0) %>%
  dplyr::select(dist, data) %>%
  distinct() %>%
  mutate(best = paste0(dist,data)) %>%
  ungroup() %>%
  dplyr::pull(., best)-> best_aic


seed_dispersal_popind %>%
  mutate(dist = str_remove(model, "_move"),
         dist = str_replace(dist, "gam", "gamma"),
         dist = str_replace(dist, "weib", "weibull"),
         combo = paste0(dist, data),
         aic_b = ifelse(combo %in% best_aic, "1", "0"),
         bic_b = ifelse(combo %in% best_bic, "1", "0")) -> with_mixed

saveRDS(with_mixed, file = paste0(directory, "seed_dispersal_popind.RDS"))

#### Family level ---------------------------------------------

popfam <- readRDS("appendixA/popfamprms.RDS")


popfam %>%
  mutate(dist_upper = dist,
         dist = str_to_lower(dist),
         dist = str_remove(dist, "onential"),
         dist = str_replace(dist, "lognormal", "lnorm"),
         param = str_to_lower(param)) -> popfam

# Population level
popSims.fam <- vector("list", nruns*6)
for(i in 1:length(popSims.fam)){
  popSims.fam[[i]] <- main_simulation(nseeds = nseeds, params = popfam %>%
                                    dplyr::filter(data == "POP"))
}

saveRDS(popSims.fam, paste0(directory, "populationSims_family.RDS"))

# These is the data for each seed that got dispersed
# Since there are 5 seeds in each run, with 4 models, each run gets 20 distances
popSeed_data_fam <- get_seed_disp_info(popSims.fam) %>%
  mutate(data = "POP")

ggplot(popSeed_data_fam, aes(x = model, y = dispersal)) +
  geom_boxplot()

popSeed_data_fam %>%
  #dplyr::filter(., dispersal > 500) %>%
  ggplot(., aes(x = dispersal, color = model)) +
  facet_grid(~model) +
  geom_line(stat = "density", lwd = 1) +
  #geom_line(aes(x = dist, group = id), stat = "density", alpha = 0.4) +
  geom_vline(xintercept = 0, color = "grey") +
  #geom_hline(yintercept = 0, color = "grey") +
  #ylim(0, 0.05) +
  labs(title = "Population") +
  #scale_x_log10() +
  theme_classic()


# Family level simulation ---------------------------------------------
famSim_fx <- function(paramList, family_group){
  sims <- vector("list", nruns)
  for(i in 1:length(sims)){
    sims[[i]] <- main_simulation(nseeds = nseeds, params = paramList %>%
                                   dplyr::filter(., data == family_group))
  }
  return(sims)
}

IDs_fam <- unique(popfam$data)[2:7]

famSeed_data <- NULL
famSims <- vector("list", 6)
for(i in 1:6){
  a <- famSim_fx(paramList = popfam, family_group = IDs_fam[i])

  b <- get_seed_disp_info(a) %>%
    mutate(data = IDs_fam[i])

  famSeed_data <- rbind(famSeed_data, b)

  famSims[[i]] <- a
}

famSims <- setNames(famSims, IDs)

saveRDS(famSims, paste0(directory, "familySims.RDS"))


seed_dispersal_popfam <- rbind.data.frame(popSeed_data_fam, famSeed_data)

# Mixed distribution model for seed dispersal

# These are the best models with BIC
popfam %>%
  #filter(., data != "POP") %>%
  filter(., deltaBIC == 0) %>%
  dplyr::select(dist, data) %>%
  distinct() %>%
  mutate(best = paste0(dist,data)) %>%
  ungroup() %>%
  dplyr::pull(., best)-> best_bic

# Find the best with AIC

popfam %>%
  #filter(., data != "POP") %>%
  filter(., deltaAIC == 0) %>%
  dplyr::select(dist, data) %>%
  distinct() %>%
  mutate(best = paste0(dist,data)) %>%
  ungroup() %>%
  dplyr::pull(., best)-> best_aic


seed_dispersal_popfam %>%
  mutate(dist = str_remove(model, "_move"),
         dist = str_replace(dist, "gam", "gamma"),
         dist = str_replace(dist, "weib", "weibull"),
         combo = paste0(dist, data),
         aic_b = ifelse(combo %in% best_aic, "1", "0"),
         bic_b = ifelse(combo %in% best_bic, "1", "0")) -> with_mixed_family

saveRDS(with_mixed_family, file = paste0(directory, "seed_dispersal_popfam.RDS"))

