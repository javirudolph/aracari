# Gut Retention time -----------------------------------------------------
# First step is to generate Gut retention times from gamma distribution to get maximum GRT and that is our simulation time

# Our simulations are with 100 seeds

gen_gut_ret_time <- function(nseeds){
  gut_retention_time <- rgamma(nseeds, shape = 4, scale = 5) + 8
  grt <- data.frame(seedID = 1:nseeds, time = round(gut_retention_time))
  grt
}

# Out of that dataframe, get the maximum gut retention time, as that is the sim_time used
# Basically the bird ran out of seeds and so we stop the simulation


# Distance moved per minute -----------------------------------------------
# Then, we simulate distances moved by the bird

# Need to create a data frame for this function that only feeds the specific estimates for each bird or population, just feed the subset, not the whole data frame.


gen_movement_data <- function(params, sim_time){
  exp_move <- rexp(sim_time, rate = params$estimate[which(params$dist == "exp")])

  gam_move <- rgamma(sim_time, shape = params$estimate[which(params$dist == "gamma") & params$param == "shape"], rate = params$estimate[which(params$dist == "gamma" & params$param == "rate")])

  weib_move <- rweibull(sim_time, shape = params$estimate[which(params$dist == "weibull" & params$param == "shape")], scale = params$estimate[which(params$dist == "weibull" & params$param == "scale")])

  lnorm_move <- rlnorm(sim_time, meanlog = params$estimate[which(params$dist == "lnorm" & params$param == "meanlog")], sdlog = params$estimate[which(params$dist == "lnorm" & params$param == "sdlog")])


  list(exp_move = exp_move, gam_move = gam_move, weib_move = weib_move,
       lnorm_move = lnorm_move)
}




# Generate seed dispersal given distances ---------------------------------
# This function will generate seed dispersal distances given the animal's distance and the simulated gut retention time

gen_seed_dispersal <- function(distance = NULL,
                               grt = NULL){
  sim_time <- max(grt$time)

  angle <- runif(sim_time, min = 0, max = 360)

  distx <- distance*cos(angle)
  xloc <- c(0, cumsum(distx))
  disty <- distance*sin(angle)
  yloc <- c(0, cumsum(disty))
  animalTraj <- data.frame(time = 0:sim_time, xloc = xloc, yloc = yloc)

  seedlocs <- left_join(grt, animalTraj)

  seedlocs$seed_disp <- sqrt(seedlocs$xloc^2 + seedlocs$yloc^2)

  info <- list(animalTraj = animalTraj,
               seedInfo = seedlocs)
  return(info)
}

# Main simulation function
# It takes all the movement models for the animals and generates seed dispersal distances for each, for a given number of simulations


main_simulation <- function(nseeds = NULL,
                            params = NULL){

  grt <- gen_gut_ret_time(nseeds)
  sim_time <- max(grt$time)

  # Generate animal movement distances
  distancesMoved <- gen_movement_data(params, sim_time)

  # Using each of the movement models, generate seed dispersal distances
  dispersal <- map(distancesMoved, gen_seed_dispersal, grt = grt)

}

# All these functions are to manage the output of the simulations:


get_seed_disp_info <- function(simulations){
  models <- c("exp_move", "gam_move", "weib_move", "lnorm_move")

  df <- NULL
  for(i in 1:4){
    oneModel <- map(simulations, models[i])
    dispInfo <- map(oneModel, "seedInfo")
    dispersal <- map(dispInfo, "seed_disp") %>%
      unlist()

    out <- data.frame(model = models[i], dispersal = dispersal)
    df <- rbind.data.frame(df, out)
  }

  return(df)
}


# Start getting seed dispersal simulations

library(tidyverse)
library(fitdistrplus)

# Start with the population level parameters
# I need to set a seed
set.seed(192)

popind <- readRDS("appendixA/popindprms.RDS")


# Population level
popSims <- vector("list", 12000)
for(i in 1:length(popSims)){
  popSims[[i]] <- main_simulation(nseeds = 5, params = popind %>%
                                    dplyr::filter(data == "pop"))
}

saveRDS(popSims, "appendixA/populationSims.RDS")

#popSims <- readRDS("appendixA/populationSims.RDS")

# These is the data for each seed that got dispersed
# Since there are 5 seeds in each run, with 4 models, each run gets 20 distances
popSeed_data <- get_seed_disp_info(popSims) %>%
  mutate(data = "pop")

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
  sims <- vector("list", 1000)
  for(i in 1:length(sims)){
    sims[[i]] <- main_simulation(nseeds = 5, params = paramList %>%
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

saveRDS(indSims, "appendixA/individualSims.RDS")

#indSims <- readRDS("appendixA/individualSims.RDS")


seed_dispersal_popind <- rbind.data.frame(popSeed_data, indSeed_data)

saveRDS(seed_dispersal_popind, file = "appendixA/seed_dispersal_popind.RDS")



# Mixed distribution for individuals
popind %>%
  filter(., data != "pop") %>%
  filter(., deltaBIC == 0)

i_exp <- NULL
for(i in c(1:1000, 3001:5000, 6001:7000, 8001:12000)){
  a <- popSims[[i]]$exp_move$seedInfo$seed_disp

  out <- data.frame(model = "mixed", dispersal = a)
  i_exp <- rbind(i_exp, out)
}

i_lnorm <- NULL
for(i in c(1001:3000, 5001:6000, 7001:8000 )){
  a <- popSims[[i]]$lnorm_move$seedInfo$seed_disp

  out <- data.frame(model = "mixed", dispersal = a)
  i_lnorm <- rbind(i_lnorm, out)
}

mixed_ind_seed_disp <- rbind.data.frame(i_exp, i_lnorm) %>%
  mutate(data = "mixed",
         type = "mixed")

saveRDS(mixed_ind_seed_disp, "appendixA/mixed_ind_seed.RDS")


#### Family level

popfam <- readRDS("appendixA/popfamprms.RDS")


# Population level
popSims.fam <- vector("list", 6000)
for(i in 1:length(popSims.fam)){
  popSims.fam[[i]] <- main_simulation(nseeds = 5, params = popfam %>%
                                    dplyr::filter(data == "pop"))
}

saveRDS(popSims.fam, "appendixA/populationSims_family.RDS")

# These is the data for each seed that got dispersed
# Since there are 5 seeds in each run, with 4 models, each run gets 20 distances
popSeed_data_fam <- get_seed_disp_info(popSims.fam) %>%
  mutate(data = "pop")

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


# Individual level simulation ---------------------------------------------
famSim_fx <- function(paramList, family_group){
  sims <- vector("list", 6000)
  for(i in 1:length(sims)){
    sims[[i]] <- main_simulation(nseeds = 5, params = paramList %>%
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

saveRDS(famSims, "appendixA/familySims.RDS")


seed_dispersal_popfam <- rbind.data.frame(popSeed_data_fam, famSeed_data)

saveRDS(seed_dispersal_popfam, file = "appendixA/seed_dispersal_popfam.RDS")

# Mixed distribution for families

popfam %>%
  filter(., data != "pop") %>%
  filter(., deltaBIC == 0)

f_exp <- NULL
for(i in c(1:2000, 4001:6000)){
  a <- popSims.fam[[i]]$exp_move$seedInfo$seed_disp

  out <- data.frame(model = "mixed", dispersal = a)
  f_exp <- rbind(f_exp, out)
}

f_lnorm <- NULL
for(i in c(2001:4000)){
  a <- popSims.fam[[i]]$lnorm_move$seedInfo$seed_disp

  out <- data.frame(model = "mixed", dispersal = a)
  f_lnorm <- rbind(f_lnorm, out)
}

mixed_fam_seed_disp <- rbind.data.frame(f_exp, f_lnorm) %>%
  mutate(data = "mixed",
         type = "mixed")

saveRDS(mixed_fam_seed_disp, "appendixA/mixed_fam_seed.RDS")



#mixed model, which ones?


# USe the new function or select data from those

