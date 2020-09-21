# Gut Retention time -----------------------------------------------------
# First step is to generate Gut retention times from gamma distribution to get maximum GRT and that is our simulation time

# Our simulations are with 100 seeds
nseeds <- 100

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

  gam_move <- rgamma(sim_time, shape = params$estimate[which(param$dist == "gamma") & params$param == "shape"], rate = params$estimate[which(params$dist == "gamma" & params$param == "rate")])

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

  # Generate GRT and get the simulation time
  #nseeds <- 10
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

popParams <- readRDS("appendixA/popindprms.RDS")

# we will do 12,000 runs to start
popSims <- vector("list", 12000)
for(i in 1:length(popSims)){
  popSims[[i]] <- main_simulation(nseeds = 5, params = popParams)
}

saveRDS(popSims, "output/populationSims.RDS")

# These is the data for each seed that got dispersed
# Since there are 5 seeds in each run, with 4 models, each run gets 20 distances
popSeed_data <- get_seed_disp_info(popSims)

ggplot(popSeed_data, aes(x = model, y = dispersal)) +
  geom_boxplot()


# Individual level simulation ---------------------------------------------

indParams <- readRDS("output/ind_lev_fit_params.RDS")


indSim_fx <- function(paramList){
  sims <- vector("list", 1000)
  for(i in 1:length(sims)){
    sims[[i]] <- main_simulation(nseeds = 5, params = paramList)
  }
  return(sims)
}


indSims <- indParams %>%
  mutate(simulation = map(pars, indSim_fx),
         dispersal = map(simulation, get_seed_disp_info))

saveRDS(indSims, "output/individualSims.RDS")


indSeed_data <- indSims %>%
  dplyr::select(., Bird_ID, dispersal) %>%
  unnest(cols = c(dispersal))

all_seed_data <- rbind.data.frame(data.frame(popSeed_data, type = "population"),
                                  data.frame(indSeed_data[,2:3], type = "individual"))

saveRDS(all_seed_data, file = "output/all_seed_data.RDS")

# To visualize this, check the Rmd 05.

# Quick workaround with the mixed model output
# Trying this and hopefully it's not wrong.
# Can I use the same information from the simulations before, but instead of taking the data from each individual based on a common distribution fit, I take the data from the best distribution?

# Load the AIC and BIC for each fitting for each individual
# These information criteria things come from the Rmd03
info_criteria_indiv <- readRDS("output/mixed_model_info_crit.RDS")

info_criteria_indiv %>%
  filter(deltaAIC == 0) %>%
  transmute(model = paste0(Bird_ID, str_replace(distribution, "_fit", "_move"))) -> aic_keep

indSeed_data %>%
  mutate(keep = paste0(Bird_ID, model)) %>%
  filter(., keep %in% aic_keep$model) %>%
  mutate(IC = "AIC") -> aic_mixed_distances

info_criteria_indiv %>%
  filter(deltaBIC == 0) %>%
  transmute(model = paste0(Bird_ID, str_replace(distribution, "_fit", "_move"))) -> bic_keep

indSeed_data %>%
  mutate(keep = paste0(Bird_ID, model)) %>%
  filter(., keep %in% bic_keep$model) %>%
  mutate(IC = "BIC") -> bic_mixed_distances

mixed_model_seed_data <- rbind.data.frame(aic_mixed_distances, bic_mixed_distances)

saveRDS(mixed_model_seed_data, file = "output/seed_data_mixed_model.RDS")





