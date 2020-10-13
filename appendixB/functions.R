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
    dispInfo <- map_df(oneModel, "seedInfo") %>%
      rename(dispersal = seed_disp)
    # dispersal <- map(dispInfo, "seed_disp") %>%
    #   unlist()

    out <- data.frame(model = models[i], dispInfo)
    df <- rbind.data.frame(df, out)
  }

  return(df)
}
