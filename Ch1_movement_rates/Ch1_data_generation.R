# DATA GENERATION
# This script is in charge of simulating all the data for seed dispersal
# We start by importing the real location data from aracaris (Holbrook 2011) and calculate movement rates.
# Dependin gon the pooling scenario considered, we simulate a different number of movement rates

# the functions for the simulation are currently on this script, but that may change in the future is this project ends up as a package.



# Script to run simulations using parameters from aracari movement rates
# LIBRARIES --------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(fitdistrplus)
library(MASS)
library(Hmisc)

set.seed(927)

# Functions --------------------------------------------------------------------------

# Fist we write the function to simulate animal movement only. Plot.it=TRUE will use base R to plot the simulated trajectory.
# prm = parameter used for the movement rate
# This functions simulates movement at every time step by assingning a random angle drawn from a uniform distribution, and a movement distance randomly sampled from an exponential distribution.

sim_movement <- function(prm, t = 1000, plot.it = TRUE, return.data.frame = FALSE){
  tru.rate <- round(prm, 3)
  movedist <- rexp(t, rate = 1/tru.rate)
  angle <- runif(t, min = 0, max = 360)
  distx <- movedist*cos(angle)
  xloc <- c(0, cumsum(distx))
  disty <- movedist*sin(angle)
  yloc <- c(0, cumsum(disty))
  animalTraj <- data.frame(time = 0:t, xloc = xloc, yloc = yloc)
  if(plot.it == TRUE){
    plot(x = xloc, y = yloc, type = "l", main = paste("Rate=",tru.rate))
  }
  if(return.data.frame == TRUE){
    return(animalTraj)
  }
}

# Simulate gut retention time for seeds
# The default parameters are based on the gamma distribution suggested in Morales & Carlo 2006
#
sim_seeds <- function(nseeds = 20, m.prms = NULL,gamma.shape = 4, gamma.scale = 5, gamma.shift = 0,...){
  grt <- round(rgamma(nseeds, shape = gamma.shape, scale = gamma.scale) + gamma.shift)
  t.grt <- max(grt)

  df <- sim_movement(m.prms, t = t.grt, plot.it = FALSE, return.data.frame = TRUE)

  df %>%
    left_join(., data.frame(s.id = 1:nseeds, time = grt), by = "time") %>%
    mutate(disp = sqrt(xloc^2+yloc^2)) -> df

  return(df)
}


# This function takes the simulated data and gives back the summary on seed dispersal
# It includes the average seed location, dispersal per run, and dispersion per run

summ_seeds <- function(df = NULL){
  df %>%
    drop_na(s.id) %>%
    mutate(xi = (mean(xloc)-xloc)^2,
           yi = (mean(yloc)-yloc)^2) %>%
    summarise(x = mean(xloc),
              y = mean(yloc),
              av.disp = mean(disp),
              se.disp = sd(disp)/sqrt(n()),
              dsprsn = sum(sqrt(xi+yi))/n()) -> s.df
  return(s.df)
}


# Load the data ------------------------------------------------------------------------

# point location data

load("data/ptpl.rda")

indiv_moverate <- ptpl %>%
  group_by(Bird_ID, fam_g) %>%
  summarise(movrate = mean(mpm))

fam_moverate <- ptpl %>%
  group_by(fam_g) %>%
  summarise(movrate = mean(mpm),
            sd = sd(mpm))

ids <- ptpl %>% distinct(., Bird_ID, fam_g)

# Gut retention time
load("data/grt_data.rda")

# Estimate parameters --------------------------------------------------------------

# Estimate gut retention time parameters

grtfit <- fitdist(grt_data$grt, "gamma")
gamma.shape <- grtfit$estimate[1]
gamma.scale <- 1/grtfit$estimate[2]

# Estimating parameters for movement rates

logfit <- fitdist(indiv_moverate$movrate, distr = 'lnorm')
logfit.mu <- logfit$estimate[1]
logfit.sigma <- as.numeric(logfit$estimate[2])

x

## Movement rates -----------------------------------------------------
### Complete pooling movrate -------------------
# For the complete pooling scenario we use the expected value of the lognormal distribution we fit to the data

movrate_cp <- as.numeric(exp(logfit$estimate[1] + (logfit$estimate[2]^2)/2))

### Partial pooling movrate -----------------------------------------
# For the partial pooling scenario, since we don't have enough data, we can't fit a lognormal to each family group and then draw movement rates from it. Instead, we use the average movement rate for each social group as the expected value of a lognormal, and use the variance from the complete pooling scenario.

# Estimate the mu's for each family group based on their movement rate.

get_mu <- function(movrate, sigma = logfit.sigma){
  mu <- log(movrate)-((sigma^2)/2)
  return(mu)
}

fam_mus <- get_mu(fam_moverate$movrate)

# Now that each family has a specific mu, we sample individuals from each family
# From field data, we know that social group size is around 6, so we will sample that

movrate_pp <- NULL

for(i in 1:length(fam_mus)){
  i.movrate <- rlnorm(6, meanlog = fam_mus[i], sdlog = logfit.sigma)
  mus.df <- data.frame(fam_g = rep(i, 6), id = paste0("f", i, "_", 1:6), movrate = i.movrate)
  movrate_pp <- rbind.data.frame(movrate_pp, mus.df)
}


### No pooling movrate --------------------------------------------

n.individuals <- 30
movrate_np <- rlnorm(n.individuals, meanlog = logfit.mu, sdlog = logfit.sigma)


## Other Simulation parameters ---------------------------------


# Set folder destination within the project:
if(dir.exists(here::here("Ch1_movement_rates", "sims_backup")) == FALSE){
  dir.create(here::here("Ch1_movement_rates", "sims_backup"))
}

# Function to save data with proper name and date

build.filename <- function(rdata.name){
  rdata.name <- paste0(as.character(rdata.name), ".RData")

  here::here("Ch1_movement_rates", "sims_backup", rdata.name)
}

# number of seeds we give every individual
nseeds <- 5

# Total number of simulation runs we want per scenario (cp, pp, np)

total.kruns <- 100000


# SIMULATED ONLY--------------
## CP Simulations -----------------------------------------------------------------------
# Generate seed dispersal data under a complete pooling scenario

kruns <- total.kruns

cp.df <- NULL
cp.summ.df <- NULL

for(k in 1:kruns){
  a <- sim_seeds(m.prms = movrate_cp, nseeds = nseeds, gamma.shape = gamma.shape, gamma.scale = gamma.scale) %>%
    mutate(run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
           model = "cp")

  b <- summ_seeds(a) %>%
    mutate(run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
           model = "cp")

  cp.df <- rbind.data.frame(cp.df, a)
  cp.summ.df <- rbind.data.frame(cp.summ.df, b)
  #print("cp_run", k)
}

save(cp.df, cp.summ.df, file = build.filename("datagen_cp"))

## PP Simulations -----------------------------------------------------------------------

# number of runs per individual (7 families, each with 6 individuals = 42 individuals)
kruns <- round(total.kruns/42)
nseeds <- 5

pp.df <- NULL
pp.summ.df <- NULL

for(j in 1:nrow(movrate_pp)){
  for(k in 1:kruns){
    a <- sim_seeds(m.prms = movrate_pp$movrate[j], nseeds = nseeds,
                   gamma.shape = gamma.shape, gamma.scale = gamma.scale) %>%
      mutate(id = movrate_pp$id[j],
             run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
             model = "pp")

    b <- summ_seeds(a) %>%
      mutate(id = movrate_pp$id[j],
             run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
             model = "pp")

    pp.df <- rbind.data.frame(pp.df, a)
    pp.summ.df <- rbind.data.frame(pp.summ.df, b)
    #print("pp_run", k, "family_", j)
  }
}

save(pp.df, pp.summ.df, file = build.filename("datagen_pp"))

## NP Simulations -----------------------------------------------

# There's 30 individuals here:
kruns <- round(total.kruns/30)
nseeds <- 5

np.df <- NULL
np.summ.df <- NULL

for(j in 1:length(movrate_np)){
  for(k in 1:kruns){
    a <- sim_seeds(m.prms = movrate_np[j], nseeds = nseeds,
                   gamma.shape = gamma.shape, gamma.scale = gamma.scale) %>%
      mutate(id = paste0("i_", j),
             run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
             model = "pp")

    b <- summ_seeds(a) %>%
      mutate(id = paste0("i_", j),
             run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
             model = "pp")

    np.df <- rbind.data.frame(np.df, a)
    np.summ.df <- rbind.data.frame(np.summ.df, b)
    #print("pp_run", k, "family_", j)
  }
}

save(np.df, np.summ.df, file = build.filename("datagen_np"))



