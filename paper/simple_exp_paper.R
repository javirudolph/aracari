

# how many seeds each bird gets on a simulation run, which is equivalent to how many fruits in a foraging session in the same tree

# nseeds <- round(runif(1, 3, 6))
nseeds <- 5

# Sample a random Gut Retention Time for each one of the seeds
grt_tt <- round(rgamma(nseeds, shape = 4, scale = 5) + 8)

# Maximum gut retention time determines the simulation time
tt <- max(grt_tt)

movrate <- 28

# Sample movement distances from an exponential distribution. These are the lengths that the animal moves at each time step in minutes.
movedist <- rexp(tt, 1/movrate)

# Movement algorithm for this random walk. Sample a random angle of movement from a uniform distribution. Determine the movement in x and y, the animal's trajectory is the cumulative sum of those. We are not imposing any limits to the area of movement.
angle <- runif(tt, min = 0, max = 360)
distx <- movedist*cos(angle)
xloc <- c(0, cumsum(distx))
disty <- movedist*sin(angle)
yloc <- c(0, cumsum(disty))
animalTraj <- data.frame(time = 0:tt, xloc = xloc, yloc = yloc)

# Determine where the seeds are dropped based on the animal's location.
seed_loc <- data.frame(seed_id = 1:nseeds, time = grt_tt)
seed_loc <- merge(animalTraj, seed_loc, by = "time")
seed_loc$dispersal <- sqrt(seed_loc$xloc^2 + seed_loc$yloc^2)

outsim <- list(lengths = movedist,
               movement = animalTraj,
              seed = seed_loc)

# Calculate the mean location of seeds
x_m <- mean(outsim$seed$xloc)
y_m <- mean(outsim$seed$yloc)
mean_dispersal <- mean(outsim$seed$dispersal)
se_dispersal <- sd(outsim$seed$dispersal)/sqrt(length(outsim$seed$dispersal))

# Calculate seed dispersion
xi <- (x_m - outsim$seed$xloc)^2
yi <- (y_m - outsim$seed$yloc)^2
seed_dispersion <- sum(sqrt(xi + yi))/length(outsim$seed$time)


par(mfrow = c(1,3))
hist(movedist, main = "MD")
plot(x = xloc, y = yloc, type = "l", main = "Movement")
points(x = seed_loc$xloc, y = seed_loc$yloc, col = "blue", pch = 16)
hist(seed_loc$dispersal, main = "Dispersal")
abline(v = mean_dispersal, col = "red")


###### Make function to run multiple simulations of this


sim_run <- function(movrate){

  nseeds <- 5
  grt_tt <- round(rgamma(nseeds, shape = 4, scale = 5) + 8)
  tt <- max(grt_tt)

  movedist <- rexp(tt, 1/movrate)

  angle <- runif(tt, min = 0, max = 360)
  distx <- movedist*cos(angle)
  xloc <- c(0, cumsum(distx))
  disty <- movedist*sin(angle)
  yloc <- c(0, cumsum(disty))
  animalTraj <- data.frame(time = 0:tt, xloc = xloc, yloc = yloc)

  seed_loc <- data.frame(seed_id = 1:nseeds, time = grt_tt)
  seed_loc <- merge(animalTraj, seed_loc, by = "time")
  seed_loc$dispersal <- sqrt(seed_loc$xloc^2 + seed_loc$yloc^2)

  x_m <- mean(seed_loc$xloc)
  y_m <- mean(seed_loc$yloc)
  mean_dispersal <- mean(seed_loc$dispersal)
  se_dispersal <- sd(seed_loc$dispersal)/sqrt(length(seed_loc$dispersal))

  # Calculate seed dispersion
  xi <- (x_m - seed_loc$xloc)^2
  yi <- (y_m - seed_loc$yloc)^2
  seed_dispersion <- sum(sqrt(xi + yi))/length(seed_loc$time)

  outsim <- list(lengths = movedist,
                 movement = animalTraj,
                 seed = seed_loc,
                 mean_dispersal = mean_dispersal,
                 se_dispersal = se_dispersal,
                 seed_dispersion = seed_dispersion)

  return(outsim)
}


onerun <- sim_run(28)

############# Simulations for the population

library(dplyr)
library(purrr)
library(ggplot2)
library(purrr)
library(cowplot)

set.seed(27)

data(ptpl)

null_moverate <- data.frame(Bird_ID = unique(ptpl$Bird_ID), movrate = mean(ptpl$mpm))

indiv_moverate <- ptpl %>%
  group_by(Bird_ID) %>%
  summarise(movrate = mean(mpm))

fam_moverate <- ptpl %>%
  group_by(fam_g) %>%
  summarise(movrate = mean(mpm))

ids <- ptpl %>% distinct(., Bird_ID, fam_g)

#### Null model
# Assume all individuals share the same movement rate, which is the average movement rate for all data points

# Run 10,000 simulation runs per individual or family group

nruns <- 10000

null_dispersal <- NULL
null_dispersion <- NULL
for(j in 1:length(null_moverate$Bird_ID)){
  moverate <- null_moverate$movrate[j]

  manysims <- vector("list", nruns)
  for(i in 1:length(manysims)){
    manysims[[i]] <- sim_run(moverate)
    manysims[[i]]$sim_run <- paste("sim_", i)
  }

  dispersal <- map_df(manysims, "seed")$dispersal
  id <- null_moverate$Bird_ID[j]
  mean_dispersal <- map_dbl(manysims, "mean_dispersal")
  se_dispersal <- map_dbl(manysims, "se_dispersal")
  seed_dispersion <- map_dbl(manysims, "seed_dispersion")

  out <- data.frame(dispersal = dispersal, id = id)
  null_dispersal <- rbind.data.frame(null_dispersal, out)

  out2 <- data.frame(id = id,
                     mean_dispersal = mean_dispersal,
                     se_dispersal = se_dispersal,
                     seed_dispersion = seed_dispersion)
  null_dispersion <- rbind.data.frame(null_dispersion, out2)

}

# Plots
plot_grid(
  null_dispersal %>%
    ggplot(., aes(x = dispersal)) +
    geom_histogram(),
  plot_grid(null_dispersion %>%
            ggplot(., aes(y = mean_dispersal)) +
            geom_boxplot(),
          null_dispersion %>%
            ggplot(., aes(y = seed_dispersion)) +
            geom_boxplot()), nrow = 2)


### Individual simulation

indiv_dispersal <- NULL
indiv_dispersion <- NULL
for(j in 1:length(indiv_moverate$Bird_ID)){
  moverate <- indiv_moverate$movrate[j]

  manysims <- vector("list", nruns)
  for(i in 1:length(manysims)){
    manysims[[i]] <- sim_run(moverate)
    manysims[[i]]$sim_run <- paste("sim_", i)
  }

  dispersal <- map_df(manysims, "seed")$dispersal
  id <- indiv_moverate$Bird_ID[j]
  mean_dispersal <- map_dbl(manysims, "mean_dispersal")
  se_dispersal <- map_dbl(manysims, "se_dispersal")
  seed_dispersion <- map_dbl(manysims, "seed_dispersion")

  out <- data.frame(dispersal = dispersal, id = id)
  indiv_dispersal <- rbind.data.frame(indiv_dispersal, out)

  out2 <- data.frame(id = id,
                     mean_dispersal = mean_dispersal,
                     se_dispersal = se_dispersal,
                     seed_dispersion = seed_dispersion)
  indiv_dispersion <- rbind.data.frame(indiv_dispersion, out2)

}

# plots
plot_grid(indiv_dispersal %>%
            ggplot(., aes(x = dispersal)) +
            geom_histogram(),
          plot_grid(indiv_dispersion %>%
                      ggplot(., aes(y = mean_dispersal)) +
                      geom_boxplot(),
                    indiv_dispersion %>%
                      ggplot(., aes(y = seed_dispersion)) +
                      geom_boxplot()), nrow = 2)



### Social or Family group simulation

fam_dispersal <- NULL
fam_dispersion <- NULL
for(j in 1:length(fam_moverate$fam_g)){
  moverate <- fam_moverate$movrate[j]

  manysims <- vector("list", nruns)
  for(i in 1:length(manysims)){
    manysims[[i]] <- sim_run(moverate)
    manysims[[i]]$sim_run <- paste("sim_", i)
  }

  dispersal <- map_df(manysims, "seed")$dispersal
  id <- fam_moverate$fam_g[j]
  mean_dispersal <- map_dbl(manysims, "mean_dispersal")
  se_dispersal <- map_dbl(manysims, "se_dispersal")
  seed_dispersion <- map_dbl(manysims, "seed_dispersion")

  out <- data.frame(dispersal = dispersal, id = id)
  fam_dispersal <- rbind.data.frame(fam_dispersal, out)

  out2 <- data.frame(id = id,
                     mean_dispersal = mean_dispersal,
                     se_dispersal = se_dispersal,
                     seed_dispersion = seed_dispersion)
  fam_dispersion <- rbind.data.frame(fam_dispersion, out2)

}

# Plots
plot_grid(
  fam_dispersal %>%
    ggplot(., aes(x = dispersal)) +
    geom_histogram(),
  plot_grid(fam_dispersion %>%
            ggplot(., aes(y = mean_dispersal)) +
            geom_boxplot(),
          fam_dispersion %>%
            ggplot(., aes(y = seed_dispersion)) +
            geom_boxplot()), nrow = 2)

save.image(file = "paper/simple_exp_runs.RData")

