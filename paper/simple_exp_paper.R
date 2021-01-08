

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


data(ptpl)

average_movrate <- mean(ptpl$mpm)

indiv_moverate <- ptpl %>%
  group_by(Bird_ID) %>%
  summarise(movrate = mean(mpm))

fam_moverate <- ptpl %>%
  group_by(fam_g) %>%
  summarise(movrate = mean(mpm))

ids <- ptpl %>% distinct(., Bird_ID, fam_g)

manysims <- vector("list", 10)
for(i in 1:length(manysims)){
  manysims[[i]] <- sim_run(28)
  #manysims[[i]]$sim_run <- paste("sim_", i)
}


dispersal <- map_df(manysims, "seed")$dispersal
hist(dispersal)
