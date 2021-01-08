

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

# Calculate the mean location of seeds
x_m <- mean(seed_loc$xloc)
y_m <- mean(seed_loc$yloc)
mean_dispersal <- mean(seed_loc$dispersal)

# Calculate seed dispersion
xi <- (x_m - seed_loc$xloc)^2
yi <- (y_m - seed_loc$yloc)^2
seed_dispersion <- sum(sqrt(xi + yi))/nseeds

par(mfrow = c(1,3))
hist(movedist, main = "MD")
plot(x = xloc, y = yloc, type = "l", main = "Movement")
points(x = seed_loc$xloc, y = seed_loc$yloc, col = "blue", pch = 16)
hist(seed_loc$dispersal, main = "Dispersal")
abline(v = mean_dispersal, col = "red")


