mus <- c(3,5,7,11,20,30)
sigsqs <- c(1,1.5,1.5,2,2,3)
xs <- seq(0,35,by=0.1)
# Pyramid age structure
pis <- c(0.4, 0.2,0.20,0.10,0.07,0.03)

samp.size <- 50000
all.samples <- rep(0,samp.size)
for(i in 1:samp.size){

  which.age <- sample(1:6, size=1, replace=TRUE, prob=pis)
  all.samples[i] <- rnorm(n=1,mean=mus[which.age], sd=sigsqs[which.age])
}


par(mfrow=c(2,1))
plot(xs,dnorm(x=xs, mean = mus[1], sd=1), type ="l", col="gray",xlab="Trait distribution per age class", ylab="Density")
points(xs,dnorm(x=xs, mean = mus[2], sd=sigsqs[2]), type ="l", col="gray")
points(xs,dnorm(x=xs, mean = mus[3], sd=sigsqs[3]), type ="l", col="gray")
points(xs,dnorm(x=xs, mean = mus[4], sd=sigsqs[4]), type ="l", col="gray")
points(xs,dnorm(x=xs, mean = mus[5], sd=sigsqs[5]), type ="l", col="gray")
points(xs,dnorm(x=xs, mean = mus[6], sd=sigsqs[6]), type ="l", col="gray")

hist(all.samples, main= paste0("Sample from proportions 0.4, 0.2,0.20,0.15,0.05"))

# JAVI - using a lognormal

meanlogs <- c(1, 1.5, 2, 2.5)
sdlogs <- c(0.5, 0.3, 1, 0.75)
xs <- seq(0,10, by = 0.1)


plot(xs, dlnorm(x = xs, meanlog = meanlogs[1], sdlog = sdlogs[1]), type = "l", col = "grey", ylab = "Density")
points(xs, dlnorm(x = xs, meanlog = meanlogs[2], sdlog = sdlogs[2]), type = "l", col = "red")
points(xs, dlnorm(x = xs, meanlog = meanlogs[3], sdlog = sdlogs[3]), type = "l", col = "blue")
points(xs, dlnorm(x = xs, meanlog = meanlogs[4], sdlog = sdlogs[4]), type = "l", col = "green")

