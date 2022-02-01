##########################################
## DATA GENERATION CH3
#########################################

# I think the overall question for this is chapter is:
#         Can we accurately estimate the true proportion of rare events produced by a complex model using a simpler hierarchical model?

# Our complex model is one that considers individual variation in animal movement, with increasing mean and variance for specific groups
# We use a mixture modeling approach, where each component distribution is associated to an individual's state or age

# The overall assumption for this model is that individuals vary in how much they move, which we determine based on step lengths only.
# With some groups having an overall larger variance and mean than other groups.
# The proportions, or weights, of the mixture, are only assumed for now to follow a bell-shaped curve.

# Libraries -----------------------------------------------
set.seed(20220201)

library(aracari)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(fitdistrplus)
library(Hmisc)

###########################
# Mixture components
###########################
# Since we will be modeling step lengths, we will use a lognormal distribution, which is limited to only positive numbers.
# We get our parameters inspiration from the rediscretized step lengths to 15 minutes from the ptpl dataset
load("Ch2_distributions/Orig_data_KH/tidy_data.RData")

# Rediscretizing the data
# Big assumptions here, we are ignoring potential correlation between consecutive relocations
# Assuming each observation is independent
ptpl %>% drop_na(rel.angle) %>% mutate(stp.len = (dist/dt)*900) %>%
  dplyr::select(id, sgroup, stp.len) -> lengths_df

# Looking at the mean and sd at the group or indidivual level

lengths_df %>%
  as_tibble() %>%
  group_by(sgroup) %>%
  summarise(across(c(stp.len), list(mean = mean, sd = sd))) %>%
  arrange(stp.len_mean) %>%
  summary()



# Building this function so that we get the parameters for a desired mean and standard deviation.
get_lnorm_params <- function(mean, sd){
  mu <- log(mean*2/(sqrt(mean*2+sd*2)))
  sigmasqrd <- log(1+(sd*/mean*2))

  return(c(mu, sigmasqrd))
}


# We will consider a range from 100-600 meters for both mean and sd to use in the lognormal distribution.
# And build 4 different groups


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

meanlogs <- c(3.011376, 2.670338, 2.842328,
              #2.830963,
              1.637819
              #,2.351636, 2.117476
              )
sdlogs <- c(1.343584, 1.114171, 1.032824,
            #1.010781,
            1.487352
            #, 1.139045, 1.150798,
            )

samepis <- 1/7
pis <- rep(samepis, 7)


samp.size <- 50000
all.samples <- rep(0,samp.size)
for(i in 1:samp.size){

  which.group <- sample(1:7, size=1, replace=TRUE, prob=pis)
  all.samples[i] <- rlnorm(n=1,mean=meanlogs[which.group], sd=sdlogs[which.group])
}

xs <- seq(0,20, by = 0.1)

par(mfrow=c(2,1))
plot(xs, dlnorm(x = xs, meanlog = meanlogs[1], sdlog = sdlogs[1]), type = "l", col = "grey", ylab = "Density")
points(xs, dlnorm(x = xs, meanlog = meanlogs[2], sdlog = sdlogs[2]), type = "l", col = "grey")
points(xs, dlnorm(x = xs, meanlog = meanlogs[3], sdlog = sdlogs[3]), type = "l", col = "grey")
points(xs, dlnorm(x = xs, meanlog = meanlogs[4], sdlog = sdlogs[4]), type = "l", col = "grey")
points(xs, dlnorm(x = xs, meanlog = meanlogs[5], sdlog = sdlogs[5]), type = "l", col = "grey")
points(xs, dlnorm(x = xs, meanlog = meanlogs[6], sdlog = sdlogs[6]), type = "l", col = "grey")
points(xs, dlnorm(x = xs, meanlog = meanlogs[7], sdlog = sdlogs[7]), type = "l", col = "grey")

hist(all.samples, main= paste0("Sample"))
