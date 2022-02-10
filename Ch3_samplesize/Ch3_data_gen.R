##########################################
## DATA GENERATION CH3
#########################################

# I think the overall question for this is chapter is:
#         Can we accurately estimate the true proportion of rare events produced by a complex model using a simpler hierarchical model?

# Our complex model is one that considers individual variation in animal movement, with increasing mean and variance for specific groups
# We use a mixture modeling approach, where each component distribution is associated to an individual's state or age.
# We assume that an individual's stage will affect their movement.
# We assume, for example, very young birds might not move as far, as will nesting or molting birds, whereas adults will move more.

# The overall assumption for this model is that individuals vary in how much they move, which we determine based on step lengths only.
# With some groups having an overall larger variance and mean than other groups.
# The proportions, or weights, of the mixture, should follow a bell curve for now
# With the average movers being the majority of the group, and low or high movement individuals being less abundant.

# Libraries -----------------------------------------------
set.seed(20220201)

library(aracari)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(purrr)
library(fitdistrplus)

## Functions ----------------------------------------------------
# Building this function so that we get the parameters for a desired mean and standard deviation.
desired_mean_sd <- function(mean, sd){

  mu <- log(mean*2/(sqrt(mean*2+sd*2)))
  sigsq <- log(1+(sd*2/mean*2))

  return(data.frame(mu = mu, sigsq = sigsq))
}

exp.val <- function(mu, sigsq){
  exp(mu + sqrt(sigsq))
}


# MIXTURE MODEL --------------------------------------------------
## Mixture proportions ----------------------------------------
# We assume four categories of individuals with increasing movement.
pis <- c( 0.1, 0.2, 0.3, 0.4)
sum(pis)
hist(pis)

## Mixture components --------------------------------------
# Since movement lengths are only positive, we use a lognormal distribution to describe them

# CASE1: same variance for all, but change expected values

pars <- desired_mean_sd(mean = c(160, 300, 600, 1000), sd = c(90, 120, 160, 200))

# Special case I am considering where for each category the mean doubles, and then sd=mean*0.6
# pars <- desired_mean_sd(mean = c(150, 300, 600, 1200), sd = c(90, 180, 360, 720))

mus <- pars$mu
sigsqs <- pars$sigsq
dens_cols <- c("#264653", "#2a9d8f", "#f4a261", "#e76f51")

### PLOT the mixture components -------------------------------
lnorm_densities <- purrr::map(1:4, function(y) stat_function(fun = dlnorm,
                                                             args = list(meanlog = mus[y], sdlog = sigsqs[y]),
                              color = dens_cols[y], size=1))
ggplot() +
  lnorm_densities +
  theme_minimal() +
  labs(y = "Density") +
  lims(x = c(0, 150)) -> densities_plot


# SAMPLING --------------------------------------------------------------

samp.size <- 50000
all.samples <- rep(0,samp.size)
for(i in 1:samp.size){

  which.cat <- sample(1:4, size=1, replace=TRUE, prob=pis)
  all.samples[i] <- rlnorm(n=1,meanlog=mus[which.cat], sdlog=sigsqs[which.cat])
}

## PLOT the samples ---------------------------------------------
# This is our f(x) from the main text
data.tail <- data.frame(values = all.samples, y = 100) %>% arrange(desc(values)) %>% filter(values >=250)
head(data.tail)

ggplot(data.frame(all.samples), aes(x = all.samples)) +
  geom_histogram(bins = 100) +
  geom_point(data = data.tail[1:50,], aes(x = values, y = y), color = "black", alpha = 0.5) +
  labs(y = "Frequency", x = "Distance") +
  theme_minimal() -> sampleshist


# VIZ -----------------------------------
plot_grid(densities_plot, sampleshist, nrow=2, labels="AUTO")
ggsave("Ch3_samplesize/Fig1.png")

summary(all.samples)

# FITTING -------------------------------------
# We use a Lomax distribution, which is actually a special case of generalized Pareto distribution
# The Lomax or Patro type II has support on x=0, and is a heavy tail distribution
# Jose's document shows that it arises from incorporating variability in the rate parameter of an exponential distribution, using a gamma distribution - which he calls the hierarchical
# It only has two parameters, alpha(shape) and lambda(scale)

### Lomax functions --------------------------

lomax.pdf <- function(x,alpha,k, log.scale=FALSE){

  if(log.scale==FALSE){out <- (k/(alpha+x))*(alpha/(alpha+x))^k
  }else{
    out <- log(k) + k*log(alpha) - (k+1)*log(alpha+x)
  }

  return(out)
}


lomax.cdf <- function(x,alpha,k){

  return(1-(alpha/(alpha+x))^k)

}

ft.nllike <- function(guess=init.betas, designmat=designmat,Y=Y){

  Xbeta         <- designmat%*%guess
  alphas        <- exp(Xbeta) # because alpha = ln(X*betas)
  n             <- length(Y)
  sumlogapy     <- sum(log(alphas+Y))
  k.hat         <- n/(sumlogapy - sum(Xbeta))
  lnft.yis     <- lomax.pdf(x=Y, alpha=alphas,k=k.hat,log.scale=TRUE)
  lnL           <- sum(lnft.yis)
  nll           <- -lnL
  return(nll)
}

lomax.glm <- function(formula=~1, my.dataf, response){

  Y           <- response
  n           <- length(Y)
  designmat   <- model.matrix(formula, data=my.dataf)
  nbetas      <- ncol(designmat)
  init.betas  <- rep(1.5,nbetas)

  opt.out <- optim(par=init.betas, fn=ft.nllike, method="BFGS",
                   designmat=designmat, Y=Y)

  mles          <- opt.out$par
  nll.hat       <- opt.out$value
  BIC.mod       <- 2*nll.hat + nbetas*log(length(Y))
  Xbeta.hat     <- designmat%*%mles
  alphas.hat    <- exp(Xbeta.hat)
  sumlogapy.hat <- sum(log(alphas.hat+Y))
  k.hat         <- n/(sumlogapy.hat - sum(Xbeta.hat))

  out.list <- list(opt.out = opt.out, designmat=designmat,Y=Y, mles=mles, nll.hat=nll.hat, BIC.mod = BIC.mod,
                   alphas.hat=alphas.hat, k.hat=k.hat,data=data)

  return(out.list)

}

# Sample 1000 movements and fit a lomax

my.df <- data.frame(x = sample(all.samples, 1000))
mod1 <- lomax.glm(formula = ~1, my.dataf = my.df, response = my.df$x)

# What are the estimated parameters?
alpha <- mod1$alphas.hat[1]
k <- mod1$k.hat

# Estimate the Pr(X <= x)
# Probability that a value drawn from this lomax is smaller or equal to 500
lomax.cdf(x = 500, alpha = alpha, k = k)

# Probability that it is greater.
lomax.st <- function(x=x, alpha=alpha, k=k){
  out <- alpha/(alpha+x)
  return(out^k)
}

lomax.st(500, alpha, k)



# calculate the mean of the fit

lomax.mean <- function(alpha=alpha, k=k){
  return(alpha/(k+1))
}

lomax.mean(alpha, k)

# Assess the fit using model quantiles
lomax.quantile <- function(alpha, k, p){
  return(alpha*(((1-p)^(1/k))-1))
}





### Try to fit lomax using a package ------------------------------
library(extRemes)

# The lomax is a special case of the Generalized Pareto

gp.fit <- fevd(my.df$x, type = "GP", threshold = 50)

# Poor fit to the tail
plot(gp.fit)

location <- 0

# GP shape is 1/alpha in lomax

scale <- gp.fit$results$par[1]
shape <- gp.fit$results$par[2]

alpha.gp <- as.numeric(1/shape)
k.gp <- as.numeric(scale*alpha)

# Compare the parameters estimated:

c(alpha, alpha.gp)
c(k, k.gp)








