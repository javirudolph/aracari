#################################################
## DATA GENERATION AND ANALYSIS - TOY MODEL CH3
################################################

# Libraries -----------------------------------------------
set.seed(20220201)

library(aracari)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(purrr)
library(fitdistrplus)

## Lnorm Param Functions ----------------------------------------------------
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
## Weights ----------------------------------------
# We assume four categories of individuals with increasing movement.
pis <- c( 0.1, 0.2, 0.3, 0.4)
# sum(pis)
# hist(pis)

## Components --------------------------------------
# Since movement lengths are only positive, we use a lognormal distribution to describe them

# CASE1

pars <- desired_mean_sd(mean = c(160, 300, 600, 1000), sd = c(90, 120, 160, 200))

# Special case I am considering where for each category the mean doubles, and then sd=mean*0.6
# pars <- desired_mean_sd(mean = c(150, 300, 600, 1200), sd = c(90, 180, 360, 720))

mus <- pars$mu
sigsqs <- pars$sigsq
dens_cols <- c("#264653", "#2a9d8f", "#f4a261", "#e76f51")

### Density Curves -------------------------------
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
categ <- rep(0,samp.size)

for(i in 1:samp.size){

  which.cat <- sample(1:4, size=1, replace=TRUE, prob=pis)
  all.samples[i] <- rlnorm(n=1,meanlog=mus[which.cat], sdlog=sigsqs[which.cat])
  categ[i] <- which.cat
}



### Histogram Samples ---------------------------------------------
# This is our f(x) from the main text
data.tail <- data.frame(values = all.samples, y = 100) %>% arrange(desc(values)) %>% filter(values >=250)
head(data.tail)

ggplot(data.frame(all.samples), aes(x = all.samples)) +
  geom_histogram(bins = 100) +
  geom_point(data = data.tail[1:50,], aes(x = values, y = y), color = "black", alpha = 0.5) +
  labs(y = "Frequency", x = "Distance") +
  theme_minimal() -> sampleshist


## Visualization -----------------------------------
plot_grid(densities_plot, sampleshist, nrow=2, labels="AUTO")
#ggsave("Ch3_samplesize/TestFig.png")

summary(all.samples)

# FITTING -------------------------------------
## Lomax functions --------------------------
### Old ------------------------------------

ft.nllike <- function(guess=init.betas, designmat=designmat,Y=Y){

  Xbeta         <- designmat%*%guess
  alphas        <- exp(Xbeta) # because alpha = ln(X*betas)
  n             <- length(Y)
  sumlogapy     <- sum(log(alphas+Y))
  k.hat         <- n/(sumlogapy - sum(Xbeta))
  lnft.yis     <- lomax.pdf(x=Y, alpha=alphas,k=k.hat,log.scale=TRUE)
  lnL           <- sum(lnft.yis)
  nll           <- -lnL

  if(is.infinite(nll)){nll <- .Machine$double.xmax}
  return(nll)
}

lomax.glm <- function(formula=~1, my.dataf, response){

  Y           <- response
  n           <- length(Y)
  designmat   <- model.matrix(formula, data=my.dataf)
  nbetas      <- ncol(designmat)
  init.betas  <- rep(1.5,nbetas)

  opt.out <- optim(par=init.betas, fn=ft.nllike, method="BFGS", control=list(trace=TRUE),
                   designmat=designmat, Y=Y)

  mles          <- opt.out$par
  nll.hat       <- opt.out$value
  BIC.mod       <- 2*nll.hat + nbetas*log(length(Y))
  Xbeta.hat     <- designmat%*%mles
  alphas.hat    <- exp(Xbeta.hat)
  sumlogapy.hat <- sum(log(alphas.hat+Y))
  k.hat         <- n/(sumlogapy.hat - sum(Xbeta.hat))

  out.list <- list(opt.out = opt.out, designmat=designmat,Y=Y, mles=mles, nll.hat=nll.hat, BIC.mod = BIC.mod,
                   alphas.hat=alphas.hat, k.hat=k.hat,data=my.dataf)

  return(out.list)

}

nllike.simp <- function(guess=c(1.5,1.5), Y=Y){

  parms         <- exp(guess) # because alpha = ln(X*betas)
  alpha         <- parms[1]
  k             <- parms[2]
  n             <- length(Y)
  #sumlogapy     <- sum(log(alphas+Y))
  #k.hat         <- n/(sumlogapy - sum(Xbeta))
  lnft.yis     <- lomax.pdf(x=Y, alpha=alpha,k=k,log.scale=TRUE)
  lnL           <- sum(lnft.yis)
  nll           <- -lnL
  #if(is.infinite(nll)){nll <- .Machine$double.xmax}
  return(nll)
}



### Current --------------------------------

## PDF ##
lomax.pdf <- function(x,alpha,k, log.scale=FALSE){

  if(log.scale==FALSE){out <- (k/(alpha+x))*(alpha/(alpha+x))^k
  }else{
    out <- log(k) + k*log(alpha) - (k+1)*log(alpha+x)
  }

  return(out)
}

## CDF ##
lomax.cdf <- function(x,alpha,k){

  return(1-(alpha/(alpha+x))^k)

}

# S(t) or P(X>=x)
lomax.st <- function(x=x, alpha=alpha, k=k){
  out <- alpha/(alpha+x)
  return(out^k)
}

## NLL ##

ft.nllike2 <- function(guess=init.betas, designmat=designmat,Y=Y){

  nbetasp1      <- length(init.betas)
  k             <- exp(guess[1])
  Xbeta         <- designmat%*%guess[2:nbetasp1]
  alphas        <- exp(Xbeta) # because alpha = ln(X*betas)
  n             <- length(Y)

  #sumlogapy     <- sum(log(alphas+Y))
  #k.hat         <- n/(sumlogapy - sum(Xbeta))
  lnft.yis     <- lomax.pdf(x=Y, alpha=alphas,k=k,log.scale=TRUE)
  lnL           <- sum(lnft.yis)
  nll           <- -lnL

  return(nll)
}

## NEW GLM ##

lomax.glm <- function(formula=~1, my.dataf, response){

  Y           <- response
  n           <- length(Y)
  designmat   <- model.matrix(formula, data=my.dataf)
  nbetas      <- ncol(designmat)
  init.betas  <- c(4,rep(5,nbetas))

  opt.out <- optim(par=init.betas, fn=ft.nllike2, method = "Nelder-Mead",
                   designmat=designmat, Y=Y)

  mles          <- opt.out$par
  nll.hat       <- opt.out$value
  BIC.mod       <- 2*nll.hat + nbetas*log(length(Y))
  Xbeta.hat     <- designmat%*%mles[-1]
  alphas.hat    <- exp(Xbeta.hat)
  #sumlogapy.hat <- sum(log(alphas.hat+Y))
  k.hat         <-  mles[1] #n/(sumlogapy.hat - sum(Xbeta.hat))

  out.list <- list(opt.out = opt.out, designmat=designmat,Y=Y, mles=mles, nll.hat=nll.hat, BIC.mod = BIC.mod,
                   alphas.hat=alphas.hat, k.hat=k.hat,data=my.dataf)

  return(out.list)

}

# TEST -------------------------------------------------------------

# Sample 1000 movements and fit a lomax

my.df <- data.frame(x = sample(all.samples, 1000))
mod1 <- lomax.glm(formula = ~1, my.dataf = my.df, response = my.df$x)

## ERROR --------------------------------------------------
# Can't find init.betas, but when you run the inside line by line it seems to work?

### Tried solution1 -----

formula = ~1
my.dataf = my.df
response = my.df$x


Y           <- response
n           <- length(Y)
designmat   <- model.matrix(formula, data=my.dataf)
nbetas      <- ncol(designmat)
init.betas  <- c(4,rep(5,nbetas))

opt.out <- optim(par=init.betas, fn=ft.nllike2, method = "Nelder-Mead",
                 designmat=designmat, Y=Y)

mles          <- opt.out$par
nll.hat       <- opt.out$value
BIC.mod       <- 2*nll.hat + nbetas*log(length(Y))
Xbeta.hat     <- designmat%*%mles[-1]
alphas.hat    <- exp(Xbeta.hat)
#sumlogapy.hat <- sum(log(alphas.hat+Y))
k.hat         <-  mles[1] #n/(sumlogapy.hat - sum(Xbeta.hat))

mod1 <- list(opt.out = opt.out, designmat=designmat,Y=Y, mles=mles, nll.hat=nll.hat, BIC.mod = BIC.mod,
                 alphas.hat=alphas.hat, k.hat=k.hat,data=my.dataf)


# What are the estimated parameters?
alpha <- mod1$alphas.hat[1] # equals exp(mod1$mles)
k <- mod1$k.hat


# QUESTION -----------------------------
# Am I doing this wrong? shouldn't the output of the CDF be a large number and the output of the S(t) be small?
# It's like they are backwards.


# Estimate the Pr(X <= x)
# Probability that a value drawn from this lomax is smaller or equal to 500
lomax.cdf(x = 500, alpha = alpha, k = k)
lomax.st(500, alpha, k)



# Calculate the mean using those parameters.

lomax.mean <- function(alpha=alpha, k=k){
  return(alpha/(k+1))
}

lomax.mean(alpha, k)
# This mean is too big, it doesn't make sense. I think we are missing an exponent or something somewhere.


# I deleted the rest of things, which correspond to lines 213 and on in the Ch3_datagen_toymodel.R script





