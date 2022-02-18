##########################################
## Javiera Rudolph                      ##
### Functions used for Ch3              ##
## February 18, 2022                    ##
##########################################



##########################################
# Data simulation ------------------------


## Mixture model -------------------------

# The weights right now are just set as values
# but this can change, and we could use a distribution
# to sample them from

### Weights ------------------------------

# This is not created yet, but on the list

generate.pis <- function(n){
  n
}


### Components --------------------------

#### Parameters -------------------------
# since I get confused with logs and such
# The following function helps gives the
# parameters for a lognormal that has the
# desired mean and variance:

desired_mean_sd <- function(mu_x, sd_x){

  sigsq <- sd_x^2

  mu <- log(mu_x^2/(sqrt(mu_x^2+sigsq)))
  sigma_sq <- log(1+(sigsq/mu_x^2))
  sigma <- sqrt(sigma_sq)

  return(data.frame(meanlog = mu, sigma_sq = sigma_sq, sdlog = sigma))
}

# Now, for a Lnorm(meanlog, sdlog), get mean and var

lnorm_mean_var <- function(mean_log, sd_log){

  lnorm_mean <- exp(mean_log + ((sd_log^2)/2))
  lnorm_var  <- (exp(sd_log^2)-1)*exp(2*mean_log+sd_log^2)
  lnorm_sd   <- sqrt(lnorm_var)

  return(data.frame(lnorm_mean, lnorm_var, lnorm_sd))
}

#### Components --------------------------

# This is a function to make the multiple curves to plot
# I always struggle to remember how to make this

lnorm_densities_fx <- function(meanlog_x, sdlog_x, color_x){

  my_curvs <- purrr::map(1:length(meanlog_x), function(y) stat_function(fun = dlnorm,
                                                               args = list(meanlog = meanlog_x[y], sdlog = sdlog_x[y]),
                                                               color = color_x[y], size=1))
  return(my_curvs)
}




# LOMAX -------------------------------------------
# All the lomax functions

# FITTING -------------------------------------
# We use a Lomax distribution, which is actually a special case of generalized Pareto distribution
# The Lomax or Patro type II has support on x=0, and is a heavy tail distribution
# Jose's document shows that it arises from incorporating variability in the rate parameter of an exponential distribution, using a gamma distribution - which he calls the hierarchical
# It only has two parameters, alpha(shape) and lambda(scale)

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

lomax.st <- function(x=x, alpha=alpha, k=k){
  out <- alpha/(alpha+x)
  return(out^k)
}

lomax.mean <- function(alpha=alpha, k=k){
  return(alpha/(k-1))
}


ft.nllike <- function(guess=init.betas, designmat=designmat,Y=Y){

  nbetasp1      <- length(guess)
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


lomax.glm <- function(formula=~1, my.dataf, response){

  Y           <- response
  n           <- length(Y)
  designmat   <- model.matrix(formula, data=my.dataf)
  nbetas      <- ncol(designmat)
  init.betas  <- c(4,rep(5,nbetas))

  opt.out <- optim(par=init.betas, fn=ft.nllike, method = "Nelder-Mead",
                   designmat=designmat, Y=Y)

  mles          <- opt.out$par
  nll.hat       <- opt.out$value
  BIC.mod       <- 2*nll.hat + nbetas*log(length(Y))
  Xbeta.hat     <- designmat%*%mles[-1]
  alphas.hat    <- log(Xbeta.hat)
  #sumlogapy.hat <- sum(log(alphas.hat+Y))
  k.hat         <-  mles[1] #n/(sumlogapy.hat - sum(Xbeta.hat))

  out.list <- list(opt.out = opt.out, designmat=designmat,Y=Y, mles=mles, nll.hat=nll.hat, BIC.mod = BIC.mod,
                   alphas.hat=alphas.hat, k.hat=k.hat,data=my.dataf)

  return(out.list)

}




























