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
                   alphas.hat=alphas.hat, k.hat=k.hat,data=my.dataf)

  return(out.list)

}
# Sample 1000 movements and fit a lomax

my.df <- data.frame(x = sample(all.samples, 1000))
mod1 <- lomax.glm(formula = ~1, my.dataf = my.df, response = my.df$x)

# What are the estimated parameters?
alpha <- mod1$alphas.hat[1] # equals exp(mod1$mles)
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

# Assess the fit using model quantiles to make a qqplot
lomax.quantile <- function(alpha, k, p){
  return(alpha*(((1-p)^(1/k))-1))
}
