

# Simulating under the lomax distribution ----

rlomax <- function(n=1000, alpha=2,k=4, plot.it=FALSE){

  hier.sims <- rep(0,n)
  for(i in 1:n){
  
    lam.rand <- rgamma(n=1, shape=k, rate=alpha)
    x        <- rexp(n=1,rate=lam.rand)
    hier.sims[i] <- x
  }

  if(plot.it==TRUE){
    range.sims <- range(hier.sims)
    ys <- seq(log(range.sims[1]), log(range.sims[2]), by=0.1)
    f.x <- function(x,alpha,k){((alpha/(x+alpha))^k)*(k/(x+alpha))}
    f.y <- f.x(x=exp(ys), alpha=alpha,k=k)*exp(ys)
    hist(log(hier.sims), freq=FALSE, xlab="log(x)", ylim=c(0,0.4),
         main= paste0(n," samples (in log scale) from the Lomax distribution"))
    points(ys,f.y, type="l", lwd=2, col="red")
  }
  return(hier.sims)
}  
  


# Lomax pdf, cdf, St=1-cdf------
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

lomax.St <- function(x,alpha,k, log.scale=FALSE){
  
  if(log.scale==FALSE){out <- (alpha/(alpha+x))^k
  }else{
    out <- k*log(alpha) - k*log(alpha+x)
  }
  
  return(out)
  
}

# Simple negative log-likelihood-----
nllike.simp <- function(guess=c(1.5,1.5), Y=Y){
  
  parms         <- exp(guess)
  alpha         <- parms[1]
  k             <- parms[2]
  n             <- length(Y)
  lnft.yis     <- lomax.pdf(x=Y, alpha=alpha,k=k,log.scale=TRUE)
  lnL           <- sum(lnft.yis)
  nll           <- -lnL
  return(nll)
}

# Lomax glm, which needs its own negative log-likelihood ftn----
ft.nllike2 <- function(guess, designmat=designmat,Y=Y){
  
  # For the glm idea to work we want ln(E[Lomax]) = ln(alpha/(k-1)) = XB
  # That is the proper link function.  Then, 
  # ln alpha = XB + ln(k-1), or
  # alpha = exp(XB + ln(k-1)) and it follows that
  # k = exp(ln(k-1)) +1.  To avoid problems with potential log(negative number) 
  # we optimize 'k-1', not 'k'
  
  nbetasp1      <- length(guess)
  lnkm1         <- guess[1]
  Xbeta         <- designmat%*%guess[2:nbetasp1]
  ln.alphas     <- Xbeta + lnkm1
  alphas        <- exp(ln.alphas)
  k             <- exp(lnkm1)+1
  n             <- length(Y)
  
  #sumlogapy     <- sum(log(alphas+Y))
  #k.hat         <- n/(sumlogapy - sum(Xbeta))
  lnft.yis     <- lomax.pdf(x=Y, alpha=alphas,k=k,log.scale=TRUE)
  lnL           <- sum(lnft.yis)
  nll           <- -lnL
  
  return(nll)
}

# NEW GLM-----
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
  BIC.mod       <- 2*nll.hat + (nbetas+1)*log(length(Y))
  Xbeta.hat     <- designmat%*%mles[-1]
  lnkm1.hat     <- mles[1]
  k.hat         <- exp(lnkm1.hat)+1
  alphas.hat    <- exp(Xbeta.hat +lnkm1.hat)

  out.list <- list(opt.out = opt.out, designmat=designmat,Y=Y, mles=mles, nll.hat=nll.hat, BIC.mod = BIC.mod,
                   alphas.hat=alphas.hat, k.hat=k.hat,data=my.dataf)
  
  return(out.list)
  
}




# Procedural section:-----
### Idea: simulate data and estimate the model parameters as well as
### the P(X>threshold) using these parameter estimates, for a fixed
### sample size.  I will do this comparing the simple estimation ftn.
### to the glm estimation function, with the same simulations


true.alpha <- 2
true.k <- 4
test.thres <- 100
samp.size <- 500
log.true.St <- lomax.St(x=test.thres,alpha=true.alpha, k=true.k,log.scale=TRUE)


B <- 2000

mles.mat <- matrix(0,nrow=B,ncol=2)
colnames(mles.mat) <- c("alpha.hat", "k.hat")
log.St.vec.hat <- rep(0,B)
sims.mat <- matrix(0, nrow=samp.size,ncol=B)

mles.mat2 <- mles.mat
log.St.vec.hat2 <- log.St.vec.hat

log.guess <- log(c(1.5,1.5))

for(i in 1:B){
  
  # Method of estimation 1
  x.star <- rlomax(n=samp.size, alpha=true.alpha, k=true.k)
  opt.out <- optim(par=log.guess, fn=nllike.simp, method="BFGS", Y=x.star)
  mles.star <- exp(opt.out$par)
  mles.mat[i,] <- mles.star 
  alpha.star <- mles.star[1]
  k.star <- mles.star[2]
  log.St.star <- lomax.St(x=test.thres,alpha=alpha.star,k=k.star,log.scale=TRUE)
  log.St.vec.hat[i] <- log.St.star 
  
  # Method of estimation 2
  df.star <- data.frame(Y=x.star)
  glm.out <- lomax.glm(formula=~1, my.dataf=df.star, response=df.star$Y)
  alpha.2 <- glm.out$alphas.hat[1]
  k.2     <- glm.out$k.hat
  log.st.star2 <- lomax.St(x=test.thres,alpha=alpha.2,k=k.2,log.scale=TRUE)
  
  mles.mat2[i,] <- c(alpha.2,k.2)
  log.St.vec.hat2[i] <- log.st.star2
  

}


par(mfrow=c(2,3), oma=c(2,3,2,1),mar=c(1,2.5,1,1),bty="n")

boxplot(mles.mat[,1]/true.alpha, outline=FALSE,main=expression(alpha),ylab="")
abline(h=1,lty=2)
boxplot(mles.mat[,2]/true.k, outline=FALSE, main=expression(k),ylab="")
abline(h=1,lty=2)
boxplot(exp(log.St.vec.hat-log.true.St), outline=FALSE, main=expression(S[t]),ylab="")
abline(h=1,lty=2)

boxplot(mles.mat2[,1]/true.alpha, outline=FALSE,main=expression(alpha),ylab="")
abline(h=1,lty=2)
boxplot(mles.mat2[,2]/true.k, outline=FALSE, main=expression(k),ylab="")
abline(h=1,lty=2)
boxplot(exp(log.St.vec.hat2-log.true.St), outline=FALSE, main=expression(S[t]),ylab="")
abline(h=1,lty=2)

mtext(text="Relative bias = Estimated/True", side=2, outer=TRUE)



# Profile likelihood function -----
prof.nllike <- function(lalpha.guess=1.5,k=4, Y=Y){
  
  alpha         <- exp(lalpha.guess)
  n             <- length(Y)
  lnft.yis     <- lomax.pdf(x=Y, alpha=alpha,k=k,log.scale=TRUE)
  lnL           <- sum(lnft.yis)
  nll           <- -lnL
  return(nll)
}

true.alpha <- 2
true.k <- 4
samp.size <- 500
trial.sim <-  rlomax(n=samp.size, alpha=true.alpha, k=true.k)
joint.opt <- optim(par=log.guess, fn=nllike.simp, method="BFGS", Y=trial.sim)
khat.joint <- exp(joint.opt$par[2])


k.vec <- seq(1.01,100,by=0.01)
nks    <- length(k.vec)
prof.loglike <- rep(0,nks)
prof.alphahats <- rep(0,nks)
for(i in 1:nks){
  
  ith.k <- k.vec[i]
  ith.opt <- optim(par=log(4), fn=prof.nllike, method="BFGS", k=ith.k, Y=trial.sim)
  prof.loglike[i] <- ith.opt$value
  prof.alphahats[i] <- exp(ith.opt$par)

}

plot(k.vec, -prof.loglike, type="l", col="red", lwd=2)

prof.alpha.hat <- k.vec[prof.loglike==min(prof.loglike)]











