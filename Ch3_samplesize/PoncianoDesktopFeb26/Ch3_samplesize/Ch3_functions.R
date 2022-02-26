
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

# Lomax functions ---------------------------------
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

# Simple negative log-likelihood
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

# BIas and boxplots ----------------------------------------------------

bxplt_data_fx <- function(nreps){
  out <- data.frame()
  
  for(j in 1:nreps){
    samp_n_tests <- c(80, 200, 500, 800, 1000, 1600)
    mles_df <- data.frame(expand.grid(thresh_tests = thresh_tests, samp_n_tests = samp_n_tests))
    mles_df$nrep <- j
    
    for(i in 1:nrow(mles_df)){
      ith_n <- mles_df$samp_n_tests[i]
      ith_thresh <- mles_df$thresh_tests[i]
      ith_samps <- data.frame(x_star = sample(truth_df$x_samps, ith_n))
      ith_tail <- length(which(ith_samps$x_star >= ith_thresh))
      mles_df$tail_p[i] <- ith_tail/ith_n
      
      # Fit Lomax
      optim_out <- optim(par=log(c(1.5, 1.5)), fn=nllike.simp, method="BFGS", Y=ith_samps$x_star)
      mles_star <- exp(optim_out$par)
      alpha_star <- mles_star[1]
      k_star <- mles_star[2]
      mles_df$alpha_star[i] <-alpha_star
      mles_df$k_star[i] <- k_star
      
      # Check with Lomax GLM
      glm_out <- lomax.glm(formula=~1, my.dataf=ith_samps, response=ith_samps$x_star)
      alpha.2 <- glm_out$alphas.hat[1]
      k.2     <- glm_out$k.hat
      mles_df$alpha_glm[i] <-alpha.2
      mles_df$k_glm[i] <- k.2
      
      # Check with GP fit
      ith_evd <- fevd(ith_samps$x_star, threshold = 0, type = "GP")
      gp_scale <- summary(ith_evd)$par[1]
      gp_shape <- summary(ith_evd)$par[2]
      mles_df$alpha_GP[i] <- gp_scale*gp_shape
      mles_df$k_GP[i] <- 1/gp_shape
      mles_df$scale[i] <- gp_scale
      mles_df$shape[i] <- gp_shape
      
      
      
      # Estimate the tail
      #Lomax simple
      log_St_star <- lomax.St(x = ith_thresh,alpha = alpha_star,k = k_star,log.scale=TRUE)
      mles_df$log_St_hat[i] <- log_St_star
      
      # Lomax glm
      log.st.star2 <- lomax.St(x=ith_thresh,alpha=alpha.2,k=k.2,log.scale=TRUE)
      mles_df$log_St_hat2[i] <- log.st.star2
      
      # GP tail
      # mles_df$log_GP[i] <- log(pextRemes(ith_evd, ith_thresh, lower.tail = FALSE))
      mles_df$GP_theta_hat[i] <- pextRemes(ith_evd, ith_thresh, lower.tail = FALSE)
      
      
    }
    
    mles_df %>%
      right_join(., thetas[, c(1,3)]) -> mles_df
    out <- rbind.data.frame(out, mles_df)
  }
  
  return(out)
  
}

nonparam_boot <- function(B=5, star_data_vec=ith_samps$x_star, samp_size=ith_n, threshold_test=ith_thresh){
  
  bth_df <- data.frame(n_B = 1:B)
  ith_n <- samp_size
  ith_thresh <- threshold_test
  
  for(b in 1:B){
    # Sample with replacement from the given sample
    bth_samps <- data.frame(x_star = sample(star_data_vec, ith_n, replace = TRUE))
    
    # Estimate params using glm Lomax
    bth_glm <- lomax.glm(formula=~1, my.dataf=bth_samps, response=bth_samps$x_star)
    bth_alpha_hat <- bth_glm$alphas.hat[1]
    bth_k_hat    <- bth_glm$k.hat
    
    # Estimate params using Generalized Pareto
    bth_evd <- fevd(bth_samps$x_star, threshold = 0, type = "GP")
    bth_gp_scale <- summary(bth_evd)$par[1]
    bth_gp_shape <- summary(bth_evd)$par[2]
    
    
    # Estimate the tail
    # Lomax glm
    bth_log_St_hat <- lomax.St(x=ith_thresh,alpha=bth_alpha_hat,k=bth_k_hat,log.scale=TRUE)
    bth_df$log_St_star[b] <- bth_log_St_hat
    
    # GP tail
    bth_df$GP_theta_star[b] <- pextRemes(bth_evd, ith_thresh, lower.tail = FALSE)
    
  }
  return(bth_df)
  
}

param_boot_lomax <- function(B=5, star_data = ith_samps, star_data_vec=ith_samps$x_star, samp_size=ith_n, threshold_test=ith_thresh){
  
  bth_df <- data.frame(n_B = 1:B)
  ith_n <- samp_size
  ith_thresh <- threshold_test
  ith_glm <- lomax.glm(formula=~1, my.dataf=ith_samps, response=ith_samps$x_star)
  ith_alpha_hat <- ith_glm$alphas.hat[1]
  ith_k_hat    <- ith_glm$k.hat
  
  for(b in 1:B){
    # Sample with replacement from the given sample
    bth_samps <- data.frame(x_star = rlomax(n = ith_n, alpha = ith_alpha_hat, k = ith_k_hat))
    
    # Estimate params using glm Lomax
    bth_glm <- lomax.glm(formula=~1, my.dataf=bth_samps, response=bth_samps$x_star)
    bth_alpha_hat <- bth_glm$alphas.hat[1]
    bth_k_hat    <- bth_glm$k.hat
    
    # Estimate the tail
    # Lomax glm
    bth_log_St_hat <- lomax.St(x=ith_thresh,alpha=bth_alpha_hat,k=bth_k_hat,log.scale=TRUE)
    bth_df$log_St_star[b] <- bth_log_St_hat
    
  }
  return(bth_df)
  
}