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
































