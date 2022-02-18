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

  mu <- log(mu_x^2/(sqrt(mu_x^2+sd_x^2)))
  sigma <- log(1+(sd_x/mu_x^2))

  return(data.frame(mu, sigma))
}

# Now, for a Lnorm(mu, sigsq), get mean and var

lnorm_mean_var <- function(mu, sigma){

  lnorm_mean <- exp(mu + ((sigma^2)/2))
  lnorm_var  <- (exp(sigma^2)-1)*exp(2*mu+sigma^2)

  return(data.frame(lnorm_mean, lnorm_var))
}

#### Components --------------------------

# This is a function to make the multiple curves to plot
# I always struggle to remember how to make this

lnorm_densities_fx <- function(mu_x, sigsq_x, color_x){

  # I'm assuming we are giving variances, so need to change into sd.
  sd_x <- sqrt(sigsq_x)

  my_curvs <- purrr::map(1:length(mu_x), function(y) stat_function(fun = dlnorm,
                                                               args = list(meanlog = mu_x[y], sdlog = sd_x[y]),
                                                               color = color_x[y], size=1))
  return(my_curvs)
}
































