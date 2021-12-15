# 20211207 Data analysis

# Libraries -----------------------------------------------
set.seed(20211207)

library(aracari)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(fitdistrplus)

# Data ----------------------------------------------------
# The script for this data is in the data-raw folder
load("data/ptpl.rda")

ptpl %>%
  filter(., mpm != 0) -> ptpl


dist.used <- c("gamma", "weibull", "lnorm")
param <- c("shape", "rate", "shape", "scale", "meanlog", "sdlog")
dist <- c("gamma", "gamma", "weibull", "weibull", "lnorm", "lnorm")
kpars <- c(2, 2, 2, 2, 2, 2)

# Complete pooling fits -------------------------------------

## Function to save the fit information and parameters -----

build_fits_df <- function(x){

  nm <- deparse(substitute(x))

  x %>%
    map(., `[`, c("estimate", "sd", "loglik", "n", "bic")) %>%
    map(as_tibble) %>%
    bind_rows()
}

# Regular fits -----------------------------------------

## CP fit ------------------------------------------------

cp.fit <- lapply(dist.used, function(x){fitdist(ptpl$mpm, distr = x)})

reg_fits_info <- build_fits_df(cp.fit) %>%
  mutate(param = param,
         dist = dist,
         kpars = kpars,
         dBIC = bic - min(bic),
         ID = "CP",
         model = "CP")

## NP fit ---------------------------------------------

ids <- unique(ptpl$id)

for(j in 1:length(ids)){

  indv.df <- ptpl %>%
    filter(., id == ids[j])

  fit <- lapply(dist.used, function(x){fitdist(indv.df$mpm, distr = x)})

  out <- build_fits_df(fit) %>%
    mutate(param = param,
           dist = dist,
           kpars = kpars,
           dBIC = bic - min(bic),
           ID = ids[j],
           model = "NP")

  reg_fits_info <- rbind.data.frame(reg_fits_info, out)

}

## PP fit -----------------------------------------------

fgs <- unique(ptpl$fam_g)

for(j in 1:length(fgs)){

  fam.df <- ptpl %>%
    filter(., fam_g == fgs[j])

  fit <- lapply(dist.used, function(x){fitdist(fam.df$mpm, distr = x)})

  out <- build_fits_df(fit) %>%
    mutate(param = param,
           dist = dist,
           kpars = kpars,
           dBIC = bic - min(bic),
           ID = fgs[j],
           model = "PP")

  reg_fits_info <- rbind.data.frame(reg_fits_info, out)

}

#################################################################

###################################################################


# Bootstrapping --------------------------------------------------
## CP bootstrapping -----------------------------------------

nboots <- 1000
nsamples <- length(ptpl$mpm)

cp.boot.fits <- NULL

for(i in 1:nboots){

  boot.df <- sample(ptpl$mpm, size = nsamples, replace = TRUE)

  fit <- lapply(dist.used, function(x){fitdist(boot.df, distr = x)})

  out <- build_fits_df(fit) %>%
    mutate(param = param,
           dist = dist,
           kpars = kpars,
           dBIC = bic - min(bic),
           boot = paste0("boot_", i))

  cp.boot.fits <- rbind.data.frame(cp.boot.fits, out)
}


# No pooling fits ----------------------------------

ids <- unique(ptpl$id)

np.boot.fits <- NULL


for(i in 1:nboots){

  boot.df <- ptpl %>%
    ungroup() %>%
    dplyr::select(.,id, mpm) %>%
    group_by(id) %>%
    nest() %>%
    ungroup() %>%
    mutate(n = map_int(data, nrow)) %>%
    mutate(samp = map2(data, n, sample_n, replace = TRUE)) %>%
    dplyr::select(-data) %>%
    unnest(samp)

  for(j in 1:length(ids)){

    indv.df <- boot.df %>%
      filter(., id == ids[j])

    fit <- lapply(dist.used, function(x){fitdist(indv.df$mpm, distr = x)})

    out <- build_fits_df(fit) %>%
      mutate(param = param,
             dist = dist,
             kpars = kpars,
             dBIC = bic - min(bic),
             individual = ids[j],
             boot = paste0("boot_", i))

    np.boot.fits <- rbind.data.frame(np.boot.fits, out)

  }

}



# Partial pooling fits ----------------------------------

fgs <- unique(ptpl$fam_g)

pp.boot.fits <- NULL


for(i in 1:nboots){

  boot.df <- ptpl %>%
    ungroup() %>%
    dplyr::select(.,fam_g, mpm) %>%
    group_by(fam_g) %>%
    nest() %>%
    ungroup() %>%
    mutate(n = map_int(data, nrow)) %>%
    mutate(samp = map2(data, n, sample_n, replace = TRUE)) %>%
    dplyr::select(-data) %>%
    unnest(samp)

  for(j in 1:length(fgs)){

    fg.df <- boot.df %>%
      filter(., fam_g == fgs[j])

    fit <- lapply(dist.used, function(x){fitdist(fg.df$mpm, distr = x)})

    out <- build_fits_df(fit) %>%
      mutate(param = param,
             dist = dist,
             kpars = kpars,
             dBIC = bic - min(bic),
             fgroup = fgs[j],
             boot = paste0("boot_", i))

    pp.boot.fits <- rbind.data.frame(pp.boot.fits, out)

  }

}



save.image(file = "Ch2_distributions/fits_df.RData")
save(nboots, dist, dist.used, ids, fgs, kpars, param, reg_fits_info,
     cp.boot.fits, pp.boot.fits, np.boot.fits,
     file = "Ch2_distributions/Ch2_fits.RData")
