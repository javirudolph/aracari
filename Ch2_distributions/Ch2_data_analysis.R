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

## CP bootstrapping -----------------------------------------

nboots <- 1000
nsamples <- length(ptpl$mpm)

cp.fits <- NULL

for(i in 1:nboots){

  boot.df <- sample(ptpl$mpm, size = nsamples, replace = TRUE)

  fit <- lapply(dist.used, function(x){fitdist(boot.df, distr = x)})

  out <- build_fits_df(fit) %>%
    mutate(param = param,
           dist = dist,
           kpars = kpars,
           dBIC = bic - min(bic),
           boot = paste0("boot_", i))

  cp.fits <- rbind.data.frame(cp.fits, out)
}

## Best fit model ---------------------------------------

cp.fits %>%
  dplyr::select(., bic, dist, dBIC, boot) %>%
  group_by(boot) %>%
  distinct() %>%
  ungroup() %>%
  filter(dBIC == 0) %>%
  count(dist) %>%
  mutate(prcnt_support = n/nboots) %>%
  mutate(model = "CP") %>%
  ggplot(., aes(x = model, y = prcnt_support, fill = dist)) +
  geom_col()


# No pooling fits ----------------------------------

ids <- unique(ptpl$id)

np.fits <- NULL


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

    np.fits <- rbind.data.frame(np.fits, out)

  }

}

## Best fit model ---------------------------------------

np.fits %>%
  dplyr::select(., bic, dist, dBIC, individual, boot) %>%
  group_by(individual, boot) %>%
  distinct() %>%
  group_by(individual) %>%
  filter(dBIC == 0) %>%
  count(dist) %>%
  mutate(prcnt_support = n/nboots) %>%
  ggplot(., aes(x = individual, y = prcnt_support, fill = dist)) +
  geom_col()


# Partial pooling fits ----------------------------------

fgs <- unique(ptpl$fam_g)

pp.fits <- NULL


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

    fit <- lapply(dist.used, function(x){fitdist(fg.df$mpm, distr = x)})

    out <- build_fits_df(fit) %>%
      mutate(param = param,
             dist = dist,
             kpars = kpars,
             dBIC = bic - min(bic),
             fgroup = fgs[j],
             boot = paste0("boot_", i))

    pp.fits <- rbind.data.frame(pp.fits, out)

  }

}

## Best fit model ---------------------------------------

pp.fits %>%
  dplyr::select(., bic, dist, dBIC, fgroup, boot) %>%
  group_by(fgroup, boot) %>%
  distinct() %>%
  group_by(fgroup) %>%
  filter(dBIC == 0) %>%
  count(dist) %>%
  mutate(prcnt_support = n/nboots) %>%
  ggplot(., aes(x = fgroup, y = prcnt_support, fill = dist)) +
  geom_col()

save(cp.fits, pp.fits, np.fits, file = "Ch2_distributions/fits_df.RData")

