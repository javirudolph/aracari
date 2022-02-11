# Correlation stuff

load("Ch2_distributions/Orig_data_KH/tidy_data.RData")

library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

raw_ptpl <- ptpl

ptpl %>%
  mutate(T_minutes = dt/60,
         Bird_ID = id,
         mpm = dist/T_minutes,
         R2n = lead(R2n)) %>%
  filter(mpm != 0) %>%
  group_by(Bird_ID) %>%
  add_tally() %>%
  filter(n >= 30) -> ptpl

S_L <- function(distance_vector){
  N <- length(distance_vector)
  d <- distance_vector

  d_sqrd <- vector()
  for(i in 1:(N-1)){
    ds <- d[i+1] - d[i]
    d_sqrd[i] <- ds^2
  }

  return(mean(d_sqrd))
}


ptpl %>%
  dplyr::select(burst, mpm) %>%
  group_by(Bird_ID, burst) %>%
  nest() %>%
  mutate(out = map_int(data, nrow),
         nperms = choose(out, out) * factorial(out),
         perm = ifelse(nperms > 999, 999, nperms)) -> a
summary(a$out)


ptpl %>%
  dplyr::select(burst, mpm) %>%
  group_by(Bird_ID, burst) %>%
  nest() %>%
  mutate(out = map_int(data, nrow),
         nperms = choose(out, out) * factorial(out),
         perm = ifelse(nperms > 999, 999, nperms)) %>%
  filter(., out > 2) %>%
  mutate(velvec = map(data, pull),
         ref_SL = map(velvec, S_L)) -> b


SL_permute <- function(d.vector, ntimes){

  out <- vector()
  for(i in 1:ntimes){
    perm <- sample(d.vector, size = length(d.vector), replace = FALSE)
    out[i] <- S_L(perm)
  }

  return(out)

}

b %>%
  mutate(perms = map(velvec, SL_permute, ntimes = perm)) -> c

ref_SLa <- b$ref_SL[[1]]
perm_SLa <- c$perms[[1]]


p_val_pos <- function(ref_SL, perm_SL){

  num <- sum(perm_SL <= ref_SL) + 1
  denom <- length(perm_SL) + 1

  return(num/denom)

}

p_val_neg <- function(ref_SL, perm_SL){

  num <- sum(perm_SL >= ref_SL) + 1
  denom <- length(perm_SL) + 1

  return(num/denom)

}


p_val_pos(ref_SLa, perm_SLa) # The null is that there is no correlation


# check positive autocorrelation

c %>%
  transmute(pval = map2(ref_SL, perms, p_val_pos)) -> e


e$perm <- c$perm
e$count <- c$out
e$pval <- as.numeric(e$pval)


e %>%
  mutate(sig = ifelse(pval <= 0.05, "yes", "no")) %>%
  ggplot(., aes(x = perm, y = pval, col = sig)) +
  geom_point()


# check negative autocorrelation

c %>%
  transmute(pval = map2(ref_SL, perms, p_val_neg)) -> f


f$perm <- c$perm
f$count <- c$out
f$pval <- as.numeric(f$pval)


f %>%
  mutate(sig = ifelse(pval <= 0.05, "yes", "no")) %>%
  ggplot(., aes(x = count, y = pval, col = sig)) +
  geom_point()

