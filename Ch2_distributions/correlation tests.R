# Correlation stuff

load("Ch2_distributions/Orig_data_KH/tidy_data.RData")

library(dplyr)
library(purrr)
library(tidyr)

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
         nperms = choose(out, out) * factorial(out)) %>%
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
  mutate(perms = map(velvec, SL_permute, ntimes = 6)) -> c


