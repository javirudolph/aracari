# Script to run simulations using parameters from aracari movement rates
# LIBRARIES --------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(fitdistrplus)
library(Hmisc)

set.seed(927)

# Functions --------------------------------------------------------------------------
sim_movement <- function(prm, t = 1000, plot.it = TRUE, return.data.frame = FALSE){
  tru.rate <- round(prm, 3)
  movedist <- rexp(t, rate = 1/tru.rate)
  angle <- runif(t, min = 0, max = 360)
  distx <- movedist*cos(angle)
  xloc <- c(0, cumsum(distx))
  disty <- movedist*sin(angle)
  yloc <- c(0, cumsum(disty))
  animalTraj <- data.frame(time = 0:t, xloc = xloc, yloc = yloc)
  if(plot.it == TRUE){
    plot(x = xloc, y = yloc, type = "l", main = paste("Rate=",tru.rate))
  }
  if(return.data.frame == TRUE){
    return(animalTraj)
  }
}

sim_seeds <- function(nseeds = 20, m.prms = NULL,...){
  grt <- round(rgamma(nseeds, shape = 4, scale = 5))
  t.grt <- max(grt)

  df <- sim_movement(m.prms, t = t.grt, plot.it = FALSE, return.data.frame = TRUE)

  df %>%
    left_join(., data.frame(s.id = 1:nseeds, time = grt), by = "time") %>%
    mutate(disp = sqrt(xloc^2+yloc^2)) -> df

  return(df)
}

summ_seeds <- function(df = NULL){
  df %>%
    drop_na(s.id) %>%
    mutate(xi = (mean(xloc)-xloc)^2,
           yi = (mean(yloc)-yloc)^2) %>%
    summarise(x = mean(xloc),
              y = mean(yloc),
              av.disp = mean(disp),
              se.disp = sd(disp)/sqrt(n()),
              dsprsn = sum(sqrt(xi+yi))/n()) -> s.df
  return(s.df)
}


# LOAD the data ------------------------------------------------------------------------

load("Ch1_movement_rates/ptpl.RData")

null_moverate <- data.frame(Bird_ID = unique(ptpl$Bird_ID), movrate = mean(ptpl$mpm))

indiv_moverate <- ptpl %>%
  group_by(Bird_ID, fam_g) %>%
  summarise(movrate = mean(mpm))

fam_moverate <- ptpl %>%
  group_by(fam_g) %>%
  summarise(movrate = mean(mpm),
            sd = sd(mpm))

ids <- ptpl %>% distinct(., Bird_ID, fam_g)

# Color palette ----------------------------------------
my.cols1 <- c("#23262f","#717492","#b3a82a","#c94f21","#980012","#0d907a","#b9bec3")


logfit <- fitdist(indiv_moverate$movrate, distr = 'lnorm')

movrate_cp <- mean(indiv_moverate$movrate)
movrate_cp_ln <- as.numeric(exp(logfit$estimate[1]))
# Which estimate to use?
# So, we are saying here that complete pooling takes the average of all the data as one. We took all the data, then fitted a lognormal distribution. So, the parameter we use, is the average movement rate that the lognormal gives us. Since we are assuming that the lognormal distribution will describe this population's movement rates across individuals. So, we use the expected value (mean) of the fit to describe the average movement rate of the population.

movrate_cp_ln

# How many seeds to use? Landon used 100 because he was taking averages.
# I guess we are focusing on a very short time frame: what does the bird do after it eats those seeds at one tree, and how does it move until it drops them all? So, only use 5 seeds.

nseeds <- 5

## CP Simulations -----------------------------------------------------------------------
# Generate seed dispersal data under a complete pooling scenario

kruns <- 10000
nseeds <- 5

cp.df <- NULL
cp.summ.df <- NULL

for(k in 1:kruns){
  a <- sim_seeds(m.prms = movrate_cp_ln, nseeds = nseeds) %>%
    mutate(run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
           model = "cp")

  b <- summ_seeds(a) %>%
    mutate(run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
           model = "cp")

  cp.df <- rbind.data.frame(cp.df, a)
  cp.summ.df <- rbind.data.frame(cp.summ.df, b)
}

save(cp.df, cp.summ.df, file = "Ch1_movement_rates/sims_backup/datagen_cp.RData")

## PP Simulations -----------------------------------------------------------------------


ptpl %>%
  dplyr::select(fam_g, Bird_ID) %>%
  distinct(Bird_ID, .keep_all = TRUE) %>%
  group_by(fam_g) %>%
  summarise(n = n()) %>%
  right_join(., fam_moverate) %>%
  mutate(logpar = log(movrate),
         fitted_sd = logfit$estimate[2]) -> fam_moverate

pp_1 <- sort(round(fam_moverate$movrate,3))

logfit_fam <- fitdist(pp_1, distr = 'lnorm')

fam_moverate %>%
  arrange(., movrate) -> fam_moverate

# PP Generate data

kruns <- 1500
nseeds <- 5

pp.df <- NULL
pp.summ.df <- NULL

for(j in 1:7){
  for(k in 1:kruns){
    a <- sim_seeds(m.prms = fam_moverate$movrate[j], nseeds = nseeds) %>%
      mutate(fam_g = fam_moverate$fam_g[j],
             run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
             model = "pp")

    b <- summ_seeds(a) %>%
      mutate(fam_g = fam_moverate$fam_g[j],
             run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
             model = "pp")

    pp.df <- rbind.data.frame(pp.df, a)
    pp.summ.df <- rbind.data.frame(pp.summ.df, b)
  }
}

save(pp.df, pp.summ.df, file = "Ch1_movement_rates/sims_backup/datagen_pp.RData")

## NP Simulations -----------------------------------------------

logfit <- fitdist(indiv_moverate$movrate, distr = 'lnorm')

n.individuals <- 20
m_1 <- sort(round(rlnorm(n.individuals, meanlog = logfit$estimate[1], sdlog = logfit$estimate[2]),3))
m_2 <- sort(round(rlnorm(n.individuals, meanlog = logfit$estimate[1], sdlog = logfit$estimate[2]),3))
m_3 <- sort(round(rlnorm(n.individuals, meanlog = logfit$estimate[1], sdlog = logfit$estimate[2]),3))

# NP Generate data
m.data <- data.frame(m_1, m_2, m_3)
kruns <- 1000
nseeds <- 5

np.df <- NULL
np.summ.df <- NULL

for(m in 1:3){
  m.0 <- m.data[m]
  for(j in 1:n.individuals){
    for(k in 1:kruns){
      a <- sim_seeds(m.prms = m.0[j,], nseeds = nseeds) %>%
        mutate(indiv = as.factor(j),
               run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
               model = paste0("np_", m))

      b <- summ_seeds(a) %>%
        mutate(indiv = as.factor(j),
               run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
               model = paste0("np_", m))

      np.df <- rbind.data.frame(np.df, a)
      np.summ.df <- rbind.data.frame(np.summ.df, b)
    }
  }
}




save(np.df, np.summ.df, file = "Ch1_movement_rates/sims_backup/datagen_np.RData")

## PP w/indiv Simulations -----------------------------------------------

# Have a lognorm for each family group = using the mean from the family group, sd from population
# Sample 6 adults (That's how many adult aracari are supposed to be in family groups)

n.individuals <- 6

f.data <- data.frame(indiv = 1:n.individuals)

for(i in 1:7){
  samp <- rlnorm(n.individuals, meanlog = fam_moverate$logpar[i], sdlog = logfit$estimate[2])
  f.data <- cbind.data.frame(f.data, samp)
  names(f.data)[i+1] <- paste(fam_moverate$fam_g[i])
}

f.data <- f.data[2:8]

# NP Generate data
kruns <- 1000
nseeds <- 5

ppi.df <- NULL
ppi.summ.df <- NULL

for(f in 1:7){
  f.0 <- f.data[f]
  for(j in 1:n.individuals){
    for(k in 1:kruns){
      a <- sim_seeds(m.prms = f.0[j,], nseeds = nseeds) %>%
        mutate(indiv = as.factor(j),
               run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
               model = paste0("ppi_", f))

      b <- summ_seeds(a) %>%
        mutate(indiv = as.factor(j),
               run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
               model = paste0("ppi_", f))

      ppi.df <- rbind.data.frame(ppi.df, a)
      ppi.summ.df <- rbind.data.frame(ppi.summ.df, b)
    }
  }
}


save(ppi.df, ppi.summ.df, file = "Ch1_movement_rates/sims_backup/datagen_ppi.RData")















