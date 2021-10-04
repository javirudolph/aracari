# Script to run simulations using parameters from aracari movement rates
# LIBRARIES --------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(fitdistrplus)
library(Hmisc)

set.seed(98)

# LOAD the data ------------------------------------------------------------------------

data(ptpl)

null_moverate <- data.frame(Bird_ID = unique(ptpl$Bird_ID), movrate = mean(ptpl$mpm))

indiv_moverate <- ptpl %>%
  group_by(Bird_ID, fam_g) %>%
  summarise(movrate = mean(mpm))

fam_moverate <- ptpl %>%
  group_by(fam_g) %>%
  summarise(movrate = mean(mpm))

ids <- ptpl %>% distinct(., Bird_ID, fam_g)

# Visualize the movement rates on top of a density distribution, so that we can "assume" these movement rates are sampled from that distribution. The case of theory is with a lognormal distribution.

indiv_moverate %>%
  arrange(., movrate) %>%
  ggplot(., aes(x = movrate)) +
  geom_histogram()


# Average movement rate for this 'population'
mean(indiv_moverate$movrate)
mean(ptpl$mpm)
sum.indivs <- summary(indiv_moverate$movrate)


# Color palette ----------------------------------------
my.cols1 <- c("#23262f","#717492","#b3a82a","#c94f21","#980012","#0d907a","#b9bec3")

# Vis Movement Rates ---------------------------------------------
# Fit a lognormal to the distribution of movement rates
# Keep in mind these are all males.

logfit <- fitdist(indiv_moverate$movrate, distr = 'lnorm')
# normfit <- fitdist(indiv_moverate$movrate, distr = "norm")
# Mean using the fit
exp(logfit$estimate[1])

# Visualize distribution and data points
vlines <- summary(indiv_moverate$movrate)[c(1,3,6)]

indiv_moverate %>%
  ggplot(., aes(x = movrate, y = 0.002)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_function(fun = dlnorm, args = list(meanlog = logfit$estimate[1], sdlog = logfit$estimate[2]),
                color = "black", alpha = 0.8, size = 1) +
  labs(x = "Movement Rate (meters/minute)", y = "Density") +
  geom_vline(xintercept = vlines, color = my.cols1[4], size = 1, lty = 2) +
  scale_x_continuous(limits = c(0, 70), breaks = c(0, vlines[1], 20, vlines[2], 40, vlines[3], 60),
                     labels = c(0, "Min.", 20, "Median.", 40, "Max.", 60)) +
  #scale_x_discrete(labels = c(0, vlines[1], 20, vlines[2], 40, vlines[3], 60)) +
  theme_bw()

ggsave2(filename = "Ch1_movement_rates/Figures/Lnorm_fit2data.png", width = 6, height = 4, units = "in")

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



# Complete pooling ---------------------------
# We've fitted this lognormal distribution to our movement rates of 12 individuals (12 movement rates). The average movement rate for this sample of the population (the 12 individuals) is the average of the movement rate for each individual:
mean(indiv_moverate$movrate)

# 27.94211 meters/minute

# In a complete pooling approach, we take this average movement rate and use it as the parameter for our exponential distribution, from which we draw the movement distance at each time step.

## Compare exponential movement distances to full dataset ----------------------------------

movrate_cp <- mean(indiv_moverate$movrate)
movrate_cp_ln <- as.numeric(exp(logfit$estimate[1]))

# test <- data.frame(mpm = rexp(500, rate = 1/movrate_cp))
#
# test %>%
#   ggplot(., aes(x = mpm)) +
#   geom_histogram(aes(y = ..density..)) +
#   stat_function(fun = dexp, args = list(rate = 1/movrate_cp), color = "red")

ptpl %>%
  ggplot(., aes(x = mpm)) +
  geom_histogram(aes(y = ..density..)) +
  stat_function(fun = dexp, args = list(rate = 1/movrate_cp), color = my.cols1[4], size = 1) +
  stat_function(fun = dexp, args = list(rate = 1/movrate_cp_ln), color = my.cols1[3], size = 1, alpha = 0.8) +
  labs(title = "Histogram of movement distances per minute from data set",
       x = "Movement distance per minute",
       y = "Density",
       caption = "Orange line is exponential density using average movement rate for population \n
       Yellow line shows exponential density using mean estimate from lognormal fit to movement rates") +
  theme_bw()

ggsave2(filename = "Ch1_movement_rates/Figures/Movement_dist_comparison_2exp.png", width = 6, height = 4, units = "in")

## Sample run -------

# Which estimate to use?
# So, we are saying here that complete pooling takes the average of all the data as one. We took all the data, then fitted a lognormal distribution. So, the parameter we use, is the average movement rate that the lognormal gives us. Since we are assuming that the lognormal distribution will describe this population's movement rates across individuals. So, we use the expected value (mean) of the fit to describe the average movement rate of the population.

movrate_cp_ln

# How many seeds to use? Landon used 100 because he was taking averages.
# I guess we are focusing on a very short time frame: what does the bird do after it eats those seeds at one tree, and how does it move until it drops them all? So, only use 5 seeds.

nseeds <- 5


# Let's visualize one run quickly

sim_movement(prm = movrate_cp_ln, plot.it = TRUE, return.data.frame = TRUE)
# Print the data frame with seed info
test_run <- sim_seeds(nseeds = nseeds, m.prms = movrate_cp_ln)
test_run

# Calculate summary for seeds in this run
summ_test_run <- summ_seeds(test_run)
# Visualize this one run:
test_run %>%
  ggplot(., aes(x = xloc, y = yloc)) +
  geom_path(color = "#2D3335") +
  geom_point(aes(x=0, y=0), color = "#87A986", size = 3) +
  geom_point(data = test_run %>% drop_na(s.id),
             aes(x = xloc, y = yloc)) +
  # geom_point(data = summ_test_run,
  #            aes(x = x, y = y), color = "#59879A") +
  labs(x = "x", y = "y",
       title = "Example of one simulation run") +
  theme_bw() -> A

test_run %>%
  ggplot(., aes(x = xloc, y = yloc)) +
  geom_path(color = "white", alpha = 0.1) +
  geom_point(aes(x=0, y=0), color = "#87A986", size = 3) +
  geom_segment(aes(x = summ_test_run$x, xend = 0, y = summ_test_run$y, yend = 0), color = "#87A986", lty = 2) +
  geom_point(data = test_run %>% drop_na(s.id),
             aes(x = xloc, y = yloc)) +
  geom_segment(data = test_run %>% drop_na(s.id),
               aes(x = xloc, y = yloc, xend = summ_test_run$x, yend = summ_test_run$y), lty = 2, color = "#E2BF80") +
  geom_point(data = summ_test_run,
             aes(x = x, y = y), color = "#E2BF80", size = 5) +
  labs(x = "x", y = "y",
       title = "Average dispersal and dispersion per run") +
  theme_bw() -> B

plot_grid(A, B, labels = "AUTO")


ggsave2(filename = "Ch1_movement_rates/Figures/Ex_one_sim_run.png", width = 6, height = 4, units = "in")

## CP Simulations -----------------------------------------------------------------------
# Generate seed dispersal data under a complete pooling scenario

kruns <- 10000
nseeds <- 5
m.prm <- movrate_cp_ln

df <- NULL
summ.df <- NULL

for(k in 1:kruns){
  a <- sim_seeds(m.prms = m.prm, nseeds = nseeds) %>%
    mutate(run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
           popu = "cp")

  b <- summ_seeds(a) %>%
    mutate(run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
           popu = "cp")

  df <- rbind.data.frame(df, a)
  summ.df <- rbind.data.frame(summ.df, b)
}



### Visualize dispersal and dispersion for complete pooling ---------------------------

# It's too many points, so sample 1000 to visualize.
summ.df %>%
  #group_by(., popu) %>%
  sample_n(., 1000) %>%
  ggplot(., aes(y = dsprsn, x = factor(popu))) +
  geom_boxplot() +
  geom_point(color = "#E2BF80", alpha = 0.2) +
  labs(title = "Dispersion per run" , y = "Seed dispersion in meters", x = "Complete pooling") +
  theme_bw() +
  theme(legend.position = "none") -> p1
# p1

summ.df %>%
  #group_by(., popu) %>%
  sample_n(., 1000) %>%
  ggplot(., aes(y = av.disp, x = factor(popu))) +
  geom_boxplot() +
  geom_point(color = "#87A986", alpha = 0.2) +
  labs(title = "Average dispersal per run", y = "Seed dispersal in meters", x = "Complete pooling") +
  theme_bw() +
  theme(legend.position = "none") -> p2
# p2

plot_grid(p1, p2)


ggsave2(filename = "Ch1_movement_rates/Figures/CP_disp_measures.png", width = 6, height = 4, units = "in")



### CP Kernel ------------------------------------------------------------

n.boots <- 10
samp.size <- 30
weib.boot.cp <- NULL

for(j in 1:n.boots){
  s.df <- df %>% drop_na(s.id) %>%
    sample_n(., samp.size) %>%
    mutate(disp = round(disp, digits = 2))

  g <- fitdist(s.df$disp, distr = "weibull", method = 'mle', lower = c(0,0))


  prms.weib <- data.frame(est.shape = as.numeric(g$estimate[1]),
                          est.scale = as.numeric(g$estimate[2]),
                          loglik = g$loglik,
                          popu = "cp")
  weib.boot.cp <- rbind.data.frame(weib.boot.cp, prms.weib %>% mutate(boot = j))
}

save.image(file = paste0("Ch1_movement_rates/workspace_", Sys.Date(), ".RData"))


weib.boot.cp %>%
  sample_n(., 100) %>%
  ggplot(., aes(x = factor(popu), y = est.shape)) +
  geom_violin() +
  # scale_color_manual(values = c("black", mycols)) +
  # geom_boxplot(width = 0.01) +
  # geom_point(color = "grey", alpha = 0.5) +
  # stat_summary(fun.data=mean_sdl, mult=1,
  #              geom="pointrange", color="black") +
  geom_jitter(position = position_jitter(0.1)) +
  labs(title = "Shape") +
  theme_bw() -> p1
# p1

weib.boot.cp %>%
  sample_n(., 100) %>%
  ggplot(., aes(x = factor(popu), y = est.scale)) +
  geom_violin() +
  # scale_color_manual(values = c("black", mycols)) +
  #geom_boxplot(width = 0.1) +
  # geom_point() +
  # stat_summary(fun.data=mean_sdl, mult=1,
  #              geom="pointrange", color="black") +
  geom_jitter(position = position_jitter(0.1)) +
  labs(title = "Scale") +
  theme_bw() -> p2

plot_grid(p1, p2)

# No pooling --------------------------------------------
# This one looks at individual variation.
# We will sample 20 individuals, 3 times from the lognormal distribution that describes movement rates for the population

n.individuals <- 20
m_1 <- sort(round(rlnorm(n.individuals, meanlog = logfit$estimate[1], sdlog = logfit$estimate[2]),3))
m_2 <- sort(round(rlnorm(n.individuals, meanlog = logfit$estimate[1], sdlog = logfit$estimate[2]),3))
m_3 <- sort(round(rlnorm(n.individuals, meanlog = logfit$estimate[1], sdlog = logfit$estimate[2]),3))


## Generate data -------------------------------------------------
m.data <- data.frame(m_1, m_2, m_3)
kruns <- 100
nseeds <- 5

df <- NULL
summ.df <- NULL

for(m in 1:3){
  m.0 <- m.data[m]
  for(j in 1:n.individuals){
    for(k in 1:kruns){
      a <- sim_seeds(m.prms = m.0[j,], nseeds = nseeds) %>%
        mutate(indiv = as.factor(j),
               run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
               popu = as.factor(m))

      b <- summ_seeds(a) %>%
        mutate(indiv = as.factor(j),
               run = factor(paste0("r_", k), levels = paste0("r_", 1:kruns)),
               popu = as.factor(m))

      df <- rbind.data.frame(df, a)
      summ.df <- rbind.data.frame(summ.df, b)
    }
  }
}

### Visualize -----
summ.df %>%
  group_by(., popu) %>%
  sample_n(., 100) %>%
  ggplot(., aes(y = dsprsn, x = factor(popu), color = factor(popu))) +
  geom_boxplot() +
  geom_point(color = "grey", alpha = 0.5) +
  labs(title = "Dispersion") +
  theme_bw() +
  theme(legend.position = "none") -> p1
# p1

summ.df %>%
  group_by(., popu) %>%
  sample_n(., 100) %>%
  ggplot(., aes(y = av.disp, x = factor(popu), color = factor(popu))) +
  geom_boxplot() +
  geom_point(color = "grey", alpha = 0.5) +
  labs(title = "Average dispersal per run") +
  theme_bw() +
  theme(legend.position = "none") -> p2
# p2

plot_grid(p1, p2)


n.boots <- 10
samp.size <- 30
weib.boot <- NULL

for(j in 1:n.boots){
  s.df <- df %>% drop_na(s.id) %>%
    group_by(popu) %>%
    sample_n(., samp.size) %>%
    mutate(disp = round(disp, digits = 2))

  # s.df %>%
  #   ggplot(., aes(x = disp, fill = popu)) +
  #   geom_histogram()

  # s.df <- df %>% drop_na(s.id)

  weib.fits <- NULL

  for(i in 1:3){
    dat <- s.df %>%
      filter(popu == i)

    g <- fitdist(dat$disp, distr = "weibull", method = 'mle', lower = c(0,0))
    #g <- fitdistr(dat$disp, densfun = "weibull", lower = c(0,0))
    prms.weib <- data.frame(est.shape = as.numeric(g$estimate[1]),
                            est.scale = as.numeric(g$estimate[2]),
                            loglik = g$loglik,
                            popu = i)
    weib.fits <- rbind.data.frame(weib.fits, prms.weib)
    # plot(g)
  }

  weib.boot <- rbind.data.frame(weib.boot, weib.fits %>% mutate(boot = j))
}


weib.boot %>%
  group_by(popu) %>%
  #sample_n(., 100) %>%
  ggplot(., aes(x = factor(popu), y = est.shape, color = factor(popu))) +
  geom_violin() +
  # scale_color_manual(values = c("black", mycols)) +
  # geom_boxplot(width = 0.01) +
  # geom_point(color = "grey", alpha = 0.5) +
  # stat_summary(fun.data=mean_sdl, mult=1,
  #              geom="pointrange", color="black") +
  geom_jitter(position = position_jitter(0.1)) +
  labs(title = "Shape") +
  theme_bw() -> p1
# p1

weib.boot %>%
  group_by(popu) %>%
  #sample_n(., 100) %>%
  ggplot(., aes(x = factor(popu), y = est.scale, color = factor(popu))) +
  geom_violin() +
  # scale_color_manual(values = c("black", mycols)) +
  #geom_boxplot(width = 0.1) +
  # geom_point() +
  # stat_summary(fun.data=mean_sdl, mult=1,
  #              geom="pointrange", color="black") +
  geom_jitter(position = position_jitter(0.1)) +
  labs(title = "Scale") +
  theme_bw() -> p2

plot_grid(p1, p2)



# Simulate movement rates -------------------------------------------------------------
n.individuals <- 20
m_1 <- sort(round(rlnorm(n.individuals, meanlog = logfit$estimate[1], sdlog = logfit$estimate[2]),3))

logfitdraw <- data.frame(indv = 1:20, movrate = m_1)

mycols <- viridis::magma(20)

indiv_moverate %>%
  ggplot(., aes(x = movrate, y = 0)) +
  geom_point(size = 5) +
  stat_function(fun = dlnorm, args = list(meanlog = logfit$estimate[1], sdlog = logfit$estimate[2]), color = "black", alpha = 0.8) +
  geom_point(data = logfitdraw, aes(x = movrate, color = factor(indv), y = 0.01), size = 5) +
  # scale_color_viridis_d(option="magma") +
  scale_color_viridis_d() +
  xlim(0, 60) +
  labs(title = "Mov rate sampling lnorm") +
  theme_bw() -> plot.logfitdraw
plot.logfitdraw

### Fist population--------------------------------------------------------------------

expcurves <- purrr::map(seq(1:n.individuals), function(y)
  stat_function(fun = dexp, args = list(rate = 1/logfitdraw$movrate[y]), color = mycols[y], alpha = 0.8))

logfitdraw %>%
  ggplot() +
  expcurves +
  xlim(-1, 125) +
  labs(title = "MD dist exp", x = 'movement distance') +
  theme_bw() -> expdraws
expdraws

cowplot::plot_grid(plot.logfitdraw, expdraws)

