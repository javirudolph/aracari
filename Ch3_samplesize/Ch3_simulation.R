##########################################
## DATA GENERATION CH3
#########################################

# I think the overall question for this is chapter is:
#         Can we accurately estimate the true proportion of rare events produced by a complex model using a simpler hierarchical model?

# Our complex model is one that considers individual variation in animal movement, with increasing mean and variance for specific groups
# We use a mixture modeling approach, where each component distribution is associated to an individual's state or age.
# We assume that an individual's stage will affect their movement.
# We assume, for example, very young birds might not move as far, as will nesting or molting birds, whereas adults will move more.

# The overall assumption for this model is that individuals vary in how much they move, which we determine based on step lengths only.
# With some groups having an overall larger variance and mean than other groups.
# The proportions, or weights, of the mixture, should follow a bell curve for now
# With the average movers being the majority of the group, and low or high movement individuals being less abundant.

# Libraries -----------------------------------------------
set.seed(20220201)

library(aracari)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(purrr)
library(fitdistrplus)

## Functions ----------------------------------------------------
# Building this function so that we get the parameters for a desired mean and standard deviation.
desired_mean_sd <- function(mean, sd){

  mu <- log(mean*2/(sqrt(mean*2+sd*2)))
  sigsq <- log(1+(sd*2/mean*2))

  return(data.frame(mu = mu, sigsq = sigsq))
}

exp.val <- function(mu, sigsq){
  exp(mu + sqrt(sigsq))
}


# MIXTURE MODEL --------------------------------------------------
## Mixture proportions ----------------------------------------
# We assume four categories of individuals with increasing movement.
pis <- c( 0.1, 0.2, 0.3, 0.4)
sum(pis)
hist(pis)

## Mixture components --------------------------------------
# Since movement lengths are only positive, we use a lognormal distribution to describe them

# CASE1: same variance for all, but change expected values

pars <- desired_mean_sd(mean = c(160, 300, 600, 1000), sd = c(90, 120, 160, 200))

# Special case I am considering where for each category the mean doubles, and then sd=mean*0.6
# pars <- desired_mean_sd(mean = c(150, 300, 600, 1200), sd = c(90, 180, 360, 720))

mus <- pars$mu
sigsqs <- pars$sigsq
dens_cols <- c("#264653", "#2a9d8f", "#f4a261", "#e76f51")

### PLOT the mixture components -------------------------------
lnorm_densities <- purrr::map(1:4, function(y) stat_function(fun = dlnorm,
                                                             args = list(meanlog = mus[y], sdlog = sigsqs[y]),
                                                             color = dens_cols[y], size=1))
ggplot() +
  lnorm_densities +
  theme_minimal() +
  labs(y = "Density") +
  lims(x = c(0, 150)) -> densities_plot


# SAMPLING --------------------------------------------------------------

samp.size <- 50000
all.samples <- rep(0,samp.size)
categ <- rep(0,samp.size)

for(i in 1:samp.size){

  which.cat <- sample(1:4, size=1, replace=TRUE, prob=pis)
  all.samples[i] <- rlnorm(n=1,meanlog=mus[which.cat], sdlog=sigsqs[which.cat])
  categ[i] <- which.cat
}



## PLOT the samples ---------------------------------------------
# This is our f(x) from the main text
data.tail <- data.frame(values = all.samples, y = 100) %>% arrange(desc(values)) %>% filter(values >=250)
head(data.tail)

ggplot(data.frame(all.samples), aes(x = all.samples)) +
  geom_histogram(bins = 100) +
  geom_point(data = data.tail[1:50,], aes(x = values, y = y), color = "black", alpha = 0.5) +
  labs(y = "Frequency", x = "Distance") +
  theme_minimal() -> sampleshist


# VIZ -----------------------------------
plot_grid(densities_plot, sampleshist, nrow=2, labels="AUTO")
#ggsave("Ch3_samplesize/TestFig.png")

summary(all.samples)

# FUNCTIONS ---------------------------------
source("Ch3_samplesize/Ch3_functions.R")

###################################
# Run the experiments with different sample sizes and tail thresholds

samp.sizes <- c(80, 200, 500, 800, 1000, 1600)
num.ns <- length(samp.sizes)
q.tests <- c(100, 250,500,1000)
num.qs <- length(q.tests)

all.sampsizes <- rep(samp.sizes,each=num.qs)
all.qtests <- rep(q.tests,num.ns)
ntests <- length(all.sampsizes)

cdfs.hat <- rep(0, ntests)


for(i in 1:ntests){

  ith.n   <- all.sampsizes[i]
  ith.samples <- data.frame(x = sample(all.samples, ith.n))
  mod1 <- lomax.glm(formula = ~1, my.dataf = ith.samples, response = ith.samples$x)
  ith.a <- mod1$alphas.hat[1]
  ith.k <- mod1$k.hat

  ith.q <- all.qtests[i]
  ith.cdf <- 1- lomax.cdf(x = ith.q, alpha = ith.a, k = ith.k)
  cdfs.hat[i] <- ith.cdf

}

true.cdfs <- rep(0,num.qs)
for(i in 1:num.qs){

  iq <- q.tests[i]

  true.cdfs[i] <- sum(all.samples>iq)/length(all.samples)

}

all.true.cdfs <- rep(true.cdfs,num.ns)


sim.test.df <- data.frame(all.sampsizes=all.sampsizes, all.qtests=all.qtests,cdfs.hat=cdfs.hat,
                          true.cdfs = all.true.cdfs)


# Now repeat but we do this 1000 times for each scenario of samples and tails.

ksamps <- 100

samp.sizes <- c(80, 200, 500, 800, 1000, 1600)
num.ns <- length(samp.sizes)
q.tests <- c(100, 250,500,1000)
num.qs <- length(q.tests)

all.sampsizes <- rep(samp.sizes,each=num.qs)
all.qtests <- rep(q.tests,num.ns)
ntests <- length(all.sampsizes)

cdfs.hat <- rep(0, ntests)

resamp.df <- NULL

for(k in 1:ksamps){

  for(i in 1:ntests){

    ith.n   <- all.sampsizes[i]
    ith.samples <- data.frame(x = sample(all.samples, ith.n))
    mod1 <- lomax.glm(formula = ~1, my.dataf = ith.samples, response = ith.samples$x)
    ith.a <- mod1$alphas.hat[1]
    ith.k <- mod1$k.hat

    ith.q <- all.qtests[i]
    ith.cdf <- 1- lomax.cdf(x = ith.q, alpha = ith.a, k = ith.k)
    cdfs.hat[i] <- ith.cdf

  }

  sim.test.df <- data.frame(all.sampsizes=all.sampsizes, all.qtests=all.qtests,cdfs.hat=cdfs.hat,
                            true.cdfs = all.true.cdfs, ksamp = paste0("samp", k))

  resamp.df <- rbind.data.frame(resamp.df, sim.test.df)
}


# Getting error: "
#
# Error in optim(par = init.betas, fn = ft.nllike, method = "BFGS", designmat = designmat,  :
#                  non-finite finite-difference value [1]

# So, change optimization method

lomax.glm <- function(formula=~1, my.dataf, response){

  Y           <- response
  n           <- length(Y)
  designmat   <- model.matrix(formula, data=my.dataf)
  nbetas      <- ncol(designmat)
  init.betas  <- rep(1.5,nbetas)

  opt.out <- optim(par=init.betas, fn=ft.nllike, method = "Nelder-Mead",
                   designmat=designmat, Y=Y)

  mles          <- opt.out$par
  nll.hat       <- opt.out$value
  BIC.mod       <- 2*nll.hat + nbetas*log(length(Y))
  Xbeta.hat     <- designmat%*%mles
  alphas.hat    <- exp(Xbeta.hat)
  sumlogapy.hat <- sum(log(alphas.hat+Y))
  k.hat         <- n/(sumlogapy.hat - sum(Xbeta.hat))

  out.list <- list(opt.out = opt.out, designmat=designmat,Y=Y, mles=mles, nll.hat=nll.hat, BIC.mod = BIC.mod,
                   alphas.hat=alphas.hat, k.hat=k.hat,data=my.dataf)

  return(out.list)

}

# Now repeat but we do this 1000 times for each scenario of samples and tails.

ksamps <- 1000

samp.sizes <- c(80, 200, 500, 800, 1000, 1600)
num.ns <- length(samp.sizes)
q.tests <- c(100, 250,500,1000)
num.qs <- length(q.tests)

all.sampsizes <- rep(samp.sizes,each=num.qs)
all.qtests <- rep(q.tests,num.ns)
ntests <- length(all.sampsizes)

cdfs.hat <- rep(0, ntests)

resamp.df <- NULL

for(k in 1:ksamps){

  for(i in 1:ntests){

    ith.n   <- all.sampsizes[i]
    ith.samples <- data.frame(x = sample(all.samples, ith.n))
    mod1 <- lomax.glm(formula = ~1, my.dataf = ith.samples, response = ith.samples$x)
    ith.a <- mod1$alphas.hat[1]
    ith.k <- mod1$k.hat

    ith.q <- all.qtests[i]
    ith.cdf <- 1- lomax.cdf(x = ith.q, alpha = ith.a, k = ith.k)
    cdfs.hat[i] <- ith.cdf

  }

  sim.test.df <- data.frame(all.sampsizes=all.sampsizes, all.qtests=all.qtests,cdfs.hat=cdfs.hat,
                            true.cdfs = all.true.cdfs, ksamp = paste0("samp", k))

  resamp.df <- rbind.data.frame(resamp.df, sim.test.df)
}

# I get results with warnings for using the Nelder-Mead estimation method. Says it's unreliable.

#save(resamp.df, file = "Ch3_samplesize/resamp.RData")

# Set working directory
load("resamp.RData")

resamp.df %>%
  mutate(diff = (cdfs.hat-true.cdfs),
         est.over.true = cdfs.hat/true.cdfs,
         sampsize = factor(all.sampsizes),
         tailthresh = factor(all.qtests))-> resamp.df

resamp.df %>%
  #filter(., all.sampsizes == "80") %>%
  ggplot(., aes(x = sampsize, y = est.over.true, color = tailthresh)) +
  #ggplot(., aes(x = sampsize, y = diff/true.cdfs, color = tailthresh)) +
  geom_boxplot() +
  labs(x = "Sample size", y = "y") +
  #scale_y_log10() +
  theme_minimal()

resamp.df %>%
  group_by(sampsize, tailthresh) %>%
  summarise(est.bias = signif(mean(est.over.true), 3)) %>%
  ggplot(., aes(x = sampsize, color = tailthresh, y = est.bias)) +
  geom_point() +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal()

