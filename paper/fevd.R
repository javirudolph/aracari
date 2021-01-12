
# Fit extreme value distributions to the dispersal kernel data

load("paper/dispersal_kernels.RData")

library(dplyr)
library(ggplot2)
library(cowplot)
library(extRemes)

#################################################################################################################
# Just testing

null_sample <- sample(null_dispersal$dispersal, 1000)
hist(null_sample)
#hist(log(null_sample))
#qqnorm(log(null_sample))
qqnorm(null_sample)
# Qq plots don't show straight lines, not consistent with a normal distribution of the data,

a <- threshrange.plot(null_sample, r = c(300, 700), nint = 200)
b <- mrlplot(null_sample, xlim = c(300, 700), nint = 100)

c <- decluster(null_sample, threshold = 500)
plot(c)

fit_D <- fevd(null_sample, threshold = 500, type = "GP")
#fit_D$results
summary(fit_D)
plot(fit_D)
t <- plot(fit_D, type = "qq")
pextRemes(fit_D, c(500, 1000, 1500, 2000, 3000, 3500, 4000), lower.tail = FALSE)

scale <- fit_D$results$par[1]
shape <- fit_D$results$par[2]
dens <- revd(10000, scale = scale, shape = shape, type = "GP")
hist(dens)
d <- density(dens)
plot(d)

ggplot(data = data.frame(x = c(0.1, 6)), aes(x)) +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = 1, shape = 0)) +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = 1, shape = 0.5))
  ylab("") +
  scale_y_continuous(breaks = NULL)


lines <- purrr::map(seq(0, 1, 0.1), function(y) stat_function(fun = devd, args = list(type = "GP", scale = 1, shape = y), color = "grey"))

ggplot(data = data.frame(x = c(0.1, 6)), aes(x)) +
  ylab("") +
  scale_y_continuous(breaks = NULL) -> p

p + lines +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = 1, shape = 0), linetype = "dashed", size = 2) +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = 1, shape = 0.5))

################################################################################################################################

orig_thres <- 500
nint <- 200
r <- c(0, 700)
#*****************************************************************************************
# NULL MODEL

# 1. Determine threshold
null_threshplot <- threshrange.plot(null_dispersal$dispersal, r = r, type = "GP", nint = nint)
null_mrl <- mrlplot(null_dispersal$dispersal, nint = nint)

null_thresh <- orig_thres
# 2. Fit Generalized Pareto Distribution

null_fit <- fevd(null_dispersal$dispersal, threshold = null_thresh, type = "GP")
null_qq <- plot(null_fit, type = "qq")
null_summ <- summary(null_fit)
null_scale <- null_summ$par[1]
null_shape <- null_summ$par[2]
null_scale_se <- null_summ$se.theta[1]
null_shape_se <- null_summ$se.theta[2]

null_scale + 1.96 * null_summ$cov.theta[1,1]
null_scale - 1.96 * null_summ$cov.theta[1,1]
null_scale + null_scale_se
null_scale - null_scale_se

null_shape + 1.96 * null_summ$cov.theta[2,2]
null_shape - 1.96 * null_summ$cov.theta[2,2]
null_shape + null_shape_se
null_shape - null_shape_se

ci(null_fit, type = parameter)

# INDIVIDUAL MODEL
# 1. Determine threshold

indiv_threshplot <- threshrange.plot(indiv_dispersal$dispersal, r = r, type = "GP", nint = nint)
indiv_mrl <- mrlplot(indiv_dispersal$dispersal, nint = nint)

indiv_thresh <- orig_thres
# 2. Fit Generalized Pareto Distribution

indiv_fit <- fevd(indiv_dispersal$dispersal, threshold = indiv_thresh, type = "GP")
indiv_qq <- plot(indiv_fit, type = "qq")
indiv_summ <- summary(indiv_fit)
indiv_scale <- indiv_summ$par[1]
indiv_shape <- indiv_summ$par[2]
indiv_scale_se <- indiv_summ$se.theta[1]
indiv_shape_se <- indiv_summ$se.theta[2]

indiv_scale + 1.96 * indiv_summ$cov.theta[1,1]
indiv_scale - 1.96 * indiv_summ$cov.theta[1,1]
indiv_scale + indiv_scale_se
indiv_scale - indiv_scale_se

indiv_shape + 1.96 * indiv_summ$cov.theta[2,2]
indiv_shape - 1.96 * indiv_summ$cov.theta[2,2]
indiv_shape + indiv_shape_se
indiv_shape - indiv_shape_se


# FAMILY MODEL
# 1. Determine threshold

fam_threshplot <- threshrange.plot(fam_dispersal$dispersal, r = r, type = "GP", nint = nint)
fam_mrl <- mrlplot(fam_dispersal$dispersal, nint = nint)

fam_thresh <- orig_thres
# 2. Fit Generalized Pareto Distribution

fam_fit <- fevd(fam_dispersal$dispersal, threshold = fam_thresh, type = "GP")
fam_qq <- plot(fam_fit, type = "qq")
fam_summ <- summary(fam_fit)
fam_scale <- fam_summ$par[1]
fam_shape <- fam_summ$par[2]
fam_scale_se <- fam_summ$se.theta[1]
fam_shape_se <- fam_summ$se.theta[2]

fam_scale + 1.96 * fam_summ$cov.theta[1,1]
fam_scale - 1.96 * fam_summ$cov.theta[1,1]
fam_scale + fam_scale_se
fam_scale - fam_scale_se

fam_shape + 1.96 * fam_summ$cov.theta[2,2]
fam_shape - 1.96 * fam_summ$cov.theta[2,2]
fam_shape + fam_shape_se
fam_shape - fam_shape_se


# Make the table
evd_table <- data.frame(Model = c("Null", "Individual", "Family"),
                        Scale = c(paste0(signif(null_scale, 4), " \u00b1 ", signif(null_scale_se, 3)),
                                  paste0(signif(indiv_scale, 4), " \u00b1 ", signif(indiv_scale_se, 3)),
                                  paste0(signif(fam_scale, 4), " \u00b1 ", signif(fam_scale_se, 3))),
                        Shape = c(paste0(signif(null_shape, 4), " \u00b1 ", signif(null_shape_se, 3)),
                                  paste0(signif(indiv_shape, 4), " \u00b1 ", signif(indiv_shape_se, 3)),
                                  paste0(signif(fam_shape, 4), " \u00b1 ", signif(fam_shape_se, 3))))
evd_table

save.image("paper/fevd.RData")







