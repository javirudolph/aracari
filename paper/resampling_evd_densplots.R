load("paper/fevd.RData")

#### THRESHOLD PLOTS
# NULL
# threshplot
par(mfrow=c(2,1))
nint <- 200
r <- quantile(null_dispersal$dispersal, probs=c(0.75,0.99))
u.i <- matrix(seq(r[1],r[2],, nint), ncol=1)
out <- null_threshplot
xlb <- "Threshold"
yl <- range(c(out[,c("low.t.scale","t.scale","up.t.scale")]), finite=TRUE)
plot(u.i, out[,"t.scale"], ylim=yl, xlab=xlb, ylab="reparameterized scale", type="b")
for(j in 1:nint) lines(c(u.i[j],u.i[j]), out[j,c("low.t.scale","up.t.scale")])

yl <- range(c(out[,c("low.shape","shape","up.shape")]), finite=TRUE)
plot(u.i, out[,"shape"], ylim=yl, xlab="Threshold", ylab="shape", type="b")
for(j in 1:nint) lines(c(u.i[j],u.i[j]), out[j,c("low.shape","up.shape")])


# mrl
par(mfrow=c(1,1))
out <- null_mrl
r <- range(null_dispersal$dispersal, finite=TRUE)
u.i <- matrix(seq(r[1], r[2] - 1,, nint), ncol=1)
xlab <- "Values"
yl <- range(c(out), finite=TRUE)
plot(u.i, out[,2], type="l", xlab=xlab, ylab="Mean Excess", ylim=yl)
lines(u.i, out[,1], lty=2, col="gray", lwd=1.5)
lines(u.i, out[,3], lty=2, col="gray", lwd=1.5)

slope <- out[1:nint-1,2]-out[2:nint,2]
yslope <- slope[1:198] - slope[2:199]
plot(u.i[1:length(yslope)], yslope, ylim = c(-5, 5))
abline(h = 0, col = "red")

# INDIV
# threshplot
par(mfrow=c(2,1))
nint <- 200
r <- quantile(indiv_dispersal$dispersal, probs=c(0.75,0.99))
u.i <- matrix(seq(r[1],r[2],, nint), ncol=1)
out <- indiv_threshplot
xlb <- "Threshold"
yl <- range(c(out[,c("low.t.scale","t.scale","up.t.scale")]), finite=TRUE)
plot(u.i, out[,"t.scale"], ylim=yl, xlab=xlb, ylab="reparameterized scale", type="b")
for(j in 1:nint) lines(c(u.i[j],u.i[j]), out[j,c("low.t.scale","up.t.scale")])

yl <- range(c(out[,c("low.shape","shape","up.shape")]), finite=TRUE)
plot(u.i, out[,"shape"], ylim=yl, xlab="Threshold", ylab="shape", type="b")
for(j in 1:nint) lines(c(u.i[j],u.i[j]), out[j,c("low.shape","up.shape")])


# mrl
par(mfrow=c(1,1))
out <- indiv_mrl
r <- range(indiv_dispersal$dispersal, finite=TRUE)
u.i <- matrix(seq(r[1], r[2] - 1,, nint), ncol=1)

yl <- range(c(out), finite=TRUE)
plot(u.i, out[,2], type="l", xlab=xlab, ylab="Mean Excess", ylim=yl)
lines(u.i, out[,1], lty=2, col="gray", lwd=1.5)
lines(u.i, out[,3], lty=2, col="gray", lwd=1.5)

slope <- out[1:nint-1,2]-out[2:nint,2]
yslope <- slope[1:198] - slope[2:199]
plot(u.i[1:length(yslope)], yslope, ylim = c(-5, 5))
abline(h = 0, col = "red")

#FAM
# threshplot
par(mfrow=c(2,1))
nint <- 200
r <- quantile(fam_dispersal$dispersal, probs=c(0.75,0.99))
u.i <- matrix(seq(r[1],r[2],, nint), ncol=1)
out <- fam_threshplot
xlb <- "Threshold"
yl <- range(c(out[,c("low.t.scale","t.scale","up.t.scale")]), finite=TRUE)
plot(u.i, out[,"t.scale"], ylim=yl, xlab=xlb, ylab="reparameterized scale", type="b")
for(j in 1:nint) lines(c(u.i[j],u.i[j]), out[j,c("low.t.scale","up.t.scale")])

yl <- range(c(out[,c("low.shape","shape","up.shape")]), finite=TRUE)
plot(u.i, out[,"shape"], ylim=yl, xlab="Threshold", ylab="shape", type="b")
for(j in 1:nint) lines(c(u.i[j],u.i[j]), out[j,c("low.shape","up.shape")])


# mrl
par(mfrow=c(1,1))
out <- fam_mrl
r <- range(fam_dispersal$dispersal, finite=TRUE)
u.i <- matrix(seq(r[1], r[2] - 1,, nint), ncol=1)

yl <- range(c(out), finite=TRUE)
plot(u.i, out[,2], type="l", xlab=xlab, ylab="Mean Excess", ylim=yl)
lines(u.i, out[,1], lty=2, col="gray", lwd=1.5)
lines(u.i, out[,3], lty=2, col="gray", lwd=1.5)

slope <- out[1:nint-1,2]-out[2:nint,2]
yslope <- slope[1:198] - slope[2:199]
plot(u.i[1:length(yslope)], yslope, ylim = c(-5, 5))
abline(h = 0, col = "red")

#*****************************************************************************************
load("paper/fevd.RData")
#*
# 3. Plot the distribution
#   a. Plot density with estimated parameters
#   b. Plot density with the 95% confidence intervals
#   c. Plot bootstrap? estimates: sample 10000 from the data, 1000 times and get those parameters. Make those density plots

ggplot(data = data.frame(x = c(100, 2000)), aes(x)) +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = null_scale, shape = null_shape), size = 1) -> null_base_plot


# Don't know if this is called bootstrapping

null_boot <- NULL
null_boot_probs <- NULL
for(i in 1:10){
  samp <- sample(null_dispersal$dispersal, 1000)
  fitD <- fevd(samp, threshold = null_thresh, type = "GP")
  scale <- fitD$results$par[1]
  shape <- fitD$results$par[2]
  out <- data.frame(boot = paste0("boot_", i), scale = scale, shape = shape)
  null_boot <- rbind.data.frame(null_boot, out)

  dist_range <- seq(null_thresh, 5000, by = 100)
  out2 <- data.frame(boot = paste0("boot_", i), distance = dist_range,
                     prob = pextRemes(null_fit, dist_range, lower.tail = FALSE))
  null_boot_probs <- rbind.data.frame(null_boot_probs, out2)
}

null_lines <- purrr::map(seq(1:1000), function(y) stat_function(fun = devd, args = list(type = "GP", scale = null_boot$scale[y], shape = null_boot$shape[i]), color = "grey"))

# CI
hist(scale(null_boot$scale))
hist(scale(null_boot$shape))

# Calculate confidence interval according to t-distribution
confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- data.frame("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}

norm.interval = function(data, variance = var(data), conf.level = 0.95) {
  z = qnorm((1 - conf.level)/2, lower.tail = FALSE)
  xbar = mean(data)
  sdx = sqrt(variance/length(data))
  c(xbar - z * sdx, xbar + z * sdx)
}

null_ci_scale <- confidence_interval(null_boot$scale, 0.95)
null_ci_shape <- confidence_interval(null_boot$shape, 0.95)

null_base_plot + null_lines +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = null_scale, shape = null_shape), size = 1) +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = null_ci_scale$lower, shape = null_ci_shape$lower),
                linetype = "dashed", size = 1) +
  stat_function(fun = devd, n = 101, args = list(type = "GP", scale = null_ci_scale$upper, shape = null_ci_shape$upper),
                linetype = "dashed", size = 1)


# 4. Calculate probability of getting those LDD events
#   a. make a figure with probability on the y axis, and distance in the x.

dist_range <- seq(null_thresh, 4000, by = 100)
null_probs <- data.frame(distance = dist_range,
                         prob = pextRemes(null_fit, dist_range, lower.tail = FALSE))


