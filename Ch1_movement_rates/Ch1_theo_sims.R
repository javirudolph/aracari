# Script to run simulations for theoretical populations with three levels of variation in movement rates

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(fitdistrplus)
library(Hmisc)

set.seed(98)

# Three populations ---------------------------------------------------------------------
av.mov.rate <- 1
sd.mov.rate <- c(0.3, 0.5, 1)
max.x <- 10

mycols <- c("red", "green", "blue")
#create density plots
curve(dlnorm(x, meanlog=av.mov.rate, sdlog=sd.mov.rate[1]),
      from=0, to=max.x,
      col=mycols[1],
      main = 'Log Normally distributed movement rates for three populations', #add title
      ylab = 'Density', #change y-axis label
)

curve(dlnorm(x, meanlog=av.mov.rate, sdlog=sd.mov.rate[2]),
      from=0, to=max.x,
      col=mycols[2], add=TRUE)

curve(dlnorm(x, meanlog=av.mov.rate, sdlog=sd.mov.rate[3]),
      from=0, to=max.x,
      col=mycols[3], add=TRUE)

#add legend
legend("topright", legend=c("sdlog=.3", "sdlog=.5", "sdlog=1"),
       col=mycols, lty=1, cex=1.2)

## Simulate movement rates -------------------------------------------------------------
n.individuals <- 30
m_1 <- sort(round(rlnorm(n.individuals, meanlog = av.mov.rate, sdlog = sd.mov.rate[1]),3))
m_2 <- sort(round(rlnorm(n.individuals, meanlog = av.mov.rate, sdlog = sd.mov.rate[2]),3))
m_3 <- sort(round(rlnorm(n.individuals, meanlog = av.mov.rate, sdlog = sd.mov.rate[3]),3))

## Visualize individual differences ----------------------------------------------------

par(mfrow=c(1,3))
mycols <- rainbow(n.individuals, s=1)
mycols <- colorRampPalette(c("#F46036", "#5B85AA", "#414770", "#372248", "#171123"))(n.individuals)
mycols <- colorRampPalette(c("#ffba08", "#d00000", "#03071e"))(n.individuals)

range <- c(-0.1,max.x+5)
#create density plots
### Fist population--------------------------------------------------------------------
curve(dexp(x, 1/m_1[1]), #notice rate for exponential is 1/movement rate.
      from=range[1], to=range[2],
      col=mycols[1],
      main = 'Distribution of movement lengths', #add title
      ylab = 'Density', #change y-axis label
      xlab = "x"
)
for(i in 2:n.individuals){
  curve(dexp(x, 1/m_1[i]),
        from=range[1], to=range[2],
        col=mycols[i], add=TRUE)
}

#add legend
legend("topright", legend= m_1,
       col=mycols, lty=1, cex=1.2, lwd = 2)


### Second population----------------------------------------------------------------
curve(dexp(x, 1/m_2[1]), #notice rate for exponential is 1/movement rate.
      from=range[1], to=range[2],
      col=mycols[1],
      main = 'Distribution of movement lengths', #add title
      ylab = 'Density', #change y-axis label
      xlab = "x"
)
for(i in 2:n.individuals){
  curve(dexp(x, 1/m_2[i]),
        from=range[1], to=range[2],
        col=mycols[i], add=TRUE)
}

#add legend
legend("topright", legend= m_2,
       col=mycols, lty=1, cex=1.2, lwd = 2)



### Third population ---------------------------------------------------------------
curve(dexp(x, 1/m_3[1]), #notice rate for exponential is 1/movement rate.
      from=range[1], to=range[2],
      col=mycols[1],
      main = 'Distribution of movement lengths', #add title
      ylab = 'Density', #change y-axis label
      xlab = "x"
)
for(i in 2:n.individuals){
  curve(dexp(x, 1/m_3[i]),
        from=range[1], to=range[2],
        col=mycols[i], add=TRUE)
}

#add legend
legend("topright", legend= m_3,
       col=mycols, lty=1, cex=1.2, lwd = 2)

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

# Simulate seed dispersal ------------------------------------------------------------
## Animal movement first ---------------------------------------------------------------

## Don't run this, too heavy of a plot and you can't really see anything.
# nruns <- 3
# m.data <- data.frame(m_1, m_2, m_3)

# df <- NULL
#
# for(m in 1:3){
#   m.0 <- m.data[m]
#   for(j in 1:n.individuals){
#     for(k in 1:nruns){
#       a <- sim_movement(m.0[j,], plot.it = FALSE, return.data.frame = TRUE) %>%
#         mutate(indiv = as.factor(j),
#                run = paste0("r_", k),
#                pop.id = paste0("p_",m))
#       df <- rbind.data.frame(df, a)
#     }
#   }
#
# }
#
#
# df %>%
#   ggplot(., aes(x = xloc, y=yloc, group = run, color = indiv)) +
#   facet_wrap(~pop.id) +
#   geom_path() +
#   scale_color_manual(values = mycols) +
#   theme_bw() +
#   theme(legend.position = "bottom") +
#   NULL

## Generate seed dispersal data -----------------------------------------------
m.data <- data.frame(m_1, m_2, m_3)
kruns <- 100
nseeds <- 20

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

save.image(file = paste0("Ch1_movement_rates/workspace_", Sys.Date(), ".RData"))

### Check with one individual ---------------------------------------------------------

df %>%
  filter(popu == 1, indiv == 1) -> df.ex

df.ex %>%
  ggplot(., aes(x = xloc, y=yloc,
                group = run,
                color = indiv)) +
  # facet_wrap(~run) +
  geom_path() +
  scale_color_manual(values = mycols) +
  geom_point(data = df.ex %>% drop_na(s.id),
             aes(x = xloc, y = yloc)) +
  geom_point(data = summ.df %>%
               filter(popu == 1, indiv == 1),
             aes(x = x, y=y), color = "black") +
  theme_bw() +
  labs(title = "Animal seed droppings", caption = "Lines show an individual's trajectory.\n Dots show all seed droppings.\n Black dots show average location of seed per run") +
  theme(legend.position = "none") -> p1
# p1

df.ex %>%
  ggplot(., aes(x = xloc, y = yloc)) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  geom_bin2d() +
  # geom_hex() +
  # scale_fill_continuous(type = "viridis") +
  scale_fill_gradient(low = "grey", high = "black") +
  # geom_point(data = df %>% drop_na(s.id),
  #            aes(x = xloc, y = yloc), color = mycols[ex.i], alpha = 0.3) +
  # geom_point(data = summ.df, aes(x = x, y=y), color = "black") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Density of animal movement") -> p2
# p2

cowplot::plot_grid(p1, p2)

## Dispersion and dispersal measures across populations ------------------------------------------------------

# It's too many points, so sample 1000 to visualize.
summ.df %>%
  group_by(., popu) %>%
  sample_n(., 1000) %>%
  ggplot(., aes(y = dsprsn, x = factor(popu), color = factor(popu))) +
  geom_boxplot() +
  geom_point(color = "grey", alpha = 0.5) +
  labs(title = "Dispersion") +
  theme_bw() +
  theme(legend.position = "none") -> p1
# p1

summ.df %>% c
group_by(., popu) %>%
  sample_n(., 1000) %>%
  ggplot(., aes(y = av.disp, x = factor(popu), color = factor(popu))) +
  geom_boxplot() +
  geom_point(color = "grey", alpha = 0.5) +
  labs(title = "Average dispersal per run") +
  theme_bw() +
  theme(legend.position = "none") -> p2
# p2

plot_grid(p1, p2)

##
