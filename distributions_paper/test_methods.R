
# Needs an ltraj object for the location data
# so, run `format_raw_dat.Rmd` to get that one (up to line 143)

library(adehabitatLT)

is.regular(aracari_ltraj)

aracari_ltraj



refda <- as.POSIXlt(Sys.Date(), "UTC")

refda <- strptime("2003-06-01 00:00", "%Y-%m-%d %H:%M", tz = "UTC")


aracari_reg <- setNA(aracari_ltraj, dt = 900, refda, units = "sec")


plotltr(aracari_ltraj, "dt/60")

dev.off()

aracariI <- typeII2typeI(aracari_ltraj)
aracariI
plot(aracariI)

aracariIr <- redisltraj(aracariI, 100)


aracari_df %>%
  filter(R2n >0) %>%
  ggplot(., aes(x = R2n)) +
  geom_histogram() +
  facet_wrap(~id)

wawotest(aracari_ltraj)
testang.ltraj(aracari_ltraj, "relative")

aracari_ltraj[3] -> a
