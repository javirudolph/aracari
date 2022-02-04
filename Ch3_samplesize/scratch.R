# Compare to the data I'm trying to recreate
load("Ch2_distributions/Orig_data_KH/tidy_data.RData")

realdist <- ptpl %>% drop_na(rel.angle) %>% transmute(realdist = dist/dt*900) %>% pull(realdist)
hist(realdist)
summary(realdist)



# DATA comparison -------------------------
## Get the data in -----------------------
load("Ch2_distributions/Orig_data_KH/tidy_data.RData")


ptpl %>%
  filter(dist != 0) %>%
  drop_na(sgroup) %>%
  group_by(id) %>%
  add_tally() %>%
  filter(n >= 30) %>%
  ungroup() %>%
  mutate(stplen = dist/(dt/900)) %>%
  dplyr::select(., dist, stplen, id, sgroup) -> ptpl

## Get some summary tables -------
# Overall summary of the step lengths
summary(ptpl$stplen)

ptpl %>%
  group_by(sgroup, id) %>%
  nest() %>%
  arrange(sgroup) %>%
  mutate(nobs = map_int(data, .f = nrow),
         summ = map(data, summary),
         stplen = map(data, pull, stplen)) %>%
  mutate(model = map(stplen, fitdist, distr = "lnorm")) -> a



ptpl %>%
  ggplot(., aes(x = stplen, group = id)) +
  facet_wrap(~sgroup) +
  geom_line(stat = "density") +
  theme_minimal()

load("Ch2_distributions/Ch2_fits.RData")

## lnorm fits to data --------------------------

reg_fits_info %>%
  filter(., dist=="lnorm" & model == "NP") %>%
  dplyr::select(estimate, param, ID) %>%
  pivot_wider(names_from = param, values_from = estimate) %>%
  mutate(color =MetBrewer::met.brewer("Hiroshige", n=12)) -> data_prms

lnorm_data_dens <- purrr::map(1:12, function(y) stat_function(fun = dlnorm,
                                                              args = list(meanlog = data_prms$meanlog[y], sdlog = data_prms$sdlog[y]),
                                                              color = data_prms$color[y]))
ggplot() +
  lnorm_data_dens +
  theme_minimal() +
  lims(x = c(0, 50)) -> data_dens_plot

