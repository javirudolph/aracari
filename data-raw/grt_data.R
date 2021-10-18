## code to prepare `grt_data` dataset goes here

# Libraries needed
library(dplyr)
library(stringr)


# Read the data csv

grt_raw <- read.csv("data/gut_ret_time.csv")


# This data is from field trials on gut retention times for four individual toucans
# The field trial seeds included a thread, and so we will keep only the data for regurgitated seeds with a thread.

grt_raw %>%
  dplyr::select(Name, X..Min, R.P, Seed.cond) %>%
  mutate(R.P = str_sub(R.P, start = 1, end = 1)) %>%
  drop_na(X..Min) %>%
  filter(R.P == 'R', Seed.cond == "THREAD") %>%
  transmute(Bird_ID = Name,
            grt = X..Min) -> grt_data

usethis::use_data(grt_data, overwrite = TRUE)
