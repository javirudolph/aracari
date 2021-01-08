## code to prepare `ptpl` dataset goes here

# Libraries needed
library(adehabitatLT)
library(lubridate)
library(dplyr)
library(tidyr)


# First, read in the raw tracking data
anim_data <- read.csv("data/raw_location_data.csv")
# Time is presented as a time ID
# Each time interval is equivalent to 15 minutes
# Interval 1 == 6AM, so if we start at midnight, 6am is 24 intervals of 15 minutes each. I'm setting 6am as the start just because I don't have actual data on it, but we know they started early morning and that's what the time IDs are for in Kimberly's work.


# create a column for number of minutes equivalent to the intervals

seconds <- (24*15 + anim_data$Time_ID * 15) * 60
td <- seconds_to_period(seconds)
anim_data$time <- sprintf("%02d:%02d:%02d", td@hour, minute(td), second(td))

date <- paste(anim_data$DATE, anim_data$time)
anim_data$date <- parse_date_time(date, "dmyHMS")
anim_data$burst <- paste(anim_data$Bird_ID, anim_data$DATE, sep="_")

# Only keep the columns we need for adehabitat

anim_data <- dplyr::select(anim_data, Bird_ID, burst,
                           date, X_Estimate, Y_Estimate)
names(anim_data) <- c("animal_id","burst_id", "date_time", "x", "y")

# we set the spatial coordinates for our data frame
coordinates(anim_data) <- c("x", "y")

#Assign a projection, I know it is UTM zone 18 because it is Ecuadorian Amazon
proj4string(anim_data) <- CRS("+proj=utm +zone=18 +datum=WGS84")

# Create ltraj object for work with adehabitat
# We set typeII = TRUE because we have variable time frames, although locations were attempted every 15 minutes, this didn't always happen.
ptpl_locs <- as.ltraj(coordinates(anim_data), burst = anim_data$burst_id,
                      date=anim_data$date_time, id=anim_data$animal_id,
                      typeII = TRUE)

# We can see the locations for all the individuals tracked
plot(ptpl_locs)

# The data frame now has the original xy locations, the time stamp, the individual
# But it has also added information about the displacement and net squared displacement.
# There is a high number of NA because we can only calculate distances moved from continous locations, happening during the same tracking session (which are a few hours per day)

ptpl <- ld(ptpl_locs)

# I want to include family groups in that data set, just in case we want to use that later
# The social groups are as follows: (1, 3, 5), (7), (13, 19), (22), (28), and (49, 84).
ptpl$id <- factor(ptpl$id, levels = c(1, 3, 5, 7, 13, 19, 22, 28, 49, 84, 10, 17, 20, 21, 24, 29, 30, 9))
ptpl %>%
  mutate(fam_g = ifelse(id %in% c(1,3,5), "f1",
                        ifelse(id == 7, "f2",
                               ifelse(id %in% c(13, 19), "f3",
                                      ifelse(id == 22, "f4",
                                             ifelse(id == 28, "f5",
                                                    ifelse(id %in% c(49,84),
                                                           "f6", "f7")
                                             )))))) -> ptpl
# Add new variables used for the selection later
# Some of the data was collected on a preliminary field season and it includes other time intervals. Here, we want to focus on the intervals that are multiples of 15, as that is how the majority of data is collected. From a biological stand point, a seed has a high probability of being regurgitated in 15-30 minutes.
# Filter by minimum number of observations. Removing the individuals with an insufficient number of observations. Since we will be fitting probability distributions, we should have at least 30 intervals per individual, so we would remove birds with ID = 17, 20, 24.
ptpl %>%
  mutate(T_minutes = dt/60,
         Bird_ID = id,
         mpm = dist/T_minutes,
         R2n = lead(R2n)) %>%
  filter(T_minutes != 0 & T_minutes %in% c(15, 30, 45, 60, 75, 90)) %>%
  group_by(Bird_ID) %>%
  add_tally() %>%
  filter(n >= 30) -> ptpl


# Last Step -----------------------------
usethis::use_data(ptpl, overwrite = TRUE)
