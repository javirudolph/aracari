# Revisiting the original data, because it's been too long.

library(readxl)
library(dplyr)
library(tidyr)
library(adehabitatLT)
library(lubridate)
library(stringr)


# Bring the datasets for the old birds
# From the Supplement table in Holbrook 2011, I'm using n.observations to identify these birds to the tag.

readxl::excel_sheets(path = "data-raw/bird1loc.xls")

# Bird 1 has 101 observations which goes along with Tag 84
bird1 <- readxl::read_excel("data-raw/bird1loc.xls", sheet = 1)
# There's no observations after 100
bird1 <- bird1[1:100,]

# This one has 108 so it's Tag 49
bird2 <- readxl::read_excel("data-raw/bird2loc.xls", sheet = 1)

# Bird 3 has 55 observations, so it is actually Tag28
bird3 <- readxl::read_excel("data-raw/bird3loc.xls", sheet = 1)

### The rest of the birds now
readxl::excel_sheets(path = "data-raw/Seed shadows PTPL 2005.xls")

# The rest of the birds

birds <- readxl::read_excel(path = "data-raw/Seed shadows PTPL 2005.xls", sheet = 6)
# Make the TIMEID numeric
birds %>%
  mutate(TIMEID = as.numeric(TIMEID)) -> birds


#############################################
# Make sure the date and time in the correct format
# I had done this before for the birds

seconds <- (24*15 + birds$TIMEID * 15) * 60
td <- seconds_to_period(seconds)
birds$TIME <- sprintf("%02d:%02d:%02d", td@hour, minute(td), second(td))

birds %>%
  rename(X_UTM = X_Estimate,
         Y_UTM = Y_Estimate) %>%
  mutate(DATE = as.POSIXct(paste(DATE, TIME)),
         TAG = `TRANS#`,
         burst = paste(TAG, DATE, sep = "_"),
         burst = sub(" .*", "", burst)) %>%
  dplyr:: select(X_UTM, Y_UTM, DATE, TAG, burst)-> xydf


# Set spatial coordinates for the dataset

coordinates(xydf) <- c("X_UTM", "Y_UTM")


#Assign a projection, I know it is UTM zone 18 because it is Ecuadorian Amazon
proj4string(xydf) <- CRS("+proj=utm +zone=18 +datum=WGS84")


# Create ltraj object for work with adehabitat
# We set typeII = TRUE because we have variable time frames, although locations were attempted every 15 minutes, this didn't always happen.
ptpl.ltraj <- as.ltraj(coordinates(xydf), burst = xydf$burst,
                       date = xydf$DATE, id = xydf$TAG, typeII = TRUE)


# We can see the locations for all the individuals tracked
plot(ptpl.ltraj)

# The data frame now has the original xy locations, the time stamp, the individual
# But it has also added information about the displacement and net squared displacement.
# There is a high number of NA because we can only calculate distances moved from continous locations, happening during the same tracking session (which are a few hours per day)

ptpl <- ld(ptpl.ltraj)





#######################################################

# Now, edit the other birds

bird1 %>%
  dplyr::select(X_UTM, Y_UTM, date, time) %>%
  mutate(date = as.Date(date),
         time = format(time, "%H:%M:%S"),
         DATE = paste(date, time)) %>%
  transmute(X_UTM,
            Y_UTM,
            DATE = as.POSIXct(DATE),
            TAG = as.character(84),
            burst = paste(TAG, DATE, sep = "_"),
            burst = sub(" .*", "", burst)) -> T84

coordinates(T84) <- c("X_UTM", "Y_UTM")
proj4string(T84) <- CRS("+proj=utm +zone=18 +datum=WGS84")
T84.traj <- as.ltraj(coordinates(T84), burst = T84$burst, date = T84$DATE, id = T84$TAG, typeII = TRUE)
T84.df <- ld(T84.traj)


bird2 %>%
  dplyr::select(X_UTM, Y_UTM, date, time) %>%
  mutate(date = as.Date(date),
         time = format(time, "%H:%M:%S"),
         DATE = paste(date, time)) %>%
  transmute(X_UTM,
            Y_UTM,
            DATE = as.POSIXct(DATE),
            TAG = as.character(49),
            burst = paste(TAG, DATE, sep = "_"),
            burst = sub(" .*", "", burst)) -> T49

coordinates(T49) <- c("X_UTM", "Y_UTM")
proj4string(T49) <- CRS("+proj=utm +zone=18 +datum=WGS84")
T49.traj <- as.ltraj(coordinates(T49), burst = T49$burst, date = T49$DATE, id = T49$TAG, typeII = TRUE)
T49.df <- ld(T49.traj)

bird3 %>%
  rename(X_UTM = X_Estimate,
         Y_UTM = Y_Estimate,
         date = date...6,
         time = time) %>%
  dplyr::select(X_UTM, Y_UTM, date, time) %>%
  mutate(date = as.Date(date),
         time = format(time, "%H:%M:%S"),
         DATE = paste(date, time)) %>%
  transmute(X_UTM,
            Y_UTM,
            DATE = as.POSIXct(DATE),
            TAG = as.character(28),
            burst = paste(TAG, DATE, sep = "_"),
            burst = sub(" .*", "", burst)) -> T28

coordinates(T28) <- c("X_UTM", "Y_UTM")
proj4string(T28) <- CRS("+proj=utm +zone=18 +datum=WGS84")
T28.traj <- as.ltraj(coordinates(T28), burst = T28$burst, date = T28$DATE, id = T28$TAG, typeII = TRUE)
T28.df <- ld(T28.traj)


#########################################################

ptpl <- ptpl %>%
  bind_rows(T28.df, T49.df, T84.df)

##################################
# MORE info

# The social groups are as follows: (1, 3, 5), (7), (13, 19), (22), (28), and (49, 84).
# Because birds 29 and 30 are from Yasuni Station, we consider those as one group as well.
# From the map figure in Holbrook 2011

ptpl %>%
  mutate(sgroup = case_when(
    id %in% c(1,3,5) ~ "G1",
    id %in% c(7) ~ "G2",
    id %in% c(13, 19) ~ "G3",
    id %in% c(22) ~ "G4",
    id %in% c(28) ~ "G5",
    id %in% c(49, 84) ~ "G6",
    id %in% c(29, 30) ~ "G7"
  )) -> ptpl

# Last Step -----------------------------
usethis::use_data(ptpl, overwrite = TRUE)






