# Revisiting the original data, because it's been too long.

library(readxl)
library(dplyr)
library(tidyr)
library(adehabitatLT)
library(lubridate)
library(stringr)


# Bring the datasets for the old birds
# From the Supplement table in Holbrook 2011, I'm using n.observations to identify these birdst to the tag.

readxl::excel_sheets(path = "Ch2_distributions/Orig_data_KH/bird1loc.xls")

# Bird 1 has 101 observations which goes along with Tag 84
bird1 <- readxl::read_excel("Ch2_distributions/Orig_data_KH/bird1loc.xls", sheet = 1)
# There's no observations after 100
bird1 <- bird1[1:100,]

# This one has 108 so it's Tag 49
bird2 <- readxl::read_excel("Ch2_distributions/Orig_data_KH/bird2loc.xls", sheet = 1)

# Bird 3 has 55 observations, so it is actually Tag28
bird3 <- readxl::read_excel("Ch2_distributions/Orig_data_KH/bird3loc.xls", sheet = 1)

### The rest of the birds now
readxl::excel_sheets(path = "Ch2_distributions/Orig_data_KH/Seed shadows PTPL 2005.xls")

# The rest of the birds

birds <- readxl::read_excel(path = "Ch2_distributions/Orig_data_KH/Seed shadows PTPL 2005.xls", sheet = 6)
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
  dplyr:: select(X_UTM, Y_UTM, DATE, TAG, burst)-> all_tags

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

## Get them all into the same tag_df

all_tags %>%
  bind_rows(., T28, T49, T84) -> ptplUTM

#########################################################

# Set spatial coordinates for the dataset

coordinates(ptplUTM) <- c("X_UTM", "Y_UTM")

#Assign a projection, I know it is UTM zone 18 because it is Ecuadorian Amazon
proj4string(ptplUTM) <- CRS("+proj=utm +zone=18 +datum=WGS84")


# Create ltraj object for work with adehabitat
# We set typeII = TRUE because we have variable time frames, although locations were attempted every 15 minutes, this didn't always happen.
ptpl.ltraj <- as.ltraj(coordinates(ptplUTM), burst = ptplUTM$burst,
                       date = ptplUTM$DATE, id = ptplUTM$TAG, typeII = TRUE)


# We can see the locations for all the individuals tracked
plot(ptpl.ltraj)

# The data frame now has the original xy locations, the time stamp, the individual
# But it has also added information about the displacement and net squared displacement.
# There is a high number of NA because we can only calculate distances moved from continous locations, happening during the same tracking session (which are a few hours per day)

ptpl <- ld(ptpl.ltraj)

##################################
# MORE info

# The social groups are as follows: (1, 3, 5), (7), (13, 19), (22), (28), and (49, 84).
# From the map figure in Holbrook 2011

ptpl %>%
  mutate(group = case_when(
    id %in% c(1,3,5) ~ "G1",
    id %in% c(7) ~ "G2",
    id %in% c(13, 19) ~ "G3",
    id %in% c(22) ~ "G4",
    id %in% c(28) ~ "G5",
    id %in% c(49, 84) ~ "G6"
  )) -> ptpl

save(ptpl, file = "Ch2_distributions/Orig_data_KH/tidy_data.RData")








