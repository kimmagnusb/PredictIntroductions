#---------------------------------------------------------------------#
# 2. Import Biotic Data Relevant to All Species
#---------------------------------------------------------------------#


#---------------------------------------------------------------------#
# 2A. Import Biotic Data ####
#
# At the moment this downloads lake area, perimeter, distance to road
# and temeprature.
# 
#---------------------------------------------------------------------#

# This script downloads biotic data that is common to all lakes.

# Import data we have so far
all_lakes_sf <- readRDS(paste0("./Data/introductions.RDS"))

pg_user=rstudioapi::askForPassword("Wallace username")
pg_password=rstudioapi::askForPassword("password")

# Connect to the NOFA database
pool <- dbPool(drv = RPostgreSQL::PostgreSQL(), dbname = 'nofa',
               host = 'vm-srv-wallace.vm.ntnu.no', user = pg_user, password = pg_password,
               idleTimeout = 36000000,
               options="-c search_path=nofa"
)
con <- poolCheckout(pool)

# Get area, shoreline, distance to road
biotic_data_import <- tbl(con, "lake") %>%
  dplyr::select(id,area_km2, perimeter_m, distance_to_road) %>%
  filter(id %in% !!all_lakes_sf$waterBodyID) %>%
  collect()

# Get location data
biotic_data_import <- tbl(con, "lake") %>%
  dplyr::select(id,area_km2, perimeter_m, distance_to_road) %>%
  filter(id %in% !!all_lakes_sf$waterBodyID) %>%
  collect()

# Have to create another connection here, as the table is in a new folder
temp_pool <- dbPool(drv = RPostgreSQL::PostgreSQL(), dbname = 'nofa',
               host = 'vm-srv-wallace.vm.ntnu.no', user = pg_user, password = pg_password,
               idleTimeout = 36000000,
               options="-c search_path=environmental"
)
temp_con <- poolCheckout(temp_pool)

temp_data_import <- tbl(temp_con, "lake_EuroLST_BioClim") %>%
  dplyr::select(waterBodyID, eurolst_bio10) %>%
  filter(waterBodyID %in% !!all_lakes_sf$waterBodyID) %>%
  collect()

print("Downloaded biotic and temperature data")

#---------------------------------------------------------------------#
# 2B. Import HFP Data ####
#
# I've uplaoded the HFP (Human Footprint Index) data to the box folder 
# in raster form. Might need a better solution in the future. For more
# info on the HFP index click the link below.
# https://wcshumanfootprint.org/
# 
#---------------------------------------------------------------------#

# Get HFP data from box download

if (HFP_download == TRUE) {
  temp <- tempdir()
  download.file("https://ntnu.box.com/shared/static/ge0jpsojqdx3nvc9drr1997ix4l1aya8.rds",destfile = paste0(temp,"/HFP.rds"))
  HFP_raster <- readRDS(paste0(temp,"/HFP.rds"))
  saveRDS(HFP_raster,"./Data/HFP_raster.rds")
} else {
  HFP_raster <- readRDS("./Data/HFP_raster.rds")
}

# Need to get site points in a format we can compare to raster data
site_points <- as.data.frame(all_lakes_sf)[,c("decimalLatitude","decimalLongitude")]
coordinates(site_points) = ~ decimalLongitude+ decimalLatitude
crs(site_points) <-  "+init=epsg:4326 +proj=longlat"

site_points_conv <- spTransform(site_points, CRSobj = crs(HFP_raster))

HFPValues <- extract(HFP_raster, site_points_conv)

print("Downloaded HFP data")

#---------------------------------------------------------------------#
# 2B. Build consolidated data frame ####
#
# This brings all the above data together and gets rid of any duplicates
# we might have.
# 
#---------------------------------------------------------------------#


all_lakes_dataBuild <- all_lakes_sf
all_lakes_dataBuild$HFP <- HFPValues                     # Import HFP data

# Area, shoreline, distance to road
all_lakes_dataBuild <- merge(all_lakes_dataBuild, biotic_data_import,by.x = "waterBodyID", by.y= "id", all.x=TRUE)  

# Aaaand temperature
all_lakes_dataBuild <- merge(all_lakes_dataBuild, temp_data_import, by = "waterBodyID", all.x=TRUE)

# Create utm values for later use
all_lakes_utm <- unlist(all_lakes_dataBuild$geometry) %>%
  matrix(ncol=2,byrow=TRUE) %>% as.data.frame()

all_lakes_dataBuild$utm_x <- all_lakes_utm[,1]
all_lakes_dataBuild$utm_y <- all_lakes_utm[,2]

# Get rid of rows with NAs, get rid of lakes less than 2 hectares and get rid of geometry column.
all_lakes_final <- all_lakes_dataBuild %>%
  filter(area_km2 > size_threshold) %>%
  as.data.frame() %>%
  dplyr::select(-geometry)

all_lakes_final_complete <- all_lakes_final[complete.cases(all_lakes_final),]

# Need to sort out our duplicated lakes
duplicated_lakes_vatnLnr <- all_lakes_final_complete[duplicated(all_lakes_final_complete$no_vatn_lnr),"no_vatn_lnr"]
lake_dupes <- all_lakes_final_complete[all_lakes_final_complete$no_vatn_lnr %in% duplicated_lakes_vatnLnr,]

# We have 4 lakes which have been duplicated, for the following reasons:
# 39447 - 2 different lakes, can just keep them both in there
# 41167 - Same lake, just delete one
lakeRow_41167 <- which(all_lakes_final_complete$no_vatn_lnr==41167)[2]
all_lakes_final_complete <- all_lakes_final_complete[-lakeRow_41167,]
# 24083 - Same lake spread across county border, merge HFP and you're fine
HFP_lake24083 <- mean(all_lakes_final_complete[all_lakes_final_complete$no_vatn_lnr == 24083, "HFP"])
all_lakes_final_complete[all_lakes_final_complete$no_vatn_lnr == 24083, "HFP"] <- HFP_lake24083
lakeRow_24083 <- which(all_lakes_final_complete$no_vatn_lnr==24083)[2:3]
all_lakes_final_complete <- all_lakes_final_complete[-lakeRow_24083,]
# 7 - This is the problem one, the one with area of 19.6 is the one that should have presences
lakeRow_7 <- which(all_lakes_final_complete$no_vatn_lnr==7 & all_lakes_final_complete$area_km2 < 19)
all_lakes_final_complete <- all_lakes_final_complete[-lakeRow_7,]


saveRDS(all_lakes_final_complete, file=paste0("./Data/all_data.RDS"))
