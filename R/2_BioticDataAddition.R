

# Import data we have so far
all_lakes_sf <- readRDS(paste0("./Data/",gsub(' ','_',species_name),"/introductions.RDS"))


pg_user=rstudioapi::askForPassword("Wallace username")
pg_password=rstudioapi::askForPassword("password")

# Connect to nofa, extract polygons
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

# Get temp
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

# Get HFP

HFP_raster <- raster('./Data/HFP/HFP1993.tif')

site_points <- as.data.frame(all_lakes_sf)[,c("decimalLatitude","decimalLongitude")]
coordinates(site_points) = ~ decimalLongitude+ decimalLatitude
crs(site_points) <-  "+init=epsg:4326 +proj=longlat"

site_points_conv<-spTransform(site_points, CRSobj = crs(HFP_raster))

HFPValues <- extract(HFP_raster, site_points_conv)
summary(HFPValues)

print("Downloaded HFP data")

# Get distance to closest poopulations
# Need to make a table with just presences
# Create utm values for all_lakes
all_lakes_utm <- unlist(all_lakes_sf$geometry) %>%
  matrix(ncol=2,byrow=TRUE) %>% as.data.frame()

all_lakes_sf$utm_x <- all_lakes_utm[,1]
all_lakes_sf$utm_y <- all_lakes_utm[,2]

just_presences <- as.data.frame(all_lakes_sf) %>%
  filter(presence == 1) %>% 
  dplyr::select(waterBodyID,utm_x,utm_y)

all_lakes_df <- as.data.frame(all_lakes_sf)


nn <- get.knnx(just_presences[c("utm_x","utm_y")],all_lakes_df[c("utm_x","utm_y")],k=2)

# If distance = 0 then we it's referring to itself, if it's not 0, means that it's got the right lake.
dist_to_presence <- ifelse(nn$nn.dist[,1] == 0, nn$nn.dist[,2], nn$nn.dist[,1])

# Last thing is to get number of populations nearby

nn2 <- get.knnx(just_presences[c("utm_x","utm_y")],all_lakes_df[c("utm_x","utm_y")],k=70)

number_nearby_pop <- nn2$nn.dist
number_nearby_pop[number_nearby_pop < 5000] <- 1
number_nearby_pop[number_nearby_pop >= 5000] <- 0
number_nearby_pop_vec <- rowSums(number_nearby_pop)

all_lakes_dataBuild <- all_lakes_sf
all_lakes_dataBuild$HFP <- HFPValues                     # Human footprint
all_lakes_dataBuild$nearby_pops <- number_nearby_pop_vec # Nearby populations
all_lakes_dataBuild$pop_dist <- dist_to_presence         # Distance to nearest population

print("Calculated number of populations nearby and distance to nearest,
      compiling all data now")


# Area, shoreline, distance to road
all_lakes_dataBuild <- merge(all_lakes_dataBuild, biotic_data_import,by.x = "waterBodyID", by.y= "id", all.x=TRUE)  

# Aaaand temperature
all_lakes_dataBuild <- merge(all_lakes_dataBuild, temp_data_import, by = "waterBodyID", all.x=TRUE)

# Get complete cases
all_lakes_final <- all_lakes_dataBuild %>%
  as.data.frame() %>%
  dplyr::select(-geometry)

all_lakes_final_complete <- all_lakes_final[complete.cases(all_lakes_final),]

print("Finished compiling data, find it in species folder under all_data.RDS")


saveRDS(all_lakes_final_complete, file=paste0("./Data/",gsub(' ','_',species_name),"/all_data.RDS"))


