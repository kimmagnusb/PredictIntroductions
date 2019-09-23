
#---------------------------------------------------------------------
# 0. load libraries and set define some parameters 
#---------------------------------------------------------------------
library(rgbif)
library(dplyr)
library(sf)
library(stringr)

initate_download <- TRUE # should script initiate GBIF download? 
download_native_range <- TRUE # should script initiate download of native range?
species_list <- c("Coregonus lavaretus", # define species of interest
                  "Esox lucius", 
                  "Rutilus rutilus",
                  "Scardinius erythrophthalmus", 
                  "Perca fluviatilis",
                  "Oncorhynchus mykiss")


#---------------------------------------------------------------------
# 1. Get occurrence data from GBIF 
#---------------------------------------------------------------------

## 1A. Initiate Download

# If this is the first time downloading a species, you will need the following code to initiate the download.
# If you already have already initiated the download and the file is waiting, you can skip to 1B
#
# Warning: don't spin off to many downloads at once (i.e. re-run this pice of script again and again). 
# GBIF will only allow three concurrent download requests per. 


if (initiate_download == TRUE) {
  # set GBIF user account credentials 
  options(gbif_user=rstudioapi::askForPassword("my gbif username"))
  options(gbif_email=rstudioapi::askForPassword("my registred gbif e-mail"))
  options(gbif_pwd=rstudioapi::askForPassword("my gbif password"))
  
  keys <- as.integer()
  for(i in 1:length(species_list)){
    keys[i] <- name_suggest(q=species_list[i], rank='species')$key[1] 
  }
 
  
  # Take the key given by the code above and initiate the download using the following command
  download_key <- occ_download(
    paste0('taxonKey = ',paste(keys[1:3],collapse=",")),
    type = "or"
  ) %>% 
    occ_download_meta
  
  
  
  stop("Download being initiated. You will have to wait for download to complete in GBIF. Suggest that
       you wait 15 minutes then run script again with initiate_download set to FALSE. n that time,
       update the download_keysM table with new doanload codes and species info.")
  
  # store download_key for re-use 
  dir.create(".Data/", showWarnings=FALSE)  
  saveRDS(download_key,"./Data/download_key.rds")
}


  ## 1B. Get the data. If you already have the data, you can skip down to 1C.
  
 # Most of our downloads have already been initiated. The following table gives you download keys for each species
  # download_keys <- c("0009277-190813142620410","0005604-190813142620410","0009282-190813142620410",
  #                      "0012533-190813142620410", "0012534-190813142620410","0014766-190813142620410")
  # download_keysM <- as.data.frame(cbind(download_keys,c("Coregonus lavaretus", "Esox lucius", "Rutilus rutilus",
  #                                                        "Scardinius erythrophthalmus", "Perca fluviatilis","Oncorhynchus mykiss"),
  #                                       c("sik", "gjedde","mort","soerv","abbor",NA)))
  # dir.create("./Data",showWarnings = TRUE)
  # save(download_keysM,file="./Data/download_keysM.rda")

download_key <- readRDS("./Data/download_key.rds")

temp <- tempdir()
download.file(url=download_key$downloadLink,
              destfile=paste0(temp,"/tmp.zip"),
              quiet=TRUE, mode="wb")

# Unzip dowloaded occurrence file
species_distribution <- rio::import(unzip(paste0(temp,"/tmp.zip"),files="occurrence.txt"))
file_name <- paste0("./Data/",gsub(' ','_',species_name),"/GBIFDownload.RDS")
saveRDS(species_distribution,file_name)
unlink(temp)

print("Finished importing raw species data")

#---------------------------------------------------------------------
# 2. Download native distribution
#
# Load native distribution range and assign establishmentMeans status
# to occurrences accordingly 
#---------------------------------------------------------------------

## 2A. Get native distribution from https://doi.org/10.21400/1mwt3950 (NB: data are in EPSG:32633)
if (download_native_range == TRUE) {
  
  url <- "https://api.loke.aws.unit.no/dlr-gui-backend-resources-content/v2/contents/links/1efa6a5a-74b7-46d1-9bd6-49a7b3c58d030e379585-05d5-42ed-b00d-35d2d2b0f9bf1fc5a283-f739-4205-8c84-9748cb6d76ae"
  temp <- tempdir()
  download.file(url,paste0(temp,"/hk_native.zip"))
  unzip(paste0(temp,"/hk_native.zip"),exdir=temp, overwrite=TRUE)
  files <- unzip(paste0(temp,"/hk_native.zip"),list = TRUE)[,1]
  shapefile <- files[str_detect(files,".shp")]
  hk_distribution_map <-  st_read(paste0(temp,"/",shapefile))
  hk_distribution_map <- st_transform(hk_distribution_map, 4326) # reproject data to lat/long wgs84
  
  dir.create("Data/",showWarnings = FALSE)
  saveRDS(hk_distribution_map,"Data/native_distribution.rds")

}



## 2B: Organise point data into spatial format

# Need to match the point data to the occurrence data now
species_dist_short <- species_distribution %>%
  dplyr::select(gbifID,occurrenceID,catalogNumber,decimalLongitude, decimalLatitude,species,taxonKey,datasetKey, locality,municipality,county,countryCode,locationID,
                eventDate,year,month,day,samplingProtocol,eventID,fieldNumber,
                recordedBy,dynamicProperties,collectionCode,datasetName,license,institutionCode)
species_dist_short$longitude <- species_dist_short$decimalLongitude
species_dist_short$latitude <- species_dist_short$decimalLatitude

# Have to convert it using the sf package first. Standard CRS uses the EPSG 32633.
dist_sf <- st_as_sf(species_dist_short, coords = c("longitude", "latitude"), 
                    crs = 4326)
dist_sf <- st_transform(dist_sf, 32633)

print("Finished importing native species distribution")


## 2C. Download lake data and match to closest lake

# Following code bring the lake map down off box, so you'll need access to the box folder.
if (download_lakes == TRUE) {
  temp <- tempdir()
  download.file("https://ntnu.box.com/shared/static/6vu4de2birf9onmej2gorexwaxposuaa.zip",destfile = paste0(temp,"/NVE_innsjodatabasen.zip"))
  unzip(paste0(temp,"/NVE_innsjodatabasen.zip"),exdir=temp)
  lakes <- st_read(paste0(temp,"/Innsjo_Innsjo.shp"))
  saveRDS(lakes,"./Data/lake_polygons.rds")
} else {
  lakes <- readRDS("./Data/lake_polygons.rds")
}

# Now we join 2 columns to our data frame which give the number of the nearest lake and 
# the distance to that lake.
occ_with_lakes <- st_join(dist_sf, lakes, join = st_nearest_feature)
index <- st_nearest_feature(x = occ_with_lakes, y = lakes) # index of closest lake
closest_lakes <- lakes %>% dplyr::slice(index) # slice based on the index
dist_to_lake <- st_distance(x = occ_with_lakes, y= closest_lakes, by_element = TRUE) # get distance
occ_with_lakes$dist_to_lake <- as.numeric(dist_to_lake) # add the distance calculations to match data

# Get rid of all coords which are more than a certain distance from the nearest lake
occ_OKlakes <- occ_with_lakes %>% filter(dist_to_lake < dist_threshold)

# Get rid of duplicates
occ_OKlakes <- occ_OKlakes[!duplicated(occ_OKlakes$vatnLnr),]

## 2D. Now figure out whether or not it's in the native range
# Check intersection
if (native_range == TRUE) {
  fish_dist <- st_transform(fish_dist, 32633)
  inOut <- st_intersects(occ_OKlakes, fish_dist)

  # This gives us a list, which gives the status of occurrence for each lake. Since we want
  # no occurrence status, the list elements with length 0 are the ones that are not in the 
  # native range.
  occ_OKlakes$inOut <- lengths(inOut)
}
occ_OKlakes$inOut <- 0

introductions <- occ_OKlakes %>% filter(inOut==0)
presences <- occ_OKlakes

print("Finished matching species distribution to native range and lakes")

#### 3. Find absences and thus get a list of introductions that is usable #####

# only thing left to do is pull all absences from NOFA. For this, we need all locations 
# outside of the native distribution
if (all_lakes ==  TRUE) {
  
  pg_user=rstudioapi::askForPassword("Wallace username")
  pg_password=rstudioapi::askForPassword("password")
  
  # Connect to nofa, extract polygons
  pool <- dbPool(drv = RPostgreSQL::PostgreSQL(), dbname = 'nofa',
                 host = 'vm-srv-wallace.vm.ntnu.no', user = pg_user, password = pg_password,
                 idleTimeout = 36000000,
                 options="-c search_path=nofa"
  )
  con <- poolCheckout(pool)
  
  # Get lakes for everything
  all_lakes <- tbl(con,"location") %>%
    filter(countryCode =='NO') %>%
    select(locationID, waterBodyID, county, no_vatn_lnr, decimalLatitude, decimalLongitude) %>%
    collect()
  saveRDS(all_lakes,file="./Data/all_lakes.RDS")
} else {
  all_lakes <- readRDS("./Data/all_lakes.RDS")
}


# Dulicate coordinates for later so that we can have a geometry column (produced when
# using the sf package) and still have lat and long columns.
all_lakes$latitude <- all_lakes$decimalLatitude
all_lakes$longitude <- all_lakes$decimalLongitude

# Introduced introductions/presence columns
all_lakes$introduced <- ifelse(all_lakes$no_vatn_lnr %in% introductions$vatnLnr, 1, 0)
all_lakes$presence <- ifelse(all_lakes$no_vatn_lnr %in% presences$vatnLnr, 1, 0)

# Convert to sf file
all_lakes_sf <- st_as_sf(all_lakes, coords = c("longitude", "latitude"), 
                         crs = 4326)
all_lakes_sf <- st_transform(all_lakes_sf, 32633)

if (native_range == TRUE) {
  inOut_all <- st_intersects(all_lakes_sf, fish_dist)
  all_lakes_sf$native <- lengths(inOut_all)
} else {
  all_lakes_sf$native <- 0
}
# If they're not present in the native range then we can't sue them in the model.
# However we'll need them for the loops later, so keep them in
all_lakes_sf$use_in_model <- ifelse(!(all_lakes_sf$native == 1 & all_lakes_sf$presence == 0),1,0)

# Get rid of all lakes north of Trondelag
if (delete_north == TRUE) {
  all_lakes_sf <- all_lakes_sf %>%
    filter(!(county %in% c("Troms","Finnmark","Nordland")))
}
print("Finished editing data")

saveRDS(all_lakes_sf, file=paste0("./Data/",gsub(' ','_',species_name),"/introductions.RDS"))

