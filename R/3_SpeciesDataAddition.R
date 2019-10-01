### SpeciesDataAddition
all_data <- readRDS("./Data/all_data.RDS")

species_to_remove <- species_list[species_list!=focal_species]

# Get only species columns
species_colnames <- colnames(all_data)[grepl(gsub(" ","_",focal_species),colnames(all_data))]
species_colnames_toKeep <- colnames(all_data)[!colnames(all_data) %in% 
                                                  grep(paste0(gsub(" ","_",species_to_remove), collapse = "|"), 
                                                       colnames(all_data), value = T)]

species_data <- all_data[,species_colnames_toKeep]

colnames(species_data)[colnames(species_data) %in% species_colnames] <- c("native","presence","introduced")

# species_colnames[1] is "native"
# species_colnames[2] is "presence"
# species_colnames[3] is "introduced"


species_introSites <- species_data[species_data$introduced == 1,"no_vatn_lnr"]

# Need to read in original raw data to get years of introduction
raw_data <- readRDS("./Data/allSpecies_GBIFDownload_Edited.RDS")
raw_data$year <- ifelse(raw_data$datasetKey == "e306fa70-381e-4330-8e68-1f447b46a850", 1918, raw_data$year)

raw_data$scientificNameShort <- word(raw_data$scientificName, 1,2, sep=" ")

get_years <- raw_data %>%
  as.data.frame() %>%
  filter(scientificNameShort == focal_species & vatnLnr %in% species_introSites & !is.na(year)) %>%
  dplyr::select(year, vatnLnr) %>%
  group_by(vatnLnr) %>%
  summarize(firstYear = min(year,na.rm=TRUE))

# I think the easiest thing to 
length(unique(get_years$firstYear))


# We split the data into two
# Get all lakes that can be used in model, get lakes with presence nearby
# and distance to nearest lake with presence by normal method.
# THen for the others we need to do things differently

# Assume that all presences in the native range were already there and not products of new intros
# Start from earliest year and work forward

# So we produce a data frame of our introductions, excluding everything in the native range
# and everything with an absence, plus everything with no identifiable year
species_intro_data <- species_data[species_data$introduced == 1 & 
                                     species_data$no_vatn_lnr %in% get_years$vatnLnr,]
species_intro_data$year <- ifelse(species_intro_data$no_vatn_lnr %in% get_years$vatnLnr,
                                  get_years$firstYear,0)

# We also produce a dataset with all presences, minus any introductions that have no identifiable
# year.
species_prese_data <- species_data[species_data$presence == 1 & !(species_data$introduced == 1 & 
  !(species_data$no_vatn_lnr %in% get_years$vatnLnr)),]
species_prese_data$year <- ifelse(species_prese_data$no_vatn_lnr %in% get_years$vatnLnr,
                                  get_years$firstYear,0)

# And produce a set of absences, since this is the easiest to calculate for
species_absen_data <- all_data[all_data[,species_colnames[1]] == 0 & all_data[,species_colnames[3]] == 0,
                               species_colnames_toKeep]

nn_absent <- get.knnx(species_prese_data[c("utm_x","utm_y")],species_absen_data[c("utm_x","utm_y")],k=1)
species_absen_data$dist_n_pop <- nn_absent$nn.dist[,1]


nn2_absent <- get.knnx(species_prese_data[c("utm_x","utm_y")],species_absen_data[c("utm_x","utm_y")],k=70)

number_nearby_pop <- nn2$nn.dist
number_nearby_pop[number_nearby_pop < 5000] <- 1
number_nearby_pop[number_nearby_pop >= 5000] <- 0
number_nearby_pop_vec <- rowSums(number_nearby_pop)







just_presences <- as.data.frame(all_lakes_sf) %>%
  filter(presence == 1) %>% 
  dplyr::select(waterBodyID,utm_x,utm_y)

all_lakes_df <- as.data.frame(all_lakes_sf)


nn <- get.knnx(just_presences[c("utm_x","utm_y")],all_lakes_df[c("utm_x","utm_y")],k=2)

# If distance = 0 then we it's referring to itself, if it's not 0, means that it's got the right lake.
dist_to_presence <- ifelse(nn$nn.dist[,1] == 0, nn$nn.dist[,2], nn$nn.dist[,1])

# Last thing is to get number of populations nearby

nn2 <- get.knnx(just_presences[c("utm_x","utm_y")],all_lakes_df[c("utm_x","utm_y")],k=70)



all_lakes_dataBuild$nearby_pops <- number_nearby_pop_vec # Nearby populations
all_lakes_dataBuild$pop_dist <- dist_to_presence         # Distance to nearest population

print("Calculated number of populations nearby and distance to nearest,
      compiling all data now")




print("Finished compiling data, find it in species folder under all_data.RDS")




