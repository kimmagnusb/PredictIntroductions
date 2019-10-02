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


# We split the data into two
# Get all lakes that can be used in model, get lakes with presence nearby
# and distance to nearest lake with presence by normal method.
# THen for the others we need to do things differently

# Assume that all presences in the native range were already there and not products of new intros
# Start from earliest year and work forward

# So we produce a data frame of our introductions, excluding everything in the native range
# and everything with an absence, plus everything with no identifiable year
species_intro_data <- species_data %>% 
  filter(introduced == 1 & no_vatn_lnr %in% get_years$vatnLnr)
species_intro_data$year <- ifelse(species_intro_data$no_vatn_lnr %in% get_years$vatnLnr,
                                  get_years$firstYear,0)



# We also produce a dataset with all presences, minus any introductions that have no identifiable
# year.
species_prese_data <- species_data %>%
  filter(presence == 1 & 
           !(introduced ==1 & !(no_vatn_lnr %in% get_years$vatnLnr)))
species_prese_data$year <- ifelse(species_prese_data$no_vatn_lnr %in% get_years$vatnLnr,
                                  get_years$firstYear,0)

# And produce a set of absences, since this is the easiest to calculate for
species_absen_data <- species_data %>%
  filter(native == 0 & presence == 0)
  
nn_nearest_abs <- get.knnx(species_prese_data[c("utm_x","utm_y")],species_absen_data[c("utm_x","utm_y")],k=1)
species_absen_data$dist_n_pop <- nn_nearest_abs$nn.dist[,1]

nn_nearby_abs <- get.knnx(species_prese_data[c("utm_x","utm_y")],species_absen_data[c("utm_x","utm_y")],k=70)

number_nearby_pop <- nn_nearby_abs$nn.dist
number_nearby_pop[number_nearby_pop < 5000] <- 1
number_nearby_pop[number_nearby_pop >= 5000] <- 0
number_nearby_pop_vec <- rowSums(number_nearby_pop)

spatial_data <- data.frame(species_absen_data$no_vatn_lnr, nn_nearest_abs$nn.dist, number_nearby_pop_vec)
colnames(spatial_data) <- c("no_vatn_lnr","dist_n_pop","no_n_pop")


time_steps <- sort(unique(species_intro_data$year)) 
for (i in 1:length(time_steps)) {
  year_step <- time_steps[i]
  species_intro_historic <- species_intro_data %>% filter(year == year_step)
  species_presence_historic <- species_prese_data %>% filter(year < year_step)
  
  nn_nearest <- get.knnx(species_presence_historic[c("utm_x","utm_y")],species_intro_historic[c("utm_x","utm_y")],k=2)
  nearest_pop_vec <- ifelse(nn_nearest$nn.dist[,1]==0,nn_nearest$nn.dist[,2],nn_nearest$nn.dist[,1])
  
  nn_nearby <- get.knnx(species_presence_historic[c("utm_x","utm_y")],species_intro_historic[c("utm_x","utm_y")],k=70)
  number_nearby_pop <- nn_nearby$nn.dist
  number_nearby_pop[number_nearby_pop < 5000] <- 1
  number_nearby_pop[number_nearby_pop >= 5000] <- 0
  number_nearby_pop_vec <- rowSums(number_nearby_pop)
  
  spatial_data_timeStep <- data.frame(species_intro_historic$no_vatn_lnr, nearest_pop_vec, number_nearby_pop_vec)
  colnames(spatial_data_timeStep) <- c("no_vatn_lnr","dist_n_pop","no_n_pop")
  
  spatial_data <- rbind(spatial_data,spatial_data_timeStep)

  }


# Now that we have these variables, need to produce a data frame full of only data that can be used in
# our model (everything outside the native range)
species_model_data <- species_data %>%
  filter(native==0)

species_model_data <- merge(species_model_data,spatial_data, all.x=TRUE, by="no_vatn_lnr")

species_model_data <- species_model_data[complete.cases(species_model_data),]

# Find duplicates
# all_data %>% group_by(no_vatn_lnr) %>%
#   tally() %>%
#   filter(n>1)
# 
# all_data %>% filter(no_vatn_lnr == 39447)
# raw_data %>% filter(vatnLnr == 39447)

duplicated_vant_Lnrs <- unique(species_model_data$no_vatn_lnr[duplicated(species_model_data$no_vatn_lnr)])
if (length(duplicated_vant_Lnrs) != 0) {
  print(paste0("Warning: You have rows with duplicated Norwegian lake numbers. Lakes ",paste(duplicated_vant_Lnrs,collapse=", "),
               " are duplicated. Use the function display_duplicates to show them."),
        header = ngettext(n, "Warning message:\n", "Warning messages:\n"))
}


print("Calculated number of populations nearby and distance to nearest,
      compiling all data now")




print("Finished compiling data, find it in species folder under all_data.RDS")




