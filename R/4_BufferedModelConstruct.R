# Import raw data
raw_data <- readRDS(paste0("./Data/",gsub(' ','_',species_name),"/all_data.RDS"))

### First thing we need to do is get down to only the data inside a buffer

# Get a polygon around our presences
get_buffer <- raw_data %>%
  filter(presence == 1) %>%
  transmute(latitude = decimalLatitude, longitude = decimalLongitude)
get_buffer <- st_as_sf(get_buffer, coords = c("longitude", "latitude"), 
                       crs = 4326)
get_buffer <-  st_transform(get_buffer, 32633)

buf <- st_buffer(get_buffer, dist = 30000)
comb_bif <- st_union(buf)

points4buffer <- st_as_sf(raw_data, coords = c("decimalLongitude", "decimalLatitude"), 
                          crs = 4326) 
points4buffer <-  st_transform(points4buffer, 32633)

# Following produces a 1 if the value is in the buffer, 0 if it's not.
inOut_buffer <- st_intersects(points4buffer, comb_bif)
raw_data$inBuffer <- lengths(inOut_buffer)

raw_data_minimised <- raw_data %>%
  filter(area_km2 > size_threshold & use_in_model == 1 & inBuffer == 1)


# Just a couple of transformation functions. Greta, and most of these Bayesian packages, requires standardised data.
# I've log-transformed all fo the variables as well, except for temperature and HFP, as the rest were heavily
# right-skewed.
log_sc_gr <- function(vector) {new_vector <- as.numeric(scale(log(vector+1)))
return(new_vector)}
sc_gr <- function(vector) {new_vector <- as.numeric(scale(vector))
return(new_vector)}


# We now filter our data to get lakes larger than 5 ha
env_data_minimised <- raw_data_minimised %>%
  transmute(areaL = log_sc_gr(area_km2),
            dist_roadL = log_sc_gr(distance_to_road),
            temp = sc_gr(eurolst_bio10),
            pop_distL = log_sc_gr(pop_dist),
            SCI = log_sc_gr((perimeter_m/1000)/(2*sqrt(pi*area_km2))),
            HFP = sc_gr(HFP),
            n_pop = log_sc_gr(nearby_pops))

intro_data_minimised <- raw_data_minimised %>%
  dplyr::select(introduced)

id_data_minimised <- raw_data_minimised %>%
  dplyr::select(locationID)

# Turn them into Greta arrays
env_dataGr <- as_data(env_data_minimised)
intro_dataGr <- as_data(intro_data_minimised)

# A couple of values that will make things easier.
n_env <- ncol(env_dataGr)
n_sites <- nrow(env_dataGr)

# Create priors
alpha <- normal(0, 10)

beta <- normal(0, 10, dim=n_env)

# Defines our equation
linear_predictor <- alpha + env_dataGr %*% beta

# These two give our equation the required transformations
p <- ilogit(linear_predictor)
distribution(intro_dataGr) <- bernoulli(p)

# These define and visualise our model.
buffered_model <- model(beta)
# plot(buffered_model)

# The following then creates our MCMC draws
buffered_draws <- greta::mcmc(buffered_model,n_samples = 500, warmup = 500)
buffered_draws_extra <- extra_samples(buffered_draws,n_samples = 1000)


buffered_model_output <- list(draws = buffered_draws_extra, beta = beta, alpha = alpha, p = p)
saveRDS(buffered_model_output, file=paste0("./Data/",gsub(' ','_',species_name),"/buffered_model_output.RDS"))


buffered_model_data <- list(raw_data = raw_data_minimised, env_data = env_data_minimised, intro_data = intro_data_minimised, id_data = id_data_minimised, buffered_plot = comb_bif)
saveRDS(buffered_model_data, file=paste0("./Data/",gsub(' ','_',species_name),"/buffered_model_data.RDS"))



