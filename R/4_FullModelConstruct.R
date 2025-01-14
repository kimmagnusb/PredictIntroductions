

# Import raw data
raw_data <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/species_model_data.RDS"))

# Just a couple of transformation functions. Greta, and most of these Bayesian packages, requires standardised data. I've log-transformed two of the variables as well, SCI and buffer_5000m_population_2006.
log_sc_gr <- function(vector) {new_vector <- as.numeric(scale(log(vector+1)))
return(new_vector)}
sc_gr <- function(vector) {new_vector <- as.numeric(scale(vector))
return(new_vector)}


# We now filter our data to get lakes larger than 2 ha
env_data <- raw_data %>%
  transmute(areaL = log_sc_gr(area_km2),
            dist_roadL = log_sc_gr(distance_to_road),
            temp = sc_gr(eurolst_bio10),
            dist_n_popL = log_sc_gr(dist_n_pop),
            sciL = log_sc_gr((perimeter_m/1000)/(2*sqrt(pi*area_km2))),
            HFP = sc_gr(HFP),
            no_n_popL = log_sc_gr(no_n_pop))

intro_data <- raw_data %>%
  dplyr::select(introduced)

id_data <- raw_data %>%
  dplyr::select(locationID)


# Turn them into Greta arrays
env_dataGr <- as_data(env_data)
intro_dataGr <- as_data(intro_data)

# A couple of values that will make things easier.
n_env <- ncol(env_dataGr)
n_sites <- nrow(env_dataGr)

print("Data is downloaded, parameters are ready")

# Create priors
alpha <- normal(0, 10)

beta <- normal(0, 10, dim=n_env)

# Defines our equation
linear_predictor <- alpha + env_dataGr %*% beta

# These two give our equation the required transformations
p <- ilogit(linear_predictor)
distribution(intro_dataGr) <- bernoulli(p)

# These define and visualise our model.
prelim_model <- model(beta)
# plot(prelim_model)

print("Constructed model, running draws now.")

# The following then creates our MCMC draws
whole_draws <- greta::mcmc(prelim_model,n_samples = 500, warmup = 500)
print("Finished running first draws, running extra draws now.")
whole_draws_extra <- extra_samples(whole_draws,n_samples = 1000)

print("Finished running draws, saving data now.")

whole_model_output <- list(draws = whole_draws_extra, beta = beta, alpha = alpha, p = p)
saveRDS(whole_model_output, file=paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/whole_model_output.RDS"))

whole_model_data <- list(raw_data = raw_data, env_data = env_data, intro_data = intro_data, id_data = id_data)
saveRDS(whole_model_data, file=paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/whole_model_data.RDS"))
