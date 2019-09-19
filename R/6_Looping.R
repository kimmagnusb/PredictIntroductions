
# Import data and draws and define components

if (use_buffered_model == TRUE) {
  model_output <- readRDS(paste0("./Data/",gsub(' ','_',species_name),"/buffered_model_output.RDS"))
  model_data_extra <- readRDS(paste0("./Data/",gsub(' ','_',species_name),"/buffered_model_data.RDS"))
} else {
  model_output <- readRDS(paste0("./Data/",gsub(' ','_',species_name),"/whole_model_output.RDS"))
  model_data_extra <- readRDS(paste0("./Data/",gsub(' ','_',species_name),"/whole_model_data.RDS"))
}

all_data <- readRDS(paste0("./Data/",gsub(' ','_',species_name),"/all_data.RDS"))

draws <- model_output$draws
beta <- model_output$beta
p <- model_output$p
alpha <- model_output$alpha

# Whatever data/model we are using for out parameters, we are now applying the parameters to all the data.
# That means we need to scale all the data using whatever means and SDs with which we scaled the data that
# was used for the moel.
data_4scaling <- model_data_extra$raw_data
#env_var <- c("HFP", "nearby_pops", "pop_dist", "area_km2", "perimeter_m", "distance_to_road", "eurolst_bio10")
data_4scaling <- data_4scaling %>%
  transmute(areaL = log(area_km2+1),
            dist_roadL = log(distance_to_road+1),
            temp = eurolst_bio10,
            pop_distL = log(pop_dist+1),
            SCI = log(((perimeter_m/1000)/(2*sqrt(pi*area_km2)))+1),
            HFP = HFP,
            n_pop = log(nearby_pops+1))

# Get the means and sds to work with
means_relData <- apply(data_4scaling, 2, mean)
SDs_relData <- apply(data_4scaling, 2, sd)

# Now we scale our full data based on these means
env_data <- all_data %>%
  transmute(areaL = log(area_km2+1),
            dist_roadL = log(distance_to_road+1),
            temp = eurolst_bio10,
            pop_distL = log(pop_dist+1),
            SCI = log(((perimeter_m/1000)/(2*sqrt(pi*area_km2)))+1),
            HFP = HFP,
            n_pop = log(nearby_pops+1))

env_data_takeMeans <- t(apply(env_data, 1, "-", means_relData))
env_data_scaled <- t(apply(env_data_takeMeans, 1, "/", SDs_relData))

print("Data is scaled.")

### Get betas and alphas
calc_beta <- apply(as.matrix(calculate(beta,draws)),2,mean)
calc_alpha <- apply(as.matrix(calculate(alpha,draws)),2,mean)

presences <- all_data$presence

# Get full index of nearby lakes that have a chance of colonisation, then narrow them down to lakes within 5000m
attempt <- as.matrix(env_data_scaled) %*% as.matrix(calc_beta)
eta <- sweep(attempt, 2, calc_alpha, "+")
expeta <- exp(eta)
init_probabilities <- expeta/(1+expeta)

nn_all <- get.knnx(all_data[which(init_probabilities > 0.005),c("utm_x","utm_y")],all_data[c("utm_x","utm_y")],k=200)
nn_all$nn.index[nn_all$nn.dist > 5000] <- 0

periods <- list()

# These variables define actions in the loop
number_reps <- 5

# This defines in incremental increase over 50 years based on a 1 degree increase in temperature over 50 years. 10 represents a 1 degree increase.
temp_step <- 21/SDs_relData["temp"]

# These variables will be used to scale projections of numnber of populations nearby and distance to nearest population
mean_n_pop <- means_relData["n_pop"]
sd_n_pop <- SDs_relData["n_pop"]

sd_dist <- SDs_relData["pop_distL"]
mean_dist <- means_relData["pop_distL"]

introductions <- matrix(NA,nrow=nrow(env_data),ncol=1)
periods <- matrix(NA,nrow=nrow(env_data),ncol=n_loops)
# introduction_probs <- list()
# populations <- list()

print("Other stuff in place, starting forecast run.")

time1 <- Sys.time()
# data <- list()
for (s in 1:n_loops) {
  for (j in 1:number_reps) {
    # So, first step is to create predictions based purely off a rise in temperature. Because this is slightly different to what we'll do in the other steps, there needs to be an if clause.
    
    if (j==1) {
      # We create a table containing our new data
      
      tempIncrease <- rnorm(nrow(env_data_scaled), temp_step/5, 0.05)
      newData <- env_data_scaled
      newData[,"temp"] <- newData[,"temp"] + tempIncrease
      
      # Unfortunately because we're predicting for 600k + lakes, it's impossible to run the preditions all at one. So the loop below runs them in chunks of 10k lakes at a time.
      
      attempt <- as.matrix(newData) %*% as.matrix(calc_beta)
      eta <- greta::sweep(attempt, 2, calc_alpha, "+")
      expeta <- exp(eta)
      probabilities <- expeta/(1+expeta)
      
      # # # Now figure out which first upstream lakes pike have spread to using the gradient you defined before. I wrote a function for this which has been loaded above.
      # introductions_by_dispersal <- introduction_by_dispersal(full_data[full_data$Esox_lucius==1,]$waterBodyID,200, connect_db)
      # new_presences <- ifelse(full_data$Esox_lucius ==1, 1, ifelse(threshold_var > threshold_int, 1, ifelse(full_data$waterBodyID %in% introductions_by_dispersal,1, 0)))
      # 
      
      # Establish which lakes now have presences by asking whether or not we have an introduction, based on a bernoulli estimate.
      new_presences <- rbinom(length(probabilities), size = 1, prob=probabilities)
      all_presences <- ifelse(presences ==1, 1, ifelse(new_presences == 1, 1, 0))
      
      newData2 <- newData
    } else {
      
      ### And now we move on to the second part of the loop, whereby we redefine those variables pertaining to closest population and numbers of close populations at each step
      
      # Introduce the new temperature thing first, it's the simplest
      tempIncrease <- rnorm(nrow(newData2), temp_step/5, 0.05)
      newData2[,"temp"] <- newData2[,"temp"] + tempIncrease
      
      # Now figure out the new distance to closest population
      pop_proximity <- cbind(all_presences, all_data[,c("utm_x","utm_y","decimalLatitude","decimalLongitude","locationID")])
      
      data_presences <- as.data.frame(pop_proximity %>% filter(all_presences == 1)
                                      %>% dplyr::select(locationID, utm_x, utm_y, decimalLongitude, decimalLatitude))
      data_all <- pop_proximity %>%
        dplyr::select(utm_x,utm_y,locationID, decimalLongitude, decimalLatitude) %>%
        distinct() %>%
        as.data.frame()
      
      # The get.knnx function returns the distance from a lake in table B to the k closest lakes in table A
      
      nn <- get.knnx(data_presences[c("utm_x","utm_y")],data_all[c("utm_x","utm_y")],k=2)
      
      # We then figure out the distance to the closest lake by taking the first column
      dist_to_closest_pop <- ifelse(nn$nn.dist[,1]==0,nn$nn.dist[,2],nn$nn.dist[,1])
      locationID <- data_all$locationID
      distance_data <- as.data.frame(cbind(dist_to_closest_pop,locationID))
      
      ordered_distances <- merge(all_data["locationID"],distance_data,all.x=TRUE,by="locationID")
      
      new_distances <- log(as.numeric(as.character(ordered_distances$dist_to_closest_pop))+1)
      
      # So now that we have the new measurements for closest population, these need to be scaled against those that we had for the initial population. So we subtract the mean and divide by the standard deviation.
      
      newData2[,"pop_distL"] <- (new_distances-mean_dist)/sd_dist
      
      # Now we need to get the new measurements for number of close populations
      
      # And we take the number of lakes within a threshold of x metres by converting all distances less than x to 1, 
      # then simply summing the columns
      # pts_presence <- data_presences %>%
      #   transmute(latitude = decimalLatitude, longitude = decimalLongitude)
      # pts_presence <- st_as_sf(pts_presence, coords = c("longitude", "latitude"), 
      #                          crs = 4326)
      # pts_presence <-  st_transform(pts_presence, 32633)
      # 
      # pts_all <- all_data %>%
      #   transmute(latitude = decimalLatitude, longitude = decimalLongitude)
      # pts_all <- st_as_sf(pts_all, coords = c("longitude", "latitude"), 
      #                     crs = 4326)
      # pts_all <-  st_transform(pts_all, 32633)
      # pts_buf <- sf::st_buffer(pts_all, 5000)
      # 
      # int <- sf::st_intersects(pts_buf, pts_presence)
      # number_nearby_pop_vec <- lengths(int)
      
      
      rowNumbers_withPresence <- which(pop_proximity$all_presences==1)
      nn_inds <- nn_all$nn.index
      nn_inds[!(nn_inds %in% rowNumbers_withPresence)] <- 0
      nn_inds[nn_inds %in% rowNumbers_withPresence] <- 1
      number_nearby_pop_vec <- log(rowSums(nn_inds)+1)
      # 
      # Scale it like we did for the last variable
      
      newData2[,"n_pop"] <- (number_nearby_pop_vec-mean_n_pop)/sd_n_pop
      
      attempt <- as.matrix(newData2) %*% as.matrix(calc_beta)
      eta <- sweep(attempt, 2, calc_alpha, "+")
      expeta <- exp(eta)
      probabilities <- expeta/(1+expeta)
      
      # Establish which lakes now have presences by asking whether or not we have an introduction, based on a bernoulli estimate.
      new_presences <- rbinom(length(probabilities), size = 1, prob=probabilities)
      
      # # Introduce upstream lakes again
      # introductions_by_dispersal2 <- introduction_by_dispersal(data_presences$waterBodyID,200, connect_db)
      # new_presences <- ifelse(pop_proximity$new_presences ==1, 1, ifelse(threshold_var > threshold_int, 1, ifelse(pop_proximity$waterBodyID %in% introductions_by_dispersal2,1, 0)))
      
      
      # Establish which lakes now have presences.
      all_presences <- ifelse(pop_proximity$all_presences ==1, 1, ifelse(new_presences == 1, 1, 0))
    }
  }
  time2 <- Sys.time()
  difftime_1 <- round(as.numeric(difftime(time2, time1,
                                          units = "mins")),3)
  if (s %% 5 == 0) {print(paste0("Run ", s, " of ",  n_loops, " finished in ",difftime_1, " minutes. Estimated ", round(difftime_1*n_loops/s-difftime_1,5), " minutes left."))}
  periods[,s] <- all_presences
}

# period_introductions <- matrix(NA, nrow=nrow(all_data), ncol=ncol(periods))
# all_period_introductions <- list()

print("Finished forecasting.")

# ## Need to reorganise tables now
# for (j in 1:ncol(periods[[1]])) {
#   for (i in 1:length(periods)) {
#     period_introductions[,i] <- periods[[i]][,j]    
#   }
#   all_period_introductions[[j]] <- period_introductions
# }

forecasts <- periods
saveRDS(forecasts, file=paste0("./Data/",gsub(' ','_',species_name),"/forecasts.RDS"))
summary(forecasts)

