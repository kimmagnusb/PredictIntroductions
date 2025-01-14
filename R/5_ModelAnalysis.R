

# This will actually be a script for checking either model

if (use_buffered_model == TRUE) {
  model_output <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/buffered_model_output.RDS"))
  model_data_extra <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/buffered_model_data.RDS"))
} else {
  model_output <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/whole_model_output.RDS"))
  model_data_extra <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/whole_model_data.RDS"))
}



draws <- model_output$draws
beta <- model_output$beta
p <- model_output$p
alpha <- model_output$alpha


### At this point let's run some quick convergence diagnostics. Values for gelman 
#   diagnostics should be close to 1, for the effective size the higher the better.
gelmans <- coda::gelman.diag(calculate(beta,draws),multivariate = FALSE)
effective_sizes <- summary(coda::effectiveSize(calculate(beta,draws)))


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
env_data <- all_data_usable %>%
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

# Get full index of nearby lakes that have a chance of colonisation, then narrow them down to lakes within 5000m
attempt <- as.matrix(env_data_scaled) %*% as.matrix(calc_beta)
eta <- sweep(attempt, 2, calc_alpha, "+")
expeta <- exp(eta)
init_probabilities <- expeta/(1+expeta)


p_dev <- init_probabilities[,1]
Y_dev <- all_data_usable[,"introduced"]
model_deviance <- deviance_yp(Y_dev,p_dev)

print("Calculated deviance. Now going through getting uncertainty.")


# Now we scale our full data based on these means
env_data <- all_data$raw_data %>%
  transmute(areaL = log(area_km2+1),
            dist_roadL = log(distance_to_road+1),
            temp = eurolst_bio10,
            pop_distL = log(pop_dist+1),
            SCI = log(((perimeter_m/1000)/(2*sqrt(pi*area_km2)))+1),
            HFP = HFP,
            n_pop = log(nearby_pops+1))

env_data_takeMeans2 <- t(apply(env_data, 1, "-", means_relData))
env_data_scaled2 <- t(apply(env_data_takeMeans2, 1, "/", SDs_relData))

# Create the matrix to stuff all the values into
intVal_mat <- data.frame(lower=numeric(),
                         mean=numeric(), 
                         upper=numeric(), 
                         stringsAsFactors=FALSE) 

if (nrow(env_data_scaled2) < 10000) {
  int_eta <- pred_env(env_data_scaled2,alpha,beta)
  int_val <- ilogit(int_eta)
  int_draws <- calculate(int_val,draws)
  
  # The following creates our probabilities of introduction.
  int_pred_ints <- as.data.frame(t(apply(as.matrix(int_draws) , 2 , quantile , probs = c(0.025,0.5,0.975) , na.rm = TRUE )))
  intVal_mat <- rbind(intVal_mat,int_pred_ints)
} else {
  loops <- ceiling(nrow(env_data_scaled2)/10000)
  
  
  time1 <- Sys.time()
  for (i in 1:(loops)) {
    if(i != (loops)) {
      analyse_df_chunk <- env_data_scaled2[(1+(i-1)*10000):(i*10000),]
    } else {
      analyse_df_chunk <- env_data_scaled2[(1+(i-1)*10000):nrow(env_data_scaled2),]
    }
    int_eta <- pred_env(analyse_df_chunk,alpha,beta)
    int_val <- ilogit(int_eta)
    int_draws <- calculate(int_val,draws)
    
    # The following creates our probabilities of introduction.
    int_pred_ints <- as.data.frame(t(apply(as.matrix(int_draws) , 2 , quantile , probs = c(0.025,0.5,0.975) , na.rm = TRUE )))
    intVal_mat <- rbind(intVal_mat,int_pred_ints)
    
    # Quick function to let us know how the loop is progressing
    time2 <- Sys.time()
    difftime_1 <- round(as.numeric(difftime(time2, time1,
                                            units = "mins")),4)
    if (i %% 10 == 0) {print(paste0("Run ", i, " finished in ",difftime_1, " minutes. Estimated ", round(difftime_1*(loops+1)/i-difftime_1,10), " minutes left.") )}
  }
  
}


## Now let's check the beta estimates
beta_ints <- get_beta_list(draws,beta_shared=beta,species_names="introduced",
                           env_names=colnames(env_data))


colnames(intVal_mat) <- c("lower","mean","upper")

intVal_mat$width <- with(intVal_mat,(upper-lower)/2)

model_analysis <- list(intervals = intVal_mat, deviance = model_deviance, 
                       conv_diags = list(gelmans, effective_sizes), parameter_estimates = beta_ints)

if (use_buffered_model == TRUE) {
  saveRDS(model_analysis,paste0("./Data/",gsub(' ','_',species_name),"/buffered_model_analytics.RDS"))
} else {
  saveRDS(model_analysis,paste0("./Data/",gsub(' ','_',species_name),"/whole_model_analytics.RDS"))
}


