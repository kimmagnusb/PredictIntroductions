# Following two functions are found at https://github.com/joyofdata/joyofdata-articles/tree/master/roc-auc.

calculate_roc <- function(df, cost_of_fp, cost_of_fn, n=100) {
  tpr <- function(df, threshold) {
    sum(df$pred >= threshold & df$survived == 1) / sum(df$survived == 1)
  }
  
  fpr <- function(df, threshold) {
    sum(df$pred >= threshold & df$survived == 0) / sum(df$survived == 0)
  }
  
  cost <- function(df, threshold, cost_of_fp, cost_of_fn) {
    sum(df$pred >= threshold & df$survived == 0) * cost_of_fp + 
      sum(df$pred < threshold & df$survived == 1) * cost_of_fn
  }
  
  roc <- data.frame(threshold = seq(0,1,length.out=n), tpr=NA, fpr=NA)
  roc$tpr <- sapply(roc$threshold, function(th) tpr(df, th))
  roc$fpr <- sapply(roc$threshold, function(th) fpr(df, th))
  roc$cost <- sapply(roc$threshold, function(th) cost(df, th, cost_of_fp, cost_of_fn))
  
  return(roc)
  
}
plot_roc <- function(roc, threshold, cost_of_fp, cost_of_fn) {
  library(gridExtra)
  
  norm_vec <- function(v) (v - min(v))/diff(range(v))
  
  idx_threshold = which.min(abs(roc$threshold-threshold))
  
  col_ramp <- colorRampPalette(c("green","orange","red","black"))(100)
  col_by_cost <- col_ramp[ceiling(norm_vec(roc$cost)*99)+1]
  p_roc <- ggplot(roc, aes(fpr,tpr)) + 
    geom_line(color=rgb(0,0,1,alpha=0.3)) +
    geom_point(color=col_by_cost, size=4, alpha=0.5) +
    coord_fixed() +
    geom_line(aes(threshold,threshold), color=rgb(0,0,1,alpha=0.5)) +
    labs(title = sprintf("ROC")) + xlab("FPR") + ylab("TPR") +
    geom_hline(yintercept=roc[idx_threshold,"tpr"], alpha=0.5, linetype="dashed") +
    geom_vline(xintercept=roc[idx_threshold,"fpr"], alpha=0.5, linetype="dashed")
  
  p_cost <- ggplot(roc, aes(threshold, cost)) +
    geom_line(color=rgb(0,0,1,alpha=0.3)) +
    geom_point(color=col_by_cost, size=4, alpha=0.5) +
    labs(title = sprintf("cost function")) +
    geom_vline(xintercept=threshold, alpha=0.5, linetype="dashed")
  
  sub_title <- sprintf("threshold at %.2f - cost of FP = %d, cost of FN = %d", threshold, cost_of_fp, cost_of_fn)
  
  grid.arrange(p_roc, p_cost, ncol=2, sub=textGrob(sub_title, gp=gpar(cex=1), just="bottom"))
}
plot_pred_type_distribution <- function(df, threshold) {
  v <- rep(NA, nrow(df))
  v <- ifelse(df$pred >= threshold & df$survived == 1, "TP", v)
  v <- ifelse(df$pred >= threshold & df$survived == 0, "FP", v)
  v <- ifelse(df$pred < threshold & df$survived == 1, "FN", v)
  v <- ifelse(df$pred < threshold & df$survived == 0, "TN", v)
  
  df$pred_type <- v
  
  ggplot(data=df, aes(x=survived, y=pred)) + 
    geom_violin(fill=rgb(1,1,1,alpha=0.6), color=NA) + 
    geom_jitter(aes(color=pred_type), alpha=0.6) +
    geom_hline(yintercept=threshold, color="red", alpha=0.6) +
    scale_color_discrete(name = "type") +
    labs(title=sprintf("Threshold at %.2f", threshold))
}


# The rest I created for working with Greta.



## assoc ----
assoc <- function(temp) {
  lambda <- lambda_int + temp * lambda_coef
  cov2cor(lambda %*% t(lambda))
}

## Produces beta estimates, means and quantiles
get_param_ints <- function(parameter,draws,species_names,env_names,ssv=FALSE,quantile_n="mean") {
  n_species <- length(species_names)
  n_env_shared <- length(env_names)
  initial_matrix <- as.matrix(calculate(parameter,draws))
  if (quantile_n=="mean") { result_matrix <- apply(initial_matrix, 2, mean)
  } else {
    result_matrix <- apply(initial_matrix, 2, quantile,
                           probs = quantile_n)
  }
  if (ssv==TRUE) {
    ssv_matrix <- matrix(result_matrix,nrow = n_species, ncol=n_species)
    betas <- matrix(diag(ssv_matrix),nrow = n_species,ncol=1)
    colnames(betas) <- "ssv"
  } else {
    betas <- matrix(result_matrix,nrow = n_species, ncol=n_env_shared)
    colnames(betas) <- env_names
  }
  return(betas)
}


### Produces cooccurrence intervals
get_cooc_ints <- function(draws,species_names,quantile_n="mean") {
  n_species <- length(species_names)
  initial_matrix <- as.matrix(calculate(R,draws))
  if (quantile_n=="mean") { result_matrix <- apply(initial_matrix, 2, mean)
  } else {
    result_matrix <- apply(initial_matrix, 2, quantile,
                           probs = quantile_n)
  }
  cooc <- matrix(result_matrix,nrow = n_species, ncol=n_species)
  rownames(cooc) <- species_names
  colnames(cooc) <- species_names
  return(cooc)
}

# Creates lambdas with contrained priors
create_lambda <- function(n_species,n_latent) {
  lambda <- zeros(n_species, n_latent)
  diag(lambda) <- normal(0, 1, dim = n_latent, truncation = c(0,Inf))
  lower_idx <- lower.tri(lambda)
  lambda[lower_idx] <- normal(0, 1, dim = sum(lower_idx))
  return(lambda)
}

# Produce series of associations matrices based on incremental average temperature increases
pred_correlation <- function(draws,X_assoc_pred,species_names, R_pred) {
  n_species <- length(species_names)
  returned <- list()
  for (i in 1:length(X_assoc_pred)) {
    temp_R <- colMeans(as.matrix(calculate(R_pred[[i]],draws)))
    temp_R_matrix <- matrix(temp_R,n_species,n_species)
    colnames(temp_R_matrix) <- species_names
    rownames(temp_R_matrix) <- species_names
    diag(temp_R_matrix) <- 0
    returned[[i]] <- temp_R_matrix
  }
  return(returned)
}

### Deviance function to calculate deviance of given model
deviance_yp <- function (y, p) {
  -2 * sum(dbinom(y, 1, p, log = TRUE))
}


# Produces list of beta CIs and mean with descriptors
get_beta_list <- function(draws,beta_shared,beta_ssv,species_names,env_names,ssv=FALSE){
  demo_betas_upper <- as.data.frame(get_param_ints(beta_shared,draws,species_names,env_names,quantile_n = 0.975))
  demo_betas_lower <- as.data.frame(get_param_ints(beta_shared,draws,species_names,env_names,quantile_n = 0.025))
  demo_betas_mean <- as.data.frame(get_param_ints(beta_shared,draws,species_names,env_names))
  if (ssv==TRUE) {
    demo_betas_ssv_upper <- get_param_ints(beta_ssv,draws,species_names,env_names,ssv=TRUE,quantile_n = 0.975)
    demo_betas_upper <- as.data.frame(cbind(demo_betas_upper,demo_betas_ssv_upper))
    demo_betas_ssv_lower <- get_param_ints(beta_ssv,draws,species_names,env_names,ssv=TRUE,quantile_n = 0.025)
    demo_betas_lower <- as.data.frame(cbind(demo_betas_lower,demo_betas_ssv_lower))
    demo_betas_ssv_mean <- get_param_ints(beta_ssv,draws,species_names,env_names,ssv=TRUE)
    demo_betas_mean <- as.data.frame(cbind(demo_betas_mean,demo_betas_ssv_mean))
  }
  demo_betas_mean$species <- species_names 
  demo_betas_lower$species <- species_names 
  demo_betas_upper$species <- species_names 
  library(reshape2)
  upper_betas_reshaped <- melt(demo_betas_upper, id.vars = "species")
  lower_betas_reshaped <- melt(demo_betas_lower, id.vars = "species")
  mean_betas_reshaped <- melt(demo_betas_mean, id.vars = "species")
  full_betas1 <- merge(lower_betas_reshaped,mean_betas_reshaped,by=c("species","variable"))
  full_betas <- merge(full_betas1,upper_betas_reshaped,by=c("species","variable"))
  names(full_betas) <- c("species","variable","lower","mean","upper")
  return(full_betas)
}


# draws <- demo_draws_extra
# X_newdata <- X_full_new

# This predicts the likelihood of presence/absense of one species at one site, given the environmental values and the presence of other species.
predict_species_perc <- function(draws,X_newdata, site_id_newdata = NULL, X_assoc_newdata = NULL,species_list) {
  
  eta_new <- pred_eta(X_newdata, site_id_newdata,X_assoc_newdata)
  p_new <- ilogit(eta_new)
  
  fill_preds <- matrix(NA,nrow=nrow(p_new),ncol=ncol(p_new))
  
  for (i in 1:length(species_list)) {
    p_new1_draws1 <- calculate(p_new[, i], draws)
    p_new1_mn1 <- c(colMeans(as.matrix(p_new1_draws1)))
    fill_preds[,i] <- p_new1_mn1
  }
  colnames(fill_preds) <- species_list
  return(fill_preds)
}


# 
# 
# 
# eta_new <-

# betas_shared <- kauto_mean_betas
# site_env=t(kauto_env[1,])
# occupancy=kauto_occ[1,]
# cooccurrence=kauto_cooccurrence
# focal_species=4

# This predicts the likelihood of presence/absense of one species at one site, given the environmental values and the presence of other species.
#pred_cond(betas_shared=nofa_mean_betas,site_env=orp_env,occupancy=orp_species,cooccurrence=one_row_correlation_temp[[1]],focal_species=2)

pred_cond <- function(betas_shared,betas_ssv=NULL,site_env,occupancy,cooccurrence,focal_species,ssv=FALSE){
  library(tmvtnorm)
  if (ssv==TRUE) {
    betas_ssv_matrix <- matrix(0,length(betas_ssv),length(betas_ssv))
    diag(betas_ssv_matrix) <- betas_ssv
    betas <- cbind(betas_shared,betas_ssv_matrix)
  } else {betas <- betas_shared}
  mean1 <- c(betas %*% site_env)
  sigma1 <- cooccurrence # Need to edit this so it generates cooccurence matrix based on our temeprature
  diag(sigma1) <- 1
  
  lower1 <- c(ifelse(occupancy==0,-Inf,0))
  lower1[focal_species] <- -Inf
  upper1 <- c(ifelse(occupancy==0,0,Inf))
  
  lowerx1 <- c(lower1)
  lowerx1[focal_species] <- 0
  upperx1 <- c(upper1)
  upperx1[focal_species] <- Inf
  
  result <- ptmvnorm(mean=mean1,sigma=sigma1,lower=lower1,upper=upper1,lowerx=lowerx1,upperx=upperx1)
  return(result[1])
}


# This predicts the likelihood of presence/absense of one species at a series of sites, given the environmental values and the presence of other species.
pred_cond_entire <- function(draws,betas_shared,betas_ssv=NULL,site_matrix,occupancy_matrix,temp_vector,focal_species,ssv=FALSE,R_pred) {
  presences <- matrix(NA,ncol=1,nrow=nrow(site_matrix))
  for (i in 1:nrow(site_matrix)) {
    site_env <- matrix(site_matrix[i,],nrow=ncol(site_matrix),ncol=1)
    occupancy <- as.vector(occupancy_matrix[i,])
    cooccurrence <- pred_correlation(draws,temp_vector[i],colnames(occupancy_matrix),R_pred)
    diag(cooccurrence[[1]]) <- 1
    
    presences[i] <- pred_cond(betas_shared,betas_ssv,site_env,occupancy,cooccurrence[[1]],focal_species,ssv)
    if (i %% 10 == 0) {print(paste0(i," runs are complete.")) }
  }
  return(presences)
}
# conditional_entire_value <- pred_cond_entire(nofa_draws,nofa_mean_betas,site_matrix=as.matrix(NOFA_Data$X_val),occupancy_matrix=as.matrix(NOFA_Data$Y_val),temp_vector = as.matrix(NOFA_Data$X_val)[,2],focal_species=4)

# This takes a species correlation matrix and accentuates the values so it's easier to see
accentuate <- function(correlation_matrix_list) {
  unlisted <- unlist(correlation_matrix_list)
  
  times_factor <- 1/max(abs(unlisted))
  accentuated <- lapply(correlation_matrix_list, "*",times_factor)
  return(list(accentuated=accentuated,factor=times_factor))
}

# This narrows down our matrices to the species we want
narrowing <- function(correlation_matrix_list, focal_species_list) {
  narrowed_correlation_matrix_list <- vector("list", length(correlation_matrix_list))
  for(i in 1:length(correlation_matrix_list)) {
    narrowed_correlation_matrix_list[[i]] <- correlation_matrix_list[[i]][focal_species_list,focal_species_list]
  }
  return(narrowed_correlation_matrix_list)
}

# This creates a matrix with the mean associations given different temepratures gradients. Plans are to expand this to include credible intervals.

#association_matrix_list <-  demo_acc_narr
#focal_species <- "Golden_perch"

create_association_gradient <- function(association_matrix_list,focal_species) {
  ascmat_length <- length(association_matrix_list)
  n_species <- nrow(association_matrix_list[[1]])
  asc_map <- matrix(NA,n_species,ascmat_length)
  for(i in 1:ascmat_length) {
    asc_ind <- association_matrix_list[[i]][focal_species,]
    asc_map[,i] <- asc_ind
  }
  rownames(asc_map) <- colnames(association_matrix_list[[1]])
  return(asc_map)
}


# Easy function to plot associations given by "create_association_gradient".
#association_gradient_change_matrix <- asc_map
plot_change_asc <- function(association_gradient_change_matrix,legend_position=c(2,1)) {
  n_temps <- ncol(association_gradient_change_matrix)
  n_species <- nrow(association_gradient_change_matrix)
  plot(seq(-1,1,length.out=n_temps) ~ c(1:n_temps),type='n',xlab="Temperature",ylab="Association")
  cl <- brewer.pal(n_species, "Set1")
  for(i in 1:n_temps) {
    lines(seq(1,n_temps,1) ,association_gradient_change_matrix[i,],col=cl[i])
  }
  legend(legend_position[1], legend_position[2], legend=rownames(association_gradient_change_matrix), fill = cl)
}


# 
# library(greta)
# X_newdata <- X_full
#  X_assoc_newdata=X_assoc
# site_id_newdata <- as_data(Spatial)

#X_newdata = X_sim,X_assoc_newdata = X_assoc_sim,alpha = alpha_wc,beta=beta_wc,lambda_int = lambda_int_wc, lambda_coef = lambda_coef_wc,z=z_wc

pred_eta <- function(X_newdata, site_id_newdata = NULL, X_assoc_newdata = NULL,alpha=alpha,beta=beta,lambda_int=lambda_int, lambda_coef = lambda_coef,z=z, catchment_effect = NULL) {
  n_sites <- nrow(X_newdata)
  n_species <- nrow(lambda_coef)
  
  
  eta <- pred_env(X_newdata,alpha,beta)
  
  if (!is.null(site_id_newdata)) {
    eta <- eta + pred_site(site_id_newdata, catchment_effect, n_species)
  }
  
  if (!is.null(X_assoc_newdata)) {
    eta <- eta + pred_assoc(X_assoc_newdata,lambda_int,lambda_coef,z)
  }
  
  return(eta)
}

# pred()


pred_env <- function(X_newdata,alpha,beta) {
  XB <- X_newdata %*% beta
  eta <- sweep(XB, 2, alpha, "+")
  return(eta)
}


pred_site <- function(site_id_newdata, catchment_effect, n_species) {
  n_sites <- length(site_id_newdata)
  catchmentB <- sweep(ones(nrow(catchment_effect),n_species), 1, catchment_effect, "+")
  return(catchmentB)
}

#X_assoc <- X_assoc_sim

pred_assoc <- function(X_assoc,lambda_int,lambda_coef,z=z) {
  n_sites <- length(X_assoc)
  n_species <- nrow(lambda_coef)
  lambda_effect_rep <- kronecker(as_data(X_assoc), lambda_coef) # (site-lambdas stacked on top of one another)
  lambda_int_rep <- kronecker(ones(n_sites, 1), lambda_int)
  lambda_rep <- lambda_int_rep + lambda_effect_rep
  
  #Do the same thing for our zs so we can perform multiplication
  z_rep <- kronecker(ones(n_species, 1), t(z))
  
  t_e <- rowSums(lambda_rep * z_rep)
  dim(t_e) <- c(n_species, n_sites)
  e <- t(t_e)
  return(e)
}

# This SHOULD give us our new R matrix, but it needs work.
#X_assoc_pred <- X_assoc_pred[1,]
# lambda_int <- create_lambda(n_species,n_latent)

pred_one_R <- function(X_assoc_pred = NULL, lower_only = FALSE, lambda_int, lambda_coef){
  
  lambda <- lambda_int
  if (!is.null(X_assoc_pred)) {
    lambda <- lambda + lambda_coef * X_assoc_pred
  }
  
  R <- greta::cov2cor(lambda %*% t(lambda))
  
  if (lower_only) {
    lower_idx_R <- lower.tri(R)
    R <- R[lower_idx_R]
  }
  
  return(R)
}


# This SHOULD give us our new R matrix, but it needs work. new versioN

#X_assoc_pred <- X_assoc_pred_sim
#lambda_int=lambda_true_int


pred_R <- function(X_assoc_pred, lower_only = FALSE, lambda_int, lambda_coef) {
  
  if (!inherits(X_assoc_pred, "greta_array")) {
    X_assoc_pred <- as.matrix(X_assoc_pred)
  }
  
  n <- nrow(X_assoc_pred)
  
  result <- list()
  
  for (i in seq_len(n)) {
    result[[i]] <- pred_one_R(X_assoc_pred[i, ], lower_only = lower_only, lambda_int, lambda_coef)
  }
  
  return(result)
  
}

