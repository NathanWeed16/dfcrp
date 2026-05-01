################################
# The Alpha Sampler
#' @param data The data which we're using to evaluate alpha
#' @param cluster_vec A vector of cluster assignments for each item (row) in the data
#' @param alpha_int The initial value of alpha for the sampler. Will default to 1 if left unspecified
#' @param alpha_prop_precision The tuning parameter to control the variance of the gamma proposal density
#' @param a The shape hyperparameter for the prior on alpha
#' @param b The rate hyperparameter for the prior on alpha
alpha_crp_sampler<-function(data, cluster_vec, alpha_int=1, alpha_prop_precision, alpha_prior_shape, alpha_prior_rate){
  # Define the number of customers and initialize the acceptance indicator
  acc<-NULL
  
  # Generate a candidate value, alpha_star
  alpha_star<-rgamma(1, alpha_prop_precision*(alpha_int)**2, alpha_prop_precision*alpha_int)
  
  # Compute terms needed in the acceptance probability (the prior on both the current and proposed values, as well as the proposal for the current and proposed values)
  star_prior<-dgamma(alpha_star, alpha_prior_shape, alpha_prior_shape, log=T)
  old_prior<-dgamma(alpha_int, alpha_prior_shape, alpha_prior_shape, log=T)
  star_proposal<-dgamma(alpha_star, alpha_prop_precision*(alpha_int)**2, alpha_prop_precision*alpha_int, log=T)
  old_proposal<-dgamma(alpha_int, alpha_prop_precision*(alpha_star)**2, alpha_prop_precision*alpha_star, log=T)
  
  # Compute the acceptance probability
  r=log_crp_pmf(data, cluster_vec, alpha=alpha_star)-log_crp_pmf(data, cluster_vec, alpha=alpha_int)+star_prior-old_prior+old_proposal-star_proposal
  
  # Accept or reject the proposed value for alpha
  if (r>log(runif(1))){
    alpha_final=alpha_star
    acc=T
  }else{
    alpha_final=alpha_int
    acc=F
  }
  return(alpha_final)
}

# The CRP PMF
#' @param data The data for which we're going to compute the probability of a specific cluster assignment (we only use the family and cluster assignment columns)
#' @param cluster_vec A vector with the true cluster assignments of the data
#' @param alpha The concentration parameter for the CRP. Will default to 1 if not specified
log_crp_pmf<-function(data, cluster_vec, alpha=1){
  # SET UP
  # Define the number of customers and initialize the log-joint probability (this will be returned)
  ncustomers <- length(cluster_vec)
  log_joint=0
  cluster_vec <- canonical(cluster_vec)
  
  # Keep track of the cluster
  clusters<-list()
  clusters[[1]]<-1
  sizes<-sapply(clusters, length)
  # Compute the likelihood
  for (i in 2:ncustomers){
    log_prob_vec <- rep(-Inf, length=length(sizes)+1)
    # Assign the probability vector a value based on the size of the cluster if it's a valid cluster for the item xj to be assigned
    for (k in seq_along(clusters)){
      if(!is.na(sizes[k])){
        log_prob_vec[k]<-log(sizes[k])
      }else{
        log_prob_vec[k]<--Inf
      }
    }
    
    # Assign the unnormalized probability of xj being assigned to a new cluster to be alpha
    log_prob_vec[length(clusters)+1] <- log(alpha)
    # Normalize the probability vector
    m <- max(log_prob_vec)
    log_denom <- m + log(sum(exp(log_prob_vec - m)))
    log_prob_vec=log_prob_vec-log_denom
    
    # Identify the actual cluster assignment of item xj and add the log probability to our current log_joint
    k_actual <- cluster_vec[i]
    log_joint=log_joint+log_prob_vec[k_actual]
    
    # Actualize the list of clusters and their sizes
    if (k_actual <= length(clusters)) {
      clusters[[k_actual]] <- c(clusters[[k_actual]], i)
      sizes[k_actual]<-sizes[k_actual]+1
    } else {
      clusters[[k_actual]] <- i
      sizes[k_actual]<-1
    }
  }
  
  # Return the accumulated log probability
  return(log_joint)
}

# Relabels cluster labels to be in canonical order
canonical <- function(partition) {
  cmatrix <- outer(partition,partition,FUN="==")
  labels <- unique(cmatrix)
  n <- length(labels[,1])
  newpart <- partition
  for (i in 1:n) {
    newpart[which(labels[i,])] <- i-1
  }
  newpart+1
}
