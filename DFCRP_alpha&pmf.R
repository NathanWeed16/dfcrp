library(compiler)
library(Rcpp)
####### DFCRP PMF ########
#' @param family_vec A vector of family assignments. 
#' @param cluster_vec A vector of current cluster assignments.
#' @param alpha The concentration parameter for the DFCRP. Will default to 1 if not specified.
#' @param permutation The permutation of the items. Will default to the identity permutation if not specified.
log_dfcrp_pmf<-function(family_vec, cluster_vec, alpha=1, permutation=c(1:nrow(data))){
  #### SET UP ####
  # Initialize quantities and vectors to be used later
  family_same_mat <- outer(family_vec, family_vec, FUN = "==")
  ncustomers <- length(cluster_vec)
  log_joint=0
  clusters<-list()
  clusters[[1]]<-permutation[1]
  sizes<-integer(0)
  sizes[1]<-1
  log_prob_vec <- numeric(length(sizes) + 1)
  
  # Relabel the cluster labels in canonical order according to the permutation
  cluster_vec<-relabel_cpp(cluster_vec, permutation)
  
  # Compute the prior
  for (i in 2:length(permutation)){
    
    # Identify which item (row in the data) comes next in the permutation and the family assignment of that item
    j<-permutation[i]
    xj<-family_vec[j]
    
    # Reset the log_prob_vec for each itration
    log_prob_vec[] <- -Inf
    
    # Identify which non-empty clusters don't already have members of item xj's family present
    valid <- logical(length(clusters))
    row_j <- family_same_mat[j, ]
    
    for (i in seq_along(clusters)) {
      valid[i] <- !any(row_j[clusters[[i]]])
    }
    
    # Assign the probability vector a value based on the size of the cluster (if it's a valid cluster)
    for (k in seq_along(clusters)){
      if (valid[k]){
        log_prob_vec[k]<-log(sizes[k])
      }
    }
    
    # Assign the unnormalized probability of j being assigned to a new cluster to be alpha
    log_prob_vec[length(clusters)+1] <- log(alpha)
    
    # Normalize the probability vector
    m <- max(log_prob_vec)
    log_denom <- m + log(sum(exp(log_prob_vec - m)))
    log_prob_vec=log_prob_vec-log_denom
    
    # Identify the actual cluster assignment of item xj and add the log probability to our current log_joint
    k_actual <- cluster_vec[j]
    log_joint=log_joint+log_prob_vec[k_actual]
    
    # Actualize the list of clusters and their sizes
    if (k_actual <= length(clusters)) {
      clusters[[k_actual]] <- c(clusters[[k_actual]], j)
      sizes[k_actual]<-sizes[k_actual]+1
    } else {
      clusters[[k_actual]] <- j
      sizes[k_actual]<-1
      length(log_prob_vec)<-length(log_prob_vec)+1
    }
  }
  
  # Return the joint log probability
  return(log_joint)
}



####### DFCRP Conditional PMF ########
#' @param family_vec A vector of family assignments. 
#' @param cluster_vec (vector) A vector of current cluster assignments.
#' @param alpha The concentration parameter for the DFCRP. Will default to 1 if not specified.
#' @param permutation The permutation of the items. Will default to the identity permutation if not specified.
#' @param start_val The point where the function should begin to compute probabilities. The function will treat all values before start_val as given.
log_dfcrp_pmf_cond<-function(family_vec, cluster_vec, alpha=1, permutation=c(1:length(cluster_vec)), start_val){
  #### SET UP ####
  # Initialize quantities and vectors to be used later
  ncustomers <- length(cluster_vec)
  log_joint=0
  sizes<-integer(length(unique(cluster_vec))+1)
  log_alpha<-log(alpha)
  # Relabel the cluster labels in canonical order according to the permutation
  cluster_vec<-relabel_cpp(cluster_vec, permutation)
  # Assign all the items before start_val to the partition
  family_in_clusters<-matrix(TRUE, nrow=max(family_vec), ncol=length(sizes))
  if (start_val==1){
    sizes[1]<-1
    start_val=2
    family_in_clusters[family_vec[permutation[1]], cluster_vec[permutation[1]]]<-F
  }
  else{
    for (i in 1:(start_val-1)){
      j<-permutation[i]
      k_actual<-cluster_vec[j]
      sizes[k_actual]<-sizes[k_actual]+1
      family_in_clusters[family_vec[j], k_actual]<-F
    }
  }
  log_prob_vec <- numeric(length(sizes))
  log_sizes<-log(sizes)
  first_empty<-min(which(sizes==0))
  # Compute the log probability of the partition beginning from start_val
  for (i in start_val:ncustomers){
    
    # Identify which item (row in the data) comes next in the permutation and the family assignment of that item
    j<-permutation[i]
    xj<-family_vec[j]
    
    # Reset the log_prob_vec for each iteration
    log_prob_vec[] <- -Inf
    
    # Identify which non-empty clusters don't already have members of item xj's family present
    valid <- family_in_clusters[xj, ]
    
    # Assign the probability vector a value based on the size of the cluster (if it's a valid cluster)
    
    log_prob_vec[valid]<-log_sizes[valid]
    
    # Assign the unnormalized probability of j being assigned to a new cluster to be alpha
    log_prob_vec[first_empty] <- log_alpha
    
    # Normalize the probability vector
    m <- max(log_prob_vec)
    log_denom <- m + log(sum(exp(log_prob_vec - m)))
    log_prob_vec=log_prob_vec-log_denom
    
    # Identify the actual cluster assignment of item xj and add the log probability to our current log_joint
    k_actual <- cluster_vec[j]
    log_joint=log_joint+log_prob_vec[k_actual]
    
    # Actualize the list of clusters and their sizes
    sizes[k_actual]<-sizes[k_actual]+1
    log_sizes[k_actual]<-log(sizes[k_actual])
    family_in_clusters[xj, k_actual]<-F
    if(k_actual == first_empty) {first_empty <- min(which(sizes==0))}
  }
  
  # Return the joint log probability
  return(log_joint)
}

# Demonstrates full and conditional PMF equivalency as well as the time speedup
# test_vec<-integer(100)
# test_vec2<-integer(100)
# 
# for (i in 1:100){
#   simdat<-simulate_family_data(K=140,  X_limits = c(0,100), Y_limits = c(0,100), a_d = 5,b_d = 2,nu = 10,Sigma = 2*diag(3),J=6,p = t(c(0.3, 0.12, 0.08, 0.08, 0.12, 0.3)))
#   family_vec<-simdat$Family
#   family_same_mat <- outer(family_vec, family_vec, FUN = "==")
#   test_vec[i]<-log_dfcrp_pmf_cond_cpp(family_same_mat, cluster_vec=c(1:nrow(simdat)), permutation = c(1:nrow(simdat)), alpha=1, start_val = 1)
#   test_vec2[i]<-log_dfcrp_pmf(simdat, "Family", c(1:nrow(simdat)))
# }
# all.equal(test_vec, test_vec2)
# 
# Rprof("Test", interval=0.001)
# log_dfcrp_pmf(simdat, "Family", c(1:nrow(simdat)))
# Rprof(NULL)
# summaryRprof("Test")
# 
# Rprof("Test", interval=0.001)
# log_dfcrp_pmf_cond_cpp(family_same_mat, cluster_vec=c(1:nrow(simdat)), permutation = c(1:nrow(simdat)), alpha=1, start_val = 1)
# Rprof(NULL)
# summaryRprof("Test")



################################
# The Alpha Sampler
#' @param family_vec A vector of family assignments for each item (row) in the data.
#' @param cluster_vec A vector of cluster assignments for each item (row) in the data.
#' @param alpha_int The initial value of alpha for the sampler. Will default to 1 if left unspecified.
#' @param alpha_prop_precision The precision of the lognormal alpha proposal density.
#' @param a The shape hyperparameter for the gamma prior on alpha.
#' @param b The rate hyperparameter for the gamma prior on alpha.
#' @param permutation The permutation to be used in the sampler. If left without an argument, will default to the identity permutation.
alpha_dfcrp_sampler<-function(family_vec, cluster_vec, alpha_int=1, alpha_prop_precision, a, b, permutation = NULL){
  # Define the number of customers and initialize the acceptance indicator
  ncustomers <- length(cluster_vec)
  if(is.null(permutation)) permutation <- seq_len(ncustomers)
  acc<-NULL
  
  # Generate a candidate value, alpha_star
  alpha_star<-rlnorm(1, log(alpha_int) - 1/(2*alpha_prop_precision), sqrt(1/alpha_prop_precision))
  # Compute terms needed in the acceptance probability (the prior on both the current and proposed values, as well as the proposal for the current and proposed values)
  star_prior<-dgamma(alpha_star, a, b, log=T)
  old_prior<-dgamma(alpha_int, a, b, log=T)
  star_proposal<-dlnorm(alpha_star, log(alpha_int) - 1/(2*alpha_prop_precision), sqrt(1/alpha_prop_precision), log = T)
  old_proposal<-dlnorm(alpha_int, log(alpha_star) - 1/(2*alpha_prop_precision), sqrt(1/alpha_prop_precision), log = T)

  # Compute the acceptance probability
  r=log_dfcrp_pmf_cond_cpp(family_vec, cluster_vec, alpha_star, permutation, start_val = 1)-log_dfcrp_pmf_cond_cpp(family_vec, cluster_vec, alpha_int, permutation, start_val = 1)+star_prior-old_prior+old_proposal-star_proposal
  # Accept or reject the proposed value for alpha
  if (r>log(runif(1))){
    alpha_final=alpha_star
    acc=T
  }else{
    alpha_final=alpha_int
    acc=F
  }
  return(list(alpha_final, acc))
}



clustmean <- function(alpha, ncust){
  alpha*(digamma(ncust + alpha)-digamma(alpha))
}
