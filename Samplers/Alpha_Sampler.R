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
