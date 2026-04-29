################################
# The Permutation Sampler
#' @param family_vec A vector of family membership for each item (row) in the data.
#' @param cluster_vec A vector of cluster assignments for each item (row) in the data.
#' @param alpha The concentration parameter for the DFCRP. Will default to 1 if left unspecified.
#' @param perm The initial permutation for the sampler. Will default to the identity permutation if left unspecified.
perm_sampler <- function(family_vec, cluster_vec, alpha=1, perm=c(1:ncustomers)){
  
  # Initialize quantities and vectors to be used later
  ncustomers <- length(cluster_vec)
  acc_prob<-0

  # To generate a proposal permutation, we swap the last item with a random item from the permutation
  # This proposal is not only symmetric, swapping just this last element enables us to use the exchangable CRP prior
  permutation_star <- perm
  i <- sample(1:length(permutation_star), 1)
  tmp <- permutation_star[i]
  permutation_star[i] <- permutation_star[length(permutation_star)]
  permutation_star[length(permutation_star)] <- tmp

  # Compute the log of the acceptance probability using the dfcrp pmf
  r1=log_dfcrp_pmf_cond_cpp(family_vec, cluster_vec, alpha, permutation_star, start_val = 1)
  r2=log_dfcrp_pmf_cond_cpp(family_vec, cluster_vec, alpha, perm, start_val = 1)
  r=r1-r2

  # Accept or reject the proposed permutation based on the acceptance ratio
  u<-log(runif(1))
  if(r>u){
    permutation=permutation_star
    acc_prob=acc_prob+1
  }else{
    permutation=perm
  }

  return(list(permutation, acc_prob))
}



