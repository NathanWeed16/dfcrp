# source("SimData.R")
################################
# The Permutation Sampler
#' @param data The data which we're using to evaluate the permutation
#' @param family_col (string) The name of the column in data where family assignments are stored. Will default to "Family" if left unspecified
#' @param cluster_vec A vector of cluster assignments for each item (row) in the data
#' @param alpha The concentration parameter for the DFCRP. Will default to 1 if left unspecified
#' @param perm_mix The number of times we want the sampler to run per iteration of the overall Gibbs sampler
#' @param pct_flips The percent of the data that should be flipped each time a proposal permutation is generated
#' @param perm The initial permutation for the sampler. Will default to the identity permutation if left unspecified
perm_sampler<-function(family_vec, cluster_vec, alpha=1, perm_mix = 1, pct_flips = .2, perm = NULL){
  # Initialize quantities and vectors to be used later
  ncustomers <- length(cluster_vec)
  if(is.null(perm)) perm = c(1:ncustomers)
  permutation<-vector("list", (perm_mix+1))
  permutation[[1]]<-perm
  acc_prob<-0
  
  # Run the sampler
  for (m in 2:(perm_mix+1)){
    # Define a proposal density
    permutation_star=flip_more(permutation[[m-1]], max(2, round(pct_flips*ncustomers)))
    
    # Compute the log of the acceptance probability using the dfcrp pmf
    r1 = log_dfcrp_pmf_cond_cpp(family_vec, cluster_vec, alpha, permutation_star, start_val = 1)
    r2 = log_dfcrp_pmf_cond_cpp(family_vec, cluster_vec, alpha, permutation[[m-1]], start_val = 1)
    r = r1 - r2
    # Accept or reject the proposed permutation
    u<-log(runif(1))
    if(r>u){
      permutation[[m]]=permutation_star
      acc_prob=acc_prob+1
    }else{
      permutation[[m]]=permutation[[m-1]]
    }
    perm_prob=(acc_prob/perm_mix)
  }
  return(list(permutation[[perm_mix+1]], perm_prob))
}

################################
###---Proposal Density---###
#' @param vec The vector of items were going to rearrange or swap
#' @param n The number of items in vec we want to swap
flip_more<-function(vec, n){
  # First, ensure that the number of items to be flipped does not exceed the number of elements in the vector
  if (n>length(vec)){
    stop("n cannot be greater than the length of the vector")
  }
  
  # Randomly select which items in the vector are to be flipped 
  flips<-sample(seq_along(vec), n)
  
  # Randomly select, among the items to be flipped, what the new values of the items will be, and assign them to those new values
  permuted_vals <- sample(vec[flips])
  
  vec[flips] <- permuted_vals
  return(vec)
}


swap_equiv_clusters <- function(clustering, family_vec, perm) {
  
  # Unique cluster ids and families (stable order)
  cluster_ids <- sort(unique(clustering))
  unique_fams <- sort(unique(family_vec))
  
  # Build code for each cluster (string of 0/1 for each family)
  cluster_codes <- vapply(cluster_ids, function(cl) {
    items <- which(clustering == cl)
    fams_present <- unique(family_vec[items])
    bits <- as.integer(unique_fams %in% fams_present)
    paste(bits, collapse = "")
  }, FUN.VALUE = character(1))
  names(cluster_codes) <- as.character(cluster_ids)
  
  # Group clusters by their code
  code_groups <- split(cluster_ids, cluster_codes)
  
  new_perm <- perm
  
  for (grp in code_groups) {
    if (length(grp) <= 1) next  # nothing to swap
    
    # For each cluster in group: items belonging to cluster, and their positions in perm
    items_list <- lapply(grp, function(cl) which(clustering == cl))
    pos_list <- lapply(items_list, function(items) {
      # positions in perm where any of these items appear, ordered as in perm
      pos <- which(new_perm %in% items)
      # sort positions to preserve order as they appear in perm
      pos[order(match(new_perm[pos], items))] # ensure stable relative order (rarely needed)
    })
    
    lens <- vapply(pos_list, length, integer(1))
    if (length(unique(lens)) != 1) {
      next
    }
    
    # Extract blocks (items in perm order) and cyclically shift them
    blocks <- lapply(pos_list, function(p) new_perm[p])
    shifted <- blocks[c(2:length(blocks), 1)]
    
    # Write shifted blocks back into their original positions
    for (i in seq_along(pos_list)) {
      new_perm[pos_list[[i]]] <- shifted[[i]]
    }
  }
  
  return(new_perm)
}




