# Required packages:
library(lubridate)
library(tmvtnorm)
library(MCMCpack)
library(mvtnorm)

source("DFCRP_alpha&pmf.R")
source("Permutation_Sampler.R")

Rcpp::sourceCpp("Optimized.cpp")
################################
#' Produces samples of the cluster/table assignment through a collapsed Gibbs sampler for the dysfunctional family Chinese restaurant process (DFCRP).
#' 
#' @param data_to_cluster (data.frame) that includes the features to cluster on and the family specification. Could contain other, unused columns. 
#' @param family (string)  name of the column in the data_to_cluster data frame that specifies the family for each observation.
#' @param features (vector of strings) Vector of the p names of the column in the data_to_cluster data.frame that specifies the features on which to cluster.
#' @param location_vars (vector of strings) Vector of the subset of features columns that specify the X and Y coordinates of the observations
#' @param radius 
#' @param starting_assignment The initial table assignment to be used to start the Gibbs sampler. An initial seating assignment will be randomly determined if this input is left unspecified. This argument was intended to be used to restart a sampler that terminated prematurely. 
#' @param alpha The concentration parameter of the Chinese Restaurant Process (CRP).
#' @param perm_mix The number of times the permutation should be proposed per iteration of the sampler
#' @param pct_flips The proportion of the data that should be flipped for each proposed permutation
#' @param alpha_prop_precision The precision (inverse variance) of the gamma proposal density for the alpha sampler
#' @param a The shape hyperparameter for the Gamma prior on alpha
#' @param b The rate hyperparameter for the Gamma prior on alpha
#' @param mu0  Prior hyperparameter specification for the p-dimensional mean vector of the Normal Inverse Wishart. Will default to a vector of 0's if unspecified.
#' @param kmeans_k This value specifies the k in kmeans to form the clusters which will be used to estimate Sigma_0.  This value must be specified if S0 is NULL.
#' @param S0 Prior hyperparameter specification of Sigma_0 (prior mean of the covariance of the Inverse Wishart). Must be a pxp positive-definite matrix.
#' @param v0 Prior hyperparameter specification of the scalar degrees of freedom parameter of the Normal Inverse Wishart distribution. Must be greater than p+1. Will default to p+2 if unspecified. 
#' @param niter Number of iterations to run the sampler.
#' @param scale (logical) Indicates if the data should be scaled before sampled. Scaling will occur after defining neighborhoods. Defaults to TRUE.
#' @param m If m < niter, the function will write the sampler results to file after every m iterations.  Defaults to no writing.
#' @param output_filename The name of the data file that will be written as the function is running.  Defaults to the name 'dfcrp_sampler', which can be re-read into R using read.table(). Note that if running in parallel, you will want to change this name for each parallelization or files will be overwritten.  
#' @param print_status (logical) indicates whether a print status should be printed to Standard I/O. If TRUE, will print status every m iteration with a time stamp.
#' 
#' 
#' @return A matrix of size (niter + 1) by number of customers, where each row represents one sample from the Gibbs sampler. The initial values are the first row of the matrix.
#' 

dfcrp_sampler <- function(data_to_cluster, family = NULL, features, location_vars = NULL, radius = NULL, starting_assignment = NULL, alpha = 1, permutation = NULL, perm_mix=10, pct_flips=0.008, alpha_prop_precision = 20, a = 1, b = 0.01, mu0 = NULL, kmeans_k = NULL, S0 = NULL, v0 =NULL, niter = 2000, scale = TRUE, print_vec = NULL, output_filename = "dfcrp_sampler", print_status = F){
  
  ########## SET UP
  ##  Compute the variables to be used later
  # number of customers & features
  ncustomers <- nrow(data_to_cluster)
  nfeatures <- length(features)
  ## holder for table assignments--this will be returned 
  table_assign <- matrix(NA,nrow=niter+1,ncol=ncustomers)
  
  ## Feature Column number determination
  feature_cols <- which(names(data_to_cluster) %in% features)
  
  # If a radius and the location variables are provided, compute the neighborhood matrix 
  if (!is.null(radius)){
    neighbor_matrix<-distance_neighbor_maker(data_to_cluster = data_to_cluster, location_vars = location_vars, radius = radius)
  }
  
  # Scale the data, if necessary
  if(scale == T){data_to_cluster[, feature_cols] <- data.frame(scale(data_to_cluster[,feature_cols]))}
  
  # The column number of the family
  family_col <- which(names(data_to_cluster) %in% family)
  largest_family <- names(which.max(table(data_to_cluster[,family_col])))
  largest_family_idx <- which(data_to_cluster[,family_col] == largest_family)
  largest_family_size <- length(largest_family_idx)
  
  ########## INITIALIZE
  if (is.null(starting_assignment)) {
    X <- as.matrix(data_to_cluster[, feature_cols])
    cov_mat <- cov(X)
    cov_inv <- solve(cov_mat)
    
    # Convert family labels to integer IDs
    family_ids <- as.integer(as.factor(data_to_cluster[, family_col]))
    max_family <- max(family_ids)
    
    members_by_table <- vector("list", ncustomers)
    for (t in seq_len(ncustomers)) {
      members_by_table[[t]] <- integer(ncustomers)
    }
    
    # Logical matrix: rows = tables, cols = family IDs
    table_has_family <- matrix(FALSE, nrow = ncustomers, ncol = max_family)
    
    centroids <- matrix(0, nrow = ncustomers, ncol = ncol(X))
    sizes <- integer(ncustomers)
    
    # 1. initialize largest family tables
    for (idx in seq_along(largest_family_idx)) {
      j <- largest_family_idx[idx]
      table_assign[1,j] <- idx
      members_by_table[[idx]][1] <- j
      sizes[idx] <- 1
      centroids[idx, ] <- X[j, ]
      
      # mark family assigned
      table_has_family[idx, family_ids[j]] <- TRUE
    }
    
    next_table_label <- length(largest_family_idx) + 1
    rest_ix <- setdiff(seq_len(ncustomers), largest_family_idx)
    rest_ix <- sample(rest_ix)
    
    # 2. assign remaining customers
    for (j in rest_ix) {
      x_j <- X[j, ]
      fam_j <- family_ids[j]
      
      tables <- integer(0)
      d2s <- numeric(0)
      
      for (t in seq_len(next_table_label - 1)) {
        # skip incompatible families
        if (table_has_family[t, fam_j]) next
        
        diff <- x_j - centroids[t, ]
        d2 <- as.numeric(crossprod(diff, cov_inv %*% diff))
        
        tables <- c(tables, t)
        d2s <- c(d2s, d2)
      }
      
      # inverse-distance sampling
      weights <- 1 / d2s
      probs <- weights / sum(weights)
      
      if (length(tables) == 1) {
        best_t <- tables
      } else {
        best_t <- sample(tables, size = 1, prob = probs)
      }
      
      # assign customer to best table
      table_assign[1,j] <- best_t
      sizes[best_t] <- sizes[best_t] + 1
      members_by_table[[best_t]][sizes[best_t]] <- j
      
      # update centroid incrementally
      centroids[best_t, ] <- centroids[best_t, ] + (x_j - centroids[best_t, ]) / sizes[best_t]
      
      # mark family assigned
      table_has_family[best_t, fam_j] <- TRUE
    }
    rm(members_by_table)
    gc()
    
  } else {
    table_assign[1, ] <- starting_assignment
  }
  
  ## compute table specific quantities (means and covariances)
  tablesize <- numeric(0)
  mus       <- list()
  sigmas    <- list()
  for(i in 1:ncustomers){
    idx <- which(table_assign[1,] == i)
    mns <- table_means_covs(data_to_cluster[idx, feature_cols])
    mus[[i]] <- mns$mu
    sigmas[[i]] <- mns$sigma
    tablesize <- c(tablesize, length(idx))
  }
  
  ###########---HYPER-PARAMETERS---########
  ## Define parameters for the Normal Inverse Wishart prior-if not supplied by the user and not defaulted in the function
  if(is.null(mu0)) mu0 <- rep(0,nfeatures)
  if(is.null(v0)) v0 <- nfeatures + 2
  if(is.null(S0)){
    # The k in the k-means clustering algorithm will be the largest family size to estimate Sigma0
    kmeans_k  <- largest_family_size
    S0 <- estimate_Sigma0(data_to_cluster[,feature_cols], kmeans_k)
  }
  
  S <- cov(data_to_cluster[,feature_cols])
  k0 <- (det(S)/det(S0))^(-1/nfeatures)
  S0 <<- S0
  k0 <<- k0
  ## Keep track of where customers are currently sitting
  current_table_assignments <- table_assign[1,]
  
  # Initialize the permutation of the data
  if (is.null(permutation)){
    permutation<-sample(1:ncustomers, ncustomers, replace = F)
  }
  else{
    permutation<-permutation
  }
  
  # Optional: Initialize the acceptance probabilities of the Metropolis-Hastings samplers
  perm_prob=0
  alpha_acc=0
  alpha_vec<-numeric(niter)
  
  family_vec<-data_to_cluster[, family]
  # permutation_mixing<-vector("numeric", niter)
  ########## DO THE GIBBS
  for(iter in 1:niter){
    
    # First, sample the current permutation from the Metropolis sampler
    permutation_res<-perm_sampler(family_vec, cluster_vec=current_table_assignments, alpha=alpha, perm_mix=perm_mix, pct_flips = pct_flips, permutation)
    
    permutation<-permutation_res[[1]]
    perm <<- permutation
    perm_prob<-permutation_res[[2]]+perm_prob

    # Then sample the current value of alpha from the Metropolis-Hastings algorithm
    alpha_res<-alpha_dfcrp_sampler(family_vec, cluster_vec=current_table_assignments, alpha_int=alpha, alpha_prop_precision=alpha_prop_precision, permutation=permutation, a = a, b = b)
    alpha<-alpha_res[[1]]
    alpha_vec[iter]<-alpha
        if(iter>=(0.2*niter)){
        alpha_acc<-alpha_res[[2]]+alpha_acc
    }

    for(i in 1:ncustomers){
      # Identify which customer to work with according to the permutation, as well as the total number of tables
      j=permutation[i]
      n_tables<-max(current_table_assignments)
      ## current customers table
      current_table <- current_table_assignments[j]
      
      ## get the tables where customer's family is sitting
      same_family_tables <- NULL
      
      ## current customers family
      current_family <- data_to_cluster[j, family_col]
      same_family_tables <- unique(current_table_assignments[setdiff(which(data_to_cluster[,family_col]==current_family),j)])
      
      ########## REMOVE CUSTOMER J AND ADJUST 
      ## remove current customer from current table assignments
      current_table_assignments[j] <- NA
      
      ## decrease size of current table
      tablesize[current_table]    <- tablesize[current_table] - 1
      
      ## update current table size
      current_table_size <- tablesize[current_table]
      
      ## update sample mean and covariance for current table
      new_mns <- table_means_covs(data_to_cluster[which(current_table_assignments==current_table), feature_cols])
      mus[[current_table]] <- new_mns$mu
      sigmas[[current_table]] <- new_mns$sigma    
      
      # Assign a table to be the first empty table
      empty_tables<-which(tablesize == 0)
      if (length(empty_tables) > 0){
        first_empty <- min(empty_tables)
      }else{
        first_empty <- n_tables+1
      }
      
      
      ########## ELIMINATE POSSIBLE TABLES IF A NEIGHBORHOOD IS USED
      ## Get the neighbors of the current customer
      neighbor_tables <- which(tablesize!= 0)
      if(!is.null(radius)){
        current_neighbors <- which(neighbor_matrix[j,]==1)
        ## Get the tables of the neighbors 
        neighbor_tables   <- current_table_assignments[current_neighbors]
      }
      ########## COMPUTE LIKELIHOOD AND PRIOR
      ## Get the tables where the customer could sit
      occupied_tables <- which(tablesize != 0)
      no_family_tables <- setdiff(occupied_tables, same_family_tables)
      if(!is.null(radius)){
        temptables <- intersect(no_family_tables, neighbor_tables)
      }else{
        temptables <- no_family_tables
      }

      ## Compute the likelihood and prior probability for all temptables (after assigning j to each temptable)
      log_lik <- rep(-Inf, n_tables)
      log_dfcrp_prob <- rep(-Inf, n_tables)
      
      for(k in temptables){ 
      log_lik[k] <- likelihood_compute(data_to_cluster[j, feature_cols], tablesize[k], k0, mu0, v0, S0, mus[[k]], sigmas[[k]]) 
      current_table_assignments[j]<-k 
      relabeled <- relabel_cpp(current_table_assignments, permutation) 
      log_dfcrp_prob[k] <- log_dfcrp_pmf_cond_cpp(family_vec, cluster_vec = relabeled, alpha = alpha, permutation = permutation, start_val = i)
      }
      
      
      ## Compute the likelihood and prior probability for j being assigned to the empty table
      log_lik[first_empty]<-likelihood_compute(data_to_cluster[j, feature_cols], 0, k0, mu0, v0, S0, mus[[first_empty]], sigmas[[first_empty]])
      current_table_assignments[j]<-first_empty
      relabeled<-relabel_cpp(current_table_assignments, permutation)
      log_dfcrp_prob[first_empty] <- log_dfcrp_pmf_cond_cpp(family_vec, cluster_vec = relabeled, alpha = alpha, permutation = permutation, start_val = i)

      

      # Convert into standard probabilities to normalize
      log_lik<-log_lik-max(log_lik)
      lik<-exp(log_lik)
      
      log_dfcrp_prob<-log_dfcrp_prob-max(log_dfcrp_prob)
      dfcrp_prob<-exp(log_dfcrp_prob)
      
      prob_product<-lik*dfcrp_prob
      prob_product<-prob_product/sum(prob_product)
      ########## DRAW TABLE ASSIGNMENT
      newtable  <- sample(1:length(prob_product),1,prob=prob_product)
      
      ## store table assignment
      table_assign[iter+1,j] <- newtable
      
      ########## UPDATE THINGS WITH NEW ASSIGNMENT
      ## update current table assignments
      current_table_assignments[j] <- newtable
      
      ## update table size
      tablesize[newtable] <- tablesize[newtable] + 1
      ## update sample mean and covariance for current table
      new_mns <- table_means_covs(data_to_cluster[which(current_table_assignments==newtable), feature_cols])
      mus[[newtable]] <- new_mns$mu
      sigmas[[newtable]] <- new_mns$sigma  
    } #--end customer loop 
    
    table_assign[(iter+1), ]<-match(table_assign[(iter+1), ], unique(table_assign[(iter+1), ]))
    
  if (!is.null(print_vec)){
    if (iter %in% print_vec){
        cat("Partition:\n", file = output_filename, append = T)
        write(t(table_assign[iter+1,]), file = output_filename, ncolumns = ncustomers, append = T)
        cat("Permutation:\n", file = output_filename, append = T)
        write(permutation, file = output_filename, ncolumns = ncustomers, append = T)
        cat("Alpha:\n", file = output_filename, append = T)
        write(alpha, file = output_filename, ncolumns = ncustomers, append = T)
    }
  }
  
  
    if(print_status == T){
      if(iter %in% print_vec){ print(paste0("Iteration ", iter, " completed at ", now())) }}
    
  } #--end iterations loop 
  cluster_assignments<-table_assign
  permutation_prop<-perm_prob/niter
  # alpha_prop<-alpha_acc/niter
  return(cluster_assignments)
} #----dfcrp_sampler----#


################################
# HELPER FUNCTIONS:
################################
# Computes the likelihood from a multivariate t distribution after updating parameters based on the table size
likelihood_compute <- function(data_k, tablek_size, k0, mu0, v0, S0, mu_k, sigma_k){
  nfeatures <- length(mu0)
  if(tablek_size == 0){
    mun <- mu0
    vn <- v0
    kn <- k0
    Sn <- S0
  } else if (tablek_size == 1) {
    # Only one point → use its mean as mu_k, but for covariance, use prior
    kn  <- k0 + 1
    vn  <- v0 + 1
    mun <- (k0*mu0 + mu_k)/(k0 + 1)  # mu_k is the single point
    Sn  <- S0  # use prior covariance
  } else {
    kn  <- k0 + tablek_size
    mun <- (k0*mu0 + tablek_size*mu_k)/(k0 + tablek_size)
    vn  <- v0 + tablek_size
    Sn  <- S0 + sigma_k*(tablek_size-1) + (k0*tablek_size)/(k0+tablek_size) * (mu_k - mu0)%*%t(mu_k-mu0)
  }
  Sigman <- (Sn*(kn+1))/(kn*(vn-nfeatures+1))
  dfn    <- vn-nfeatures+1
  return(dmvt(data_k,delta=mun,sigma=Sigman,df=dfn,log=T))
  
}
#----likelihood_compute-----#

################################
# Updates table means and covariances after removing customer j
table_means_covs <- function(table_data){
  table_size <- nrow(table_data)
  nfeatures <- ncol(table_data)
  
  if(table_size<=1){
    sigma <- matrix(0, nrow=nfeatures, ncol=nfeatures)
  }else{
    sigma  <- cov(table_data)
  }
  
  if(table_size==0){
    mu <- rep(0, times = nfeatures)
  }else{
    mu <- colMeans(table_data)
  }
  
  return(list(mu = mu, sigma = sigma))
} #-----table_means_covs-------#

################################
# Returns an empirical covariance matrix estimate after running a k-means algorithm 
estimate_Sigma0 <- function(data_to_analyze, kmeans_k){
  k_means_results <- kmeans(data_to_analyze, kmeans_k)
  data_to_analyze$k_means_cluster <- k_means_results$cluster
  cov.list <- lapply(sort(unique(data_to_analyze$k_means_cluster)), function(x){
    subdata <- data_to_analyze[data_to_analyze$k_means_cluster==x, -ncol(data_to_analyze)]
    if(nrow(subdata) >= 2){
      cov(subdata)
    } else {
      diag(ncol(subdata)) * 1e-6   # tiny PD matrix for singleton cluster
    }
  })
  
  
  return(apply(simplify2array(cov.list), 1:2, mean, na.rm = T))
}  #-----estimate_Sigma0-------#

################################
# Produces a neighborhood matrix of 1s and 0s based on the location variables of the data and the desired radius
distance_neighbor_maker <- function(data_to_cluster, location_vars, radius){
  location_cols = which(names(data_to_cluster) %in% location_vars)
  ncustomers = nrow(data_to_cluster)
  
  ## 1 = neighbor; 0 = not neighbor
  distance_matrix <- as.matrix(dist(data_to_cluster[,location_cols]))
  distance_matrix_idx <- matrix(0, nrow = ncustomers, ncol = ncustomers)
  
  # If radius is 0-force all customers to only be neighbors by themselves--even if two points have distance 0
  if(radius == 0){
    diag(distance_matrix_idx) <- 1
  } else {
    distance_matrix_idx[distance_matrix <= radius] <- 1
    diag(distance_matrix_idx) <- 0
  }
  
  return(distance_matrix_idx)
} #------distance_neighbor_maker--------#

