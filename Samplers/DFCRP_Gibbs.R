# Required packages:
library(lubridate)
library(mvtnorm)
library(Matrix)

source("Samplers/Alpha_Sampler.R")
source("Samplers/Perm_Sampler.R")

Rcpp::sourceCpp("Samplers/DFCRP_Prior_Optimized.cpp")
################################
#' Produces samples of the cluster/table assignment through a collapsed Gibbs sampler for the dysfunctional family Chinese restaurant process (DFCRP).
#' @param data_to_cluster (data.frame) that includes the features to cluster on and the family specification. Could contain other, unused columns. 
#' @param family (string)  name of the column in the data_to_cluster data frame that specifies the family for each observation.
#' @param features (vector of strings) Vector of the p names of the column in the data_to_cluster data.frame that specifies the features on which to cluster.
#' @param location_vars (vector of strings) Vector of the subset of features columns that specify the X and Y coordinates of the observations
#' @param radius Value of roe to be used in the neighborhood modification. Any cluster with no craters within this distance of the crater being reassigned is given a likelihood of 0.
#' @param starting_assignment The initial table assignment to be used to start the Gibbs sampler. The algorithm will start each observation in its own cluster if left unspecified. This argument was intented to be used to restart a sampler that terminated prematurely. 
#' @param alpha The concentration parameter of the Chinese Restaurant Process (CRP).
#' @param alpha_prop_precision The precision of the gamma proposal density for the alpha sampler.
#' @param alpha_prior_shape The shape hyperparameter for the gamma prior on alpha.
#' @param alpha_prior_rate The rate hyperparameter for the gamma prior on alpha.
#' @param permutation The starting value for the permutation. Will randomly select a permutation if left unspecified.
#' @param niter Number of iterations to run the sampler. For full scans, multiply the number of desired full scans by the number of craters to assign.
#' @param print_vec Vector of iterations for which the results will be outputted to a CSV. If initial values are desired, include the value 0. Defaults to no output.
#' @param output_filename The name of the data file that will be written as the function is running.  Defaults to the name 'dfcrp_sampler', which can be re-read into R using read.csv(). Note that if running in parallel, you will want to change this name for each parallelization or files will be overwritten.  
#' @param print_status (logical) indicates whether a print status should be printed to Standard I/O. If TRUE, will print status for every iteration in print_vec with a time stamp.
#' @param mu0 The multivariate normal mean parameter for the prior on mu, the cluster-specific mean vector. This should be specific to the size of the dataset and the application.
#' @param sigma0 The multivariate normal covariance matrix parameter for the prior on mu, the cluster-specific mean vector. This should be specific to the size of the dataset and the application.
#' @param sigma_prop_mat The covariance matrix for the multivariate normal proposal distribution on Sigma, the cluster-specific covariance matrix.
#' @param init_mus This parameter allows the user to specify the initial values of the cluster-specific mean vectors (mu). This argument was intented to be used to restart a sampler that terminated prematurely. 
#' @param init_sigmas This parameter allows the user to specify the initial values of the cluster-specific covariance matrices (Sigma). This argument was intented to be used to restart a sampler that terminated prematurely. 
#' 
#' @return A vector of length (number of customers), where each value is the cluster assignment for that item as of the last iteration, ordered canonically. For more draws, use the print_vec argument.

dfcrp_sampler <- function(data_to_cluster, family = NULL, features, location_vars = NULL, radius = NULL, starting_assignment = NULL, alpha = 1, alpha_prop_precision = 20, alpha_prior_shape = 1, alpha_prior_rate = 0.01, permutation = NULL, niter, print_vec = NULL, output_filename = "dfcrp_sampler", print_status = F, mu0 = c(2068.171, -1105.497, 3.8), sigma0 = matrix(c(920^2, 0, 0, 0, 600^2, 0, 0, 0, 0.65^2), nrow = 3), sigma_prop_mat = matrix(c(3000, 0, 0, 0, 0.9, 0, 0, 0, 0.2), nrow = 3), init_mus = NULL, init_sigmas = NULL){
  
  ########## SET UP
  ##  Compute the variables to be used later
  ncustomers <- nrow(data_to_cluster)
  nfeatures <- length(features)
  
  # Stores the current partition, will be returned iteration by iteration in print_vec
  current_table_assignments <- rep(NA_integer_,length=ncustomers)
  
  ## Column number determination
  feature_cols <- which(names(data_to_cluster) %in% features)
  diameter_col <- which(names(data_to_cluster) %in% c("Diameter"))
  family_col <- which(names(data_to_cluster) %in% family)
  
  # If a radius and the location variables are provided, compute the neighborhood matrix 
  if (!is.null(radius)){
    neighbor_matrix<-distance_neighbor_maker(data_to_cluster = data_to_cluster, location_vars = location_vars, radius = radius)
  }
  
  # Convert to logDiameter
  data_to_cluster[, diameter_col] <- log(data_to_cluster[, diameter_col])
  
  
  ########## INITIALIZE
  # If no starting assignment is provided, initialize the sampler with each crater assigned to its own cluster
  if (is.null(starting_assignment) || all(starting_assignment == 1:ncustomers)) {
    current_table_assignments <- c(1:ncustomers)
    
    # Compute table specific quantities (means and covariances) and table sizes
    tablesize <- numeric(length = ncustomers)
    mus       <- list()
    sigmas    <- list()
    for(i in 1:max(current_table_assignments)){
      idx <- which(current_table_assignments == i)
      logD <- mean(data_to_cluster[idx, diameter_col])
      mns <- draw_table_means_covs(logD, mu0 = mu0, sigma0 = sigma0)
      mus[[i]] <- as.numeric(data_to_cluster[i, c(1:3)])
      sigmas[[i]] <- mns$sigma
      tablesize[i] <- length(idx)
    }
  } 
  # Otherwise, if the init_mus, init_sigmas, and starting assignment arguments are specified, start the sampler at that location
  else if (!is.null(init_mus) && !is.null(init_sigmas)){
    current_table_assignments <- starting_assignment
    
    # Compute table sizes and store the initial position of the table specific means and covariances
    tablesize <- numeric(length = ncustomers)
    mus       <- init_mus
    sigmas    <- init_sigmas
    for(i in 1:max(current_table_assignments)){
      idx <- which(current_table_assignments == i)
      tablesize[i] <- length(idx)
    }
  } 
  # Otherwise, begin the partition at the starting_assignment value and draw both the cluster-specific means and covariances from the prior
  else {
    current_table_assignments <- starting_assignment
    
    # Compute table specific quantities (means and covariances) and table sizes
    tablesize <- numeric(length = ncustomers)
    mus       <- list()
    sigmas    <- list()
    for(i in 1:max(current_table_assignments)){
      idx <- which(current_table_assignments == i)
      logD <- mean(data_to_cluster[idx, diameter_col])
      mns <- draw_table_means_covs(logD, mu0 = mu0, sigma0 = sigma0)
      mus[[i]] <- mns$mu
      sigmas[[i]] <- mns$sigma
      tablesize[i] <- length(idx)
    }
  }
  
  # Initialize the permutation of the data
  if (is.null(permutation)){
    permutation<-sample(1:ncustomers, ncustomers, replace = F)
  }
  else{
    permutation<-permutation
  }

  # Optional: Initialize the acceptance probabilities of the Metropolis-Hastings samplers
  perm_prob <- 0
  alpha_acc <- 0
  alpha_count <- 0
  sigma_acc <- 0
  sigma_count <- 0
  
  # Store the vector of family assignments for use later
  family_vec <- data_to_cluster[, family]
  
  # Initialize the CSV to store the results if specified via the print_vec argument
  if (!is.null(print_vec)){
    csv_file <- paste0(output_filename, ".csv")
    header <- c("iter", "type", paste0("v", seq_len(ncustomers)))
    write.table(
      t(header),
      file = csv_file,
      sep = ",",
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
    
    if (0 %in% print_vec){
    write_mu_draw(mus, current_table_assignments, paste0("mus ", csv_file))
    write_sigma_draw(sigmas, current_table_assignments, paste0("sigmas ", csv_file))
    
    write_row(0,
              "partition",
              match(current_table_assignments, unique(current_table_assignments)),
              csv_file)
    
    write_row(0,
              "permutation",
              permutation,
              csv_file)
    
    write_row(0,
              "alpha",
              rep(alpha, ncustomers),
              csv_file)
    }
  }
  
  ########## DO THE GIBBS
  # Iteration loop
  for(iter in 1:niter){
    # If the iteration is the last in a full scan, sample alpha
    if (iter %% ncustomers == 0){
      alpha_res<-alpha_dfcrp_sampler(family_vec, cluster_vec=current_table_assignments, alpha_int=alpha, alpha_prop_precision=alpha_prop_precision, alpha_prior_shape, alpha_prior_rate, permutation=permutation)
      alpha<-alpha_res[[1]]
      # If the iteration is outside a 20% burn-in window, record the result of the sampler for the acceptance ratio
      if(iter>=(0.2*niter)){
        alpha_acc<-alpha_res[[2]]+alpha_acc
        alpha_count <- alpha_count + 1
      }
    }
    
    # Sample the permutation
    permutation_res<-perm_sampler(family_vec, cluster_vec=current_table_assignments, alpha=alpha, permutation)
    permutation<-permutation_res[[1]]
    perm_prob<-permutation_res[[2]]+perm_prob
    
    # Customer "loop": because of our modified permutation proposal, we only update the last item in the partition
    for(i in ncustomers){ 
      
      # Identify which item is the last according to the permutation
      j=permutation[i]
      
      # Compute table quantities
      n_tables <- max(current_table_assignments)
      current_table <- current_table_assignments[j]
      
      # Compute family quantities
      same_family_tables <- NULL
      current_family <- data_to_cluster[j, family_col]
      same_family_tables <- unique(current_table_assignments[setdiff(which(data_to_cluster[,family_col]==current_family),j)])
      
      ########## REMOVE CUSTOMER J AND ADJUST 
      ## remove current customer from current table assignments
      current_table_assignments[j] <- NA
      
      ## decrease size of current table
      tablesize[current_table]    <- tablesize[current_table] - 1
      
      ## update current table size
      current_table_size <- tablesize[current_table]
      
      # Assign a table to be the first empty table
      if (tablesize[current_table] == 0) {
        first_empty <- current_table
      } else {
        empty_tables<-which(tablesize == 0)
        first_empty <- min(empty_tables)
      }
      
      ########## ELIMINATE POSSIBLE TABLES IF A NEIGHBORHOOD IS USED
      ## Get the neighbors of the current customer
      neighbor_tables <- which(tablesize!= 0)
      if(!is.null(radius)){
        current_neighbors <- which(neighbor_matrix[j,]==T)
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
      
      ########## COMPUTE LIKELIHOOD    
      ## temp_crp_prior is the current tablesize, with alpha put in at the first empty table
      temp_crp_prior<-rep(0, (length(tablesize)+1))
      for (k in temptables){
        temp_crp_prior[k] <- tablesize[k]
      }
      temp_crp_prior[first_empty] <- alpha
      
      ## compute the likelihood for all temptables
      templik <- rep(-Inf, (length(tablesize)+1))
      
      for(k in temptables){
        templik[k] <- dmvnorm(as.matrix(data_to_cluster[j, feature_cols]), mus[[k]], sigmas[[k]], log = T)
      }
      
      # If the item was the only item present in the cluster, the likelihood of that assignment comes from the current cluster-specific mean and covariance matrix
      if (tablesize[current_table] == 0) {
        
        # Assign the likelihood of a new cluster
        templik[first_empty] <- dmvnorm(data_to_cluster[j, feature_cols], mus[[current_table]], sigmas[[current_table]], log = T)
        
        # Normalize the log-probabilities
        m <- max(templik)
        likprobs <- exp(templik - m)
        likprobs <- likprobs / sum(likprobs)
        
        # Sample a new table
        newtable  <- sample(1:length(temp_crp_prior),1,prob=temp_crp_prior*likprobs)
        
        # Assign the current customer to the sampled table
        current_table_assignments[j] <- newtable
        
        # Increase the tablesize
        tablesize[newtable] <- tablesize[newtable] + 1
        
        # If the item doesn't return to its original table, remove that mean and covariance matrix
        if (newtable != current_table){
          sigmas[current_table] <- list(NULL)
          mus[current_table]    <- list(NULL)
        }
      } 
      # Otherwise, draw the cluster-specific means and covariance matrix from the prior
      else {
        # Draw the cluster-specific mean vector and covariance matrix from the prior
        mns <- draw_table_means_covs(data_to_cluster[j, diameter_col], mu0 = mu0, sigma0 = sigma0)
        new_mu <- mns$mu
        new_sigma <- mns$sigma
        
        # Assign the likelihood of a new cluster
        templik[first_empty] <- dmvnorm(data_to_cluster[j, feature_cols], new_mu, new_sigma, log = T)
        
        # Normalize the log-probabilities
        m <- max(templik)
        likprobs <- exp(templik - m)
        likprobs <- likprobs / sum(likprobs)
        
        # Draw a new table assignment
        newtable  <- sample(1:length(temp_crp_prior),1,prob=temp_crp_prior*likprobs)
        
        # Assign the current customer to the sampled table
        current_table_assignments[j] <- newtable
        
        # Increase the tablesize
        tablesize[newtable] <- tablesize[newtable] + 1
        
        # If the item is assigned to a new cluster, store the mean vector and covariance matrix
        if (newtable == first_empty){
          mus[[newtable]] <- new_mu
          sigmas[[newtable]] <- new_sigma
        }
      }
    } #--end customer loop 
    
    # If the iteration is the last in a full scan, update the cluster-specific parameters
    if (iter %% ncustomers == 0){
      # Loop over clusters
      for (i in which(tablesize != 0)){
        
        # Get the data of the items assigned to the current cluster
        idx <- which(current_table_assignments == i)
        data_list <- vector("list", length = length(idx))
        for (j in seq_along(idx)){
          data_list[[j]] <- as.matrix(data_to_cluster[idx[j], feature_cols])
        }      
        
        # Obtain values for the components of the proposed cluster-specific covariance matrix
        sigma_prop <- rmvnorm(1, c(sigmas[[i]][1,1], sigmas[[i]][3,3], sigmas[[i]][1,3]), sigma_prop_mat)
        Vxy_prop <- sigma_prop[1]
        Vd_prop  <- sigma_prop[2]
        Cv_prop  <- sigma_prop[3]
        
        # After ensuring the proposed matrix is PSD, compute the MH acceptance ratio
        if (Vxy_prop < 0 | Vd_prop < 0) {
          r <- -Inf
        } else if (abs(Cv_prop) > sqrt((Vxy_prop * Vd_prop)/2)) {
          r <- -Inf
        } else {
          sigma_prop <- matrix(c(Vxy_prop, 0, Cv_prop, 0, Vxy_prop, Cv_prop, Cv_prop, Cv_prop, Vd_prop), nrow = 3)
          r <- posterior_compute(data_list, mus[[i]], sigma_prop, mu0 = mu0, sigma0 = sigma0)-posterior_compute(data_list, mus[[i]], sigmas[[i]], mu0 = mu0, sigma0 = sigma0)
        }
        
        # Accept or reject the proposed value of the covariance matrix according to the acceptance ratio
        if (log(runif(1)) < r) {
          sigmas[[i]] <- sigma_prop
          sigma_acc = sigma_acc + 1
          sigma_count = sigma_count + 1
        } else {
          sigmas[[i]] <- sigmas[[i]]
          sigma_count = sigma_count + 1
        }
        
        # Compute inverses
        Sigma0_inv <- solve(sigma0)
        Sigma_inv  <- solve(sigmas[[i]])
        
        # Data summaries
        n    <- length(idx)
        xbar <- colMeans(data_to_cluster[idx, feature_cols])
        
        # Posterior covariance
        Sigma_star <- solve(Sigma0_inv + n * Sigma_inv)
        
        # Posterior mean
        m <- Sigma_star %*% (
          Sigma0_inv %*% mu0 +
            n * Sigma_inv %*% xbar
        )
        
        # Draw new cluster-specific mean vector
        mus[[i]] <- rmvnorm(1, m, Sigma_star)
      } #--- end cluster loop
    }
    
    # If the iteration is in the print_vec argument, write the results of the iteration to the .csv
    if (!is.null(print_vec)){
      if (iter %in% print_vec){
        write_mu_draw(mus, current_table_assignments, paste0("mus ", csv_file))
        write_sigma_draw(sigmas, current_table_assignments, paste0("sigmas ", csv_file))
        
        write_row(iter,
                  "partition",
                  match(current_table_assignments, unique(current_table_assignments)),
                  csv_file)
        
        write_row(iter,
                  "permutation",
                  permutation,
                  csv_file)
        
        write_row(iter,
                  "alpha",
                  rep(alpha, ncustomers),
                  csv_file)
      }
    }
    
    if(print_status == T){
      if(iter %in% print_vec){ print(paste0("Iteration ", iter, " completed at ", now())) }}
  } #--end iterations loop 
  
  # Optional: compute acceptance probabilities and write them to the .csv
  permutation_prop<-perm_prob/niter
  alpha_prop<-alpha_acc/alpha_count
  sigma_prop<-sigma_acc/sigma_count
  if (!is.null(print_vec)){
    write_row(iter,
              "acceptance",
              c(alpha_prop, permutation_prop, sigma_prop, rep(NA, ncustomers-3)),
              csv_file)
  }

  # Return the most recent vector of table assignments in canonical order
  return(match(current_table_assignments, unique(current_table_assignments)))
} #----dfcrp_sampler----#


################################
# HELPER FUNCTIONS:
################################
# Computes the posterior for use in the Metropolis ratio when assigning sigma matrices to a cluster
posterior_compute <- function(data, mu, sigma, Vxy_precision = 1, Vd_precision = 100, mu0, sigma0){
  
  # Extract variance values from covariance matrix
  Vxy <- sigma[1,1]
  Vd <- sigma[3,3]
  Cv <- sigma[1,3]
  logD <- mean(sapply(data, function(x) x[3]))
  
  
  # Incorporate the prior on these variances in the log posterior
  log_posterior <- 0
  log_posterior <- 2*dgamma(Vxy, shape = (logD^4.5 * 0.08)*Vxy_precision, scale = 1/Vxy_precision, log = T) + log_posterior
  log_posterior <- dgamma(Vd, shape = (logD^-0.8 * 0.124)*Vd_precision, scale = 1/Vd_precision, log = T) + log_posterior
  
  # Incorporate the covariance prior term contribution
  log_posterior <- dbeta(((Cv/(sqrt((Vxy * Vd)/2))) + 1)/2, 100, 100, log=T) - log(2) - log(sqrt((Vxy * Vd)/2)) + log_posterior
  
  # Incorporate the prior on mu in the log posterior
  log_posterior <- dmvnorm(mu, mu0, sigma0, log = T) + log_posterior
  
  # For each data point present in the cluster, evaluate it's fit based on the cluster mean vector and covariance matrix
  for (x in data){
    log_posterior <- dmvnorm(as.numeric(unname(x)), mu, sigma, log = T) + log_posterior
  }
  
  return(log_posterior)
}

################################
# Draws from the prior on cluster mus and sigmas
draw_table_means_covs <- function(logD, Vxy_precision = 1, Vd_precision = 100, mu0, sigma0){
  
  # Draw elements of the sigma prior draw
  Vxy <- rgamma(1, shape = (logD^4.5 * 0.08)*Vxy_precision, scale = 1/Vxy_precision)
  Vd <- rgamma(1, shape = (logD^-0.8 * 0.124)*Vd_precision, scale = 1/Vd_precision)
  rho <- (2 * rbeta(1, 100, 100)) - 1
  Cv  <- rho * sqrt((Vxy * Vd)/2)
  
  # Compile the elements into a single matrix
  sigma <- matrix(c(Vxy, 0, Cv, 0, Vxy, Cv, Cv, Cv, Vd), nrow = 3)
  
  mu <- rmvnorm(1, mu0, sigma0)
  
  return(list(mu = mu, sigma = sigma))
} #-----table_means_covs-------#

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
  distance_matrix_idx <- distance_matrix_idx == 1
  distance_matrix_idx <- Matrix(distance_matrix_idx, sparse = TRUE)
  
  return(distance_matrix_idx)
} #------distance_neighbor_maker--------#

#################################
# Writes a row to the csv file
write_row <- function(iter, type, values, csv_file) {
  df <- data.frame(
    iter = iter,
    type = type,
    t(values)
  )
  write.table(
    df,
    file = csv_file,
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE
  )
} #------------write_row--------------#

#################################
# Flattens a cluster covariance matrix and writes it to a .csv file
write_sigma_draw <- function(sigmas, current_table_assignments, file){
  
  # tables in order of appearance (same order used by match())
  tables <- unique(current_table_assignments)
  
  # grab sigmas in that same order
  sigmas_ordered <- lapply(tables, function(k) sigmas[[k]])
  
  # flatten each covariance matrix row-wise
  sigmas_matrix <- do.call(rbind, lapply(sigmas_ordered, function(s) as.vector(t(s))))
  
  # write cluster rows
  write.table(sigmas_matrix,
              file,
              sep = ",",
              row.names = FALSE,
              col.names = FALSE,
              append = TRUE)
  
  # blank separator row
  write.table(matrix(NA, nrow = 1, ncol = ncol(sigmas_matrix)),
              file,
              sep = ",",
              row.names = FALSE,
              col.names = FALSE,
              append = TRUE)
} #------write_sigma_draw--------#

#################################
# Writes a cluster mean vector to a .csv file
write_mu_draw <- function(mus, current_table_assignments, file){
  
  # tables in order of appearance (same as match())
  tables <- unique(current_table_assignments)
  
  # grab mus in that same order
  mus_ordered <- lapply(tables, function(k) mus[[k]])
  
  mus_matrix <- do.call(rbind, lapply(mus_ordered, as.vector))
  
  write.table(mus_matrix, file, sep=",",
              row.names=FALSE, col.names=FALSE, append=TRUE)
  
  # blank separator row
  write.table(matrix(NA, 1, ncol(mus_matrix)), file, sep=",",
              row.names=FALSE, col.names=FALSE, append=TRUE)
  
} #------write_mu_draw--------#

