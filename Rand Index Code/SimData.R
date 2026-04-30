library(mvtnorm)
library(MCMCpack)
################################
###---For Data---###
simulate_data <- function(K, J, X_limits = c(0, 1), Y_limits = c(0, 1), a_diameter = 5, b_diameter = 5, true_rate, error_rate, Sigma_meas = diag(c(5, 5, 0.01))) {
  
  ## ---- true crater parameters ----
  true_X <- runif(K, X_limits[1], X_limits[2])
  true_Y <- runif(K, Y_limits[1], Y_limits[2])
  true_log_D <- rgamma(K, a_diameter, b_diameter)
  
  mus <- cbind(true_X, true_Y, true_log_D)
  
  simulated_data <- list()
  false_cluster_idx <- K + 1
  
  for (j in seq_len(J)) {
    
    ## ---- real crater detections ----
    detected <- sample(
      seq_len(K),
      size = round(true_rate[j] * K),
      replace = FALSE
    )
    
    real_obs <- lapply(detected, function(i) {
      obs <- rmvnorm(1, mus[i, ], Sigma_meas)
      c(obs, Family = j, Cluster = i)
    })
    
    ## ---- false positives ----
    n_false <- round(error_rate[j] * K)
    
    false_obs <- lapply(seq_len(n_false), function(i) {
      fx <- runif(1, X_limits[1], X_limits[2])
      fy <- runif(1, Y_limits[1], Y_limits[2])
      fd <- rgamma(1, a_diameter, b_diameter)
      c(fx, fy, fd, Family = j, Cluster = false_cluster_idx + i - 1)
    })
    
    false_cluster_idx <- false_cluster_idx + n_false
    
    simulated_data[[j]] <- rbind(
      do.call(rbind, real_obs),
      do.call(rbind, false_obs)
    )
  }
  
  simulated_data <- do.call(rbind, simulated_data)
  simulated_data <- data.frame(simulated_data)
  
  names(simulated_data) <- c("X", "Y", "Diameter", "Family", "Cluster")
  
  simulated_data$Family  <- as.integer(simulated_data$Family)
  simulated_data$Cluster <- as.integer(simulated_data$Cluster)
  
  simulated_data <- simulated_data[order(simulated_data$Family), ]
  
  simulated_data[,3] <- exp(simulated_data[,3])
  
  return(simulated_data)
}
