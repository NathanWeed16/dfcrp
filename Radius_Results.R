library(Matrix)

# R CMD BATCH Radius_Results.R &

set.seed(3102026)
#### Helper Functions

count_radius_violations <- function(data_to_cluster, cluster_assignments, neighbor_matrix) {
  
  n <- length(cluster_assignments)
  clusters <- split(seq_len(n), cluster_assignments)
  
  n_violations <- 0L
  
  for (cl in clusters) {
    
    # Singletons are always valid
    if (length(cl) <= 1) next
    
    subgraph <- neighbor_matrix[cl, cl, drop = FALSE]
    
    visited <- rep(FALSE, length(cl))
    stack <- 1L
    visited[1] <- TRUE
    
    while (length(stack) > 0) {
      v <- stack[[1]]
      stack <- stack[-1]
      
      neighbors <- which(subgraph[v, ] & !visited)
      if (length(neighbors)) {
        visited[neighbors] <- TRUE
        stack <- c(stack, neighbors)
      }
    }
    
    if (!all(visited)) {
      n_violations <- n_violations + 1L
    }
  }
  
  return(n_violations)
}

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
}

partitions <- vector("list", length = 160200)

for (i in 1:200){
  name <- paste0("Draws", i, ".csv")
  draws <- read.csv(name)
  for (j in seq(1, 2403, by = 3)){
    partitions[[((i-1)*801)+((j+2)/3)]] <- as.numeric(draws[j, 3:ncol(draws)])
  }
}

crater_data<-read.table("Crater_Meas_data.txt", header = T, sep = " ")
which_con <- which(crater_data$Observer == "Concensus")
crater_data <- crater_data[-which_con, ]
too_small <- which(crater_data$Diameter < 18)
crater_data <- crater_data[-too_small, ]
crater_data$Observer <- match(crater_data$Observer, c("Antonenko1", "Antonenko2", "Antonenko3", "Chapman", "Fassett", "Herrick", "Kirchoff", "Robbins1", "Robbins2", "Singer", "Zanetti"))

# Subset (X,Y) for testing the effect of the radius
Xmin      <- 1700
Xmax      <- 2400
Ymin      <- -300
Ymax      <- 0
data_to_cluster <- subset(crater_data, Image=="NAC" & X>=Xmin & X<=Xmax & Y>= Ymin & Y <= Ymax)

vals <- seq(10, 80, by = 10)
n_parts <- length(partitions)
n_vals  <- length(vals)
radius_mat <- matrix(NA, nrow = n_parts, ncol = n_vals)

denoms <- numeric(n_parts)
for (i in 1:n_parts){
  denoms[i] <- count_radius_violations(data_to_cluster, partitions[[i]], distance_neighbor_maker(data_to_cluster, c("X", "Y"), 0))/max(partitions[[i]])
}

for (i in seq_along(vals)){
  neighbor_matrix <- distance_neighbor_maker(data_to_cluster, c("X", "Y"), vals[i])
  for (j in 1:n_parts){
    radius_mat[j, i] <- (count_radius_violations(data_to_cluster, partitions[[j]], neighbor_matrix)/max(partitions[[j]])/denoms[j])
  }
}

write.csv(radius_mat, "Radius_Results.csv")

# See Violating_Clusters.R for cleanup and analysis