library(ggforce)



##### Helper Functions #####
radius_violations <- function(cluster_assignments, neighbor_matrix) {
  
  n <- length(cluster_assignments)
  clusters <- split(seq_len(n), cluster_assignments)
  
  violating_clusters <- vector("list", 0)
  
  for (cl_label in names(clusters)) {
    
    cl <- clusters[[cl_label]]
    
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
      violating_clusters[[cl_label]] <- cl
    }
  }
  
  list(
    n_violations = length(violating_clusters),
    violating_clusters = violating_clusters
  )
}

plot_craters <- function(x, y, diameter, highlight_idx = integer(0)) {
  
  df <- data.frame(
    x = x,
    y = y,
    radius = diameter / 2,
    highlight = seq_along(x) %in% highlight_idx
  )
  
  ggplot(df) +
    
    # background craters
    geom_circle(
      data = subset(df, !highlight),
      aes(x0 = x, y0 = y, r = radius),
      fill = "grey70",
      color = NA,
      alpha = 0.1
    ) +
    
    # highlighted craters
    geom_circle(
      data = subset(df, highlight),
      aes(x0 = x, y0 = y, r = radius),
      fill = "black",
      color = NA,
      alpha = 1
    ) +
    
    coord_equal() +
    theme_void()
}

add_craters <- function(x, y, diameter,
                        alpha = 1) {
  
  df <- data.frame(
    x0 = x,
    y0 = y,
    r = diameter / 2,
    crater_id = factor(seq_along(x))
  )
  
  geom_circle(
    data = df,
    aes(x0 = x0, y0 = y0, r = r, fill = crater_id),
    fill = "black",
    alpha = alpha,
    inherit.aes = FALSE, 
    show.legend = FALSE
  )
}

add_craters2 <- function(x, y, diameter,
                        alpha = 1) {
  
  df <- data.frame(
    x0 = x,
    y0 = y,
    r = diameter / 2,
    crater_id = factor(seq_along(x))
  )
  
  geom_circle(
    data = df,
    aes(x0 = x0, y0 = y0, r = r, fill = crater_id),
    color = NA,
    alpha = alpha,
    inherit.aes = FALSE, 
    show.legend = FALSE
  )
}






rad_res <- read.csv("Radius_Results.csv")
xs <- data_to_cluster[,1]
ys <- data_to_cluster[,2]
diameters <- data_to_cluster[,3]

radius <- 20
burnin <- 0.198
to_keep <- floor(1000*(1-burnin))-1
index <- 1
vals <- list(2)

draws_per_chain <- list(200)
i = 0
j = 1
for (j in 1:200){
  draws_per_chain[[j]] <- rad_res[(i+1):(i+801), ]
  i = i + 801
}

for (j in 1:200){
  draws_per_chain[[j]] <- draws_per_chain[[j]][(nrow(draws_per_chain[[j]])-to_keep):nrow(draws_per_chain[[j]]), ]
}

rad_res <- do.call(rbind, draws_per_chain)
rad_res <- rad_res[,-1]
colnames(rad_res) <- seq(10, 80, by = 10)
# write.csv(rad_res, "Radius_Results309cleaned.csv")

violating_indices <- which(rad_res[,(radius/10)] != 0)
vals[[1]] <- ceiling(violating_indices[index] / (to_keep))
vals[[2]] <- 3*(violating_indices[index] %% (to_keep)) - 2
if (vals[[2]] == 0) vals[[2]] <- 3*(to_keep + 1) - 2

name <- paste0("Draws", vals[[1]], ".csv")
draw <- read.csv(name)
partition <- as.numeric(draw[vals[[2]], 3:ncol(draw)])


nonos <- radius_violations(partition, distance_neighbor_maker(data_to_cluster, c("X", "Y"), radius = radius))


violation_plot <- plot_craters(xs, ys, diameters, highlight_idx = as.numeric(unlist(nonos[[2]])))

ggsave(filename = "Violation.pdf", plot = violation_plot, width = 10, height = 8, units = "in")

# # Helpful for finding indices on the big dataset from those of the reduced dataset
# idx <- unname(apply(data_to_cluster[as.numeric(unlist(nonos[[2]])),1:3], 1, function(x) {
#   which(
#     crater_data[,1] == x[1] &
#       crater_data[,2] == x[2] &
#       crater_data[,3] == x[3]
#   )[1]
# }))

# Additionally, here is the code for the craters which actually dictate a radius of 75:
plot_craters(crater_data[,1], crater_data[,2], crater_data[,3], highlight_idx = which(crater_data[,3] < 400 & crater_data[,3] > 280 & crater_data[,2] < -2000))

max(dist(crater_data[which(crater_data[,3] < 400 & crater_data[,3] > 280 & crater_data[,2] < -2000), c(1:2)]))
