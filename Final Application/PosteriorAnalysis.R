library(flexclust)
library(reshape2)
library(ggplot2)
library(ggforce)

crater_data<-read.table("Crater_Meas_data.txt", header = T, sep = " ")
which_con <- which(crater_data$Observer == "Concensus")
crater_data <- crater_data[-which_con, ]
too_small <- which(crater_data$Diameter < 18)
crater_data <- crater_data[-too_small, ]
crater_data$Observer <- match(crater_data$Observer, c("Antonenko1", "Antonenko2", "Antonenko3", "Chapman", "Fassett", "Herrick", "Kirchoff", "Robbins1", "Robbins2", "Singer", "Zanetti"))

# List of experts sorted by number of craters identified
by_size <- c(1, 10, 3, 6, 5, 7, 11, 2, 9, 8, 4)

partitions <- fread("used_partitions.csv")

n <- nrow(partitions)

partitions_by_expert <- list()

for (i in 1:11){
  idx <- which(crater_data[,4] == i)
  partitions_by_expert[[i]] <- partitions[, ..idx]
}


jaccard_mat <- matrix(NA, 11, 11)

unique_partitions <- lapply(partitions_by_expert, function(mat) {
  apply(mat, 1, function(row) unique(row))
})

for(i in 1:10) {
  for(j in (i+1):11) {
    # Compute all rowwise Jaccards at once
    jacc_vec <- sapply(1:100000, function(n) {
      set_i <- unique_partitions[[i]][,n]
      set_j <- unique_partitions[[j]][,n]
      length(intersect(set_i, set_j)) / length(unique(c(set_i, set_j)))
    })
    
    jaccard_mat[i, j] <- mean(jacc_vec)
  }
}

jaccard_mat[lower.tri(jaccard_mat)] <- t(jaccard_mat)[lower.tri(jaccard_mat)]
diag(jaccard_mat) <- 1

avg_jacc <- sapply(1:11, function(i) mean(jaccard_mat[i, -i]))

# Average Jaccard similarity index sorted by size of expert
round(avg_jacc[by_size], 2)




# Jaccard similarity heatmap
pdf("Jaccard_Heatmap.pdf", width = 8, height = 6)

expert_letters <- LETTERS[by_size]

df <- melt(jaccard_mat)
colnames(df) <- c("Expert1", "Expert2", "Jaccard")

# Convert to factor with correct order
df$Expert1 <- factor(df$Expert1, levels = as.character(by_size), labels = expert_letters)
df$Expert2 <- factor(df$Expert2, levels = as.character(by_size), labels = expert_letters)

# Remove diagonal
df$Jaccard[df$Expert1 == df$Expert2] <- NA

ggplot(df, aes(Expert1, Expert2, fill = Jaccard)) +
  geom_tile(color = "gray80") +
  scale_fill_gradient(low = "white", high = "black", na.value = "white", name = "Jaccard\nSimilarity\nCoefficient") +
  theme_minimal() +
  labs(x = "Expert", y = "Expert")

dev.off()






# Preallocate
size1_mat <- matrix(0, nrow(partitions), ncol = 11)

# Expert indices list
expert_list <- lapply(1:11, function(i) which(crater_data[,4] == i))

# Function to count singleton intersections for one row
count_singletons <- function(row) {
  # Find positions that occur exactly once
  tbl <- tabulate(row)
  single_vals <- which(tbl == 1)   # values that appear once
  single_pos <- which(row %in% single_vals)
  
  # Count for each expert
  sapply(expert_list, function(expert) length(intersect(expert, single_pos)))
}

# Apply over all rows
size1_mat <- t(apply(partitions, 1, count_singletons))

size1_mat_prop <- size1_mat / rowSums(size1_mat)

round(colMeans(size1_mat_prop[,by_size])*100, 3)
round(apply(size1_mat_prop[,by_size], 2, function(x) quantile(x, c(0.025, 0.975)))*100, 3)




# Preallocate
size10_mat <- matrix(0, nrow(partitions), ncol = 11)

# Function to compute size-10 counts per row
count_size10_not_in_expert <- function(row) {
  # Count occurrences
  tbl <- table(row)
  vals10 <- as.numeric(names(tbl[tbl == 10]))  # actual values appearing exactly 10 times
  pos10 <- which(row %in% vals10)              # positions in row with these values

  # Count how many of these positions are NOT in each expert
  sapply(expert_list, function(expert) length(vals10) - length(intersect(pos10, expert)))
}

# Apply over all rows (fast enough for 25,000 rows)
size10_mat <- t(apply(partitions, 1, count_size10_not_in_expert))

size10_mat_prop <- size10_mat / rowSums(size10_mat)

# Column means
round(colMeans(size10_mat_prop[,by_size])*100, 3)

# 2.5% and 97.5% quantiles per expert
round(apply(size10_mat_prop[,by_size], 2, function(x) quantile(x, c(0.025, 0.975)))*100, 3)





# Experts B and D size 1 cluster plot
plot_craters2 <- function(x, y, diameter, highlight_idx = integer(0)) {
  
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
      alpha = 0.5
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

size1_list <- lapply(seq_len(nrow(partitions)), function(i) {
  vec <- as.numeric(partitions[i, ])
  
  counts <- tabulate(vec)
  singleton_clusters <- which(counts == 1)
  
  which(match(vec, singleton_clusters, nomatch = 0) > 0)
})

all_indices <- unlist(size1_list)

index_counts <- table(all_indices)

threshold <- 0.5 * length(size1_list)

size1 <- as.numeric(
  names(index_counts[index_counts >= threshold])
)

expert_d1 <- as.numeric(intersect(size1, expert_list[[4]]))
expert_b1 <- intersect(size1, expert_list[[2]])


pdf("SingletonsBD.pdf", width = 8, height = 6)

plot_craters2(crater_data[c(expert_b1, expert_d1),1], crater_data[c(expert_b1, expert_d1),2], crater_data[c(expert_b1, expert_d1),3], highlight_idx = 1:length(expert_b1))

dev.off()      






# Average permutation acceptance probability
perm_acc <- numeric(400)
for (i in seq(1, 400, by = 2)){
  
  file <- paste0("Final Application/Final Draws_cont/Final_Draws_2cont", (i+1)/2, ".csv")
  draws <- read.csv(file)
  perm_acc[i] <- draws[751, 4]
  
  file <- paste0("Final Application/Final Draws_cont/Final_Draws_4cont", (i+1)/2, ".csv")
  draws <- read.csv(file)
  perm_acc[i+1] <- draws[751, 4]
}

# 0.3354854
mean(perm_acc)





# Identified craters by size by expert (black histogram dots)
small_counts <- numeric()
for (i in 1:11){
  small_counts[i] <- length(which(crater_data[,3] >= 18 & crater_data[,3] < 50 & crater_data[,4] == i))
}

mid_counts <- numeric()
for (i in 1:11){
  mid_counts[i] <- length(which(crater_data[,3] >= 50 & crater_data[,3] < 100 & crater_data[,4] == i))  
}

large_counts <- numeric()
for (i in 1:11){
  large_counts[i] <- length(which(crater_data[,3] >= 100 & crater_data[,4] == i))  
}

# DBSCAN size 5 or greater values (diamonds on histogram)
dbscan <- c(754, 94, 41)




n4_small2 <- integer(200*500)
n4_medium2 <- integer(200*500)
n4_large2 <- integer(200*500)

n5_small2 <- integer(200*500)
n5_medium2 <- integer(200*500)
n5_large2 <- integer(200*500)

n6_small2 <- integer(200*500)
n6_medium2 <- integer(200*500)
n6_large2 <- integer(200*500)
idx <- 1
for (draw in 1:200){
  
  # Read in partitions and cluster means
  partitions1 <- read.csv(paste0("Final Application/Used Partitions/Final_Draws_2cont", draw, ".csv"))
  partitions2 <- read.csv(paste0("Final Application/Used Partitions/Final_Draws_4cont", draw, ".csv"))
  partitions1 <- partitions1[seq(1, 750, by = 3), c(3:ncol(partitions1))]
  partitions2 <- partitions2[seq(1, 750, by = 3), c(3:ncol(partitions2))]
  partitions <- rbind(partitions1, partitions2)
  partitions <- as.matrix(partitions)
  
  cluster_means1 <- read.csv(paste0("Final Application/Used Mus/mus Final_Draws_2cont", draw, ".csv"), header = F)
  cluster_means2 <- read.csv(paste0("Final Application/Used Mus/mus Final_Draws_4cont", draw, ".csv"), header = F)
  cluster_means <- rbind(cluster_means1, cluster_means2)
  
  is_blank <- apply(cluster_means, 1, function(x) all(is.na(x) | x == ""))
  
  # Indices of blank rows
  blank_rows <- which(is_blank)
  
  # Add start and end boundaries
  starts <- c(1, blank_rows + 1)
  ends <- c(blank_rows - 1, nrow(cluster_means))
  
  # Remove any empty ranges
  valid <- which(starts <= ends)
  starts <- starts[valid]
  ends <- ends[valid]
  
  # Now we can create a list of partitions
  n_parts <- length(starts)
  partitions_means <- vector("list", n_parts)
  for (j in seq_len(n_parts)) {
    partitions_means[[j]] <- cluster_means[starts[j]:ends[j], , drop = FALSE]
  }
  partitions_means <- tail(partitions_means, 500)
  
  # Define your thresholds for mean classification
  small_thresh <- 50
  medium_thresh <- 100
  
  n4_small <- integer(500)
  n4_medium <- integer(500)
  n4_large <- integer(500)
  
  n5_small <- integer(500)
  n5_medium <- integer(500)
  n5_large <- integer(500)
  
  n6_small <- integer(500)
  n6_medium <- integer(500)
  n6_large <- integer(500)
  
  for (i in 1:500) {
    # Get the current partition as a vector of cluster IDs
    part_i <- partitions[i, ]
    
    counts <- tabulate(part_i)
    clust_ids <- which(counts > 0)
    sizes <- counts[clust_ids]
    
    # Get mean values for the clusters in this partition
    means_i <- partitions_means[[i]]
    
    means_vec <- exp(means_i[clust_ids, 3])
    
    valid4 <- sizes >= 4
    valid5 <- sizes >= 5
    valid6 <- sizes >= 6
    
    n4_small[i]  <- sum(valid4 & means_vec < small_thresh)
    n4_medium[i] <- sum(valid4 & means_vec >= small_thresh & means_vec < medium_thresh)
    n4_large[i]  <- sum(valid4 & means_vec >= medium_thresh)
    
    n5_small[i]  <- sum(valid5 & means_vec < small_thresh)
    n5_medium[i] <- sum(valid5 & means_vec >= small_thresh & means_vec < medium_thresh)
    n5_large[i]  <- sum(valid5 & means_vec >= medium_thresh)
    
    n6_small[i]  <- sum(valid6 & means_vec < small_thresh)
    n6_medium[i] <- sum(valid6 & means_vec >= small_thresh & means_vec < medium_thresh)
    n6_large[i]  <- sum(valid6 & means_vec >= medium_thresh)
  }
  range <- idx:(idx + 499)
  n4_small2[range] <- n4_small
  n4_medium2[range] <- n4_medium
  n4_large2[range] <- n4_large
  
  n5_small2[range] <- n5_small
  n5_medium2[range] <- n5_medium
  n5_large2[range] <- n5_large
  
  n6_small2[range] <- n6_small
  n6_medium2[range] <- n6_medium
  n6_large2[range] <- n6_large
  idx <- idx + 500
  print(draw)
}

mean(n4_small2)
quantile(n4_small2, c(0.025, 0.975))
mean(n4_medium2)
quantile(n4_medium2, c(0.025, 0.975))
mean(n4_large2)
quantile(n4_large2, c(0.025, 0.975))

mean(n5_small2)
quantile(n5_small2, c(0.025, 0.975))
mean(n5_medium2)
quantile(n5_medium2, c(0.025, 0.975))
mean(n5_large2)
quantile(n5_large2, c(0.025, 0.975))

mean(n6_small2)
quantile(n6_small2, c(0.025, 0.975))
mean(n6_medium2)
quantile(n6_medium2, c(0.025, 0.975))
mean(n6_large2)
quantile(n6_large2, c(0.025, 0.975))



min(n5_small2)
max(n5_small2)



# Histogram data
hist_df <- data.frame(
  value = c(n5_small2, n5_medium2, n5_large2),
  group = factor(
    c(
      rep("18 px ≤ Diameter < 50 px", length(n5_small2)),
      rep("50 px ≤ Diameter < 100 px", length(n5_medium2)),
      rep("Diameter ≥ 100 px", length(n5_large2))
    ),
    levels = c(
      "18 px ≤ Diameter < 50 px",
      "50 px ≤ Diameter < 100 px",
      "Diameter ≥ 100 px"
    )
  )
)


dot_df <- data.frame(
  value = c(small_counts, mid_counts, large_counts),
  group = factor(
    c(
      rep("18 px ≤ Diameter < 50 px", length(small_counts)),
      rep("50 px ≤ Diameter < 100 px", length(mid_counts)),
      rep("Diameter ≥ 100 px", length(large_counts))
    ),
    levels = levels(hist_df$group)
  )
)

diamond_df <- data.frame(
  value = dbscan,
  group = factor(
    c(
      "18 px ≤ Diameter < 50 px",
      "50 px ≤ Diameter < 100 px",
      "Diameter ≥ 100 px"
    ),
    levels = levels(hist_df$group)
  )
)

cairo_pdf("ConsensusHistogram.pdf", width = 8, height = 3)

# Create a list of plots, one per group
plots <- lapply(unique(hist_df$group), function(g) {
  ggplot(subset(hist_df, group == g), aes(x = value)) +
    geom_histogram(binwidth = 1, fill = "gray39", color = NA) +
    geom_point(data = subset(dot_df, group == g), aes(x = value, y = 0),
               inherit.aes = FALSE, shape = 16) +
    geom_point(data = subset(diamond_df, group == g), aes(x = value, y = 0),
               inherit.aes = FALSE, shape = 18, size = 6, color = "gray71") +
    labs(x = "Number of Craters", y = NULL, title = g) +
    theme_bw() +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5)  # center the title
    )
})

# Combine horizontally
wrap_plots(plots, nrow = 1)

dev.off()


