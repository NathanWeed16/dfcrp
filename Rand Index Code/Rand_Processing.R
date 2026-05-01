library(data.table)

# ---------- Helper function ----------
process_partition_file <- function(name, row_seq){
  
  draws <- fread(name)
  
  rand_partitions <- draws[row_seq(nrow(draws)), 3:ncol(draws)]
  
  n <- nrow(rand_partitions)
  
  out <- matrix(0, n, 6)
  
  for (k in 1:n){
    
    tab <- table(as.numeric(rand_partitions[k,]))
    
    # Faster than looping l = 1:6
    out[k, ] <- colSums(outer(tab, 1:6, ">="))
    
  }
  
  return(out)
}

# ---------- DFCRP ----------
counts_list <- vector("list", 125 * 4)
idx <- 1

for (i in 1:125){
  for (j in 1:4){
    
    name <- paste0(
      "DFCRP Partitions/DFCRP_Rand_job",
      i,
      "_run",
      j,
      ".csv"
    )
    
    counts_list[[idx]] <- process_partition_file(
      name,
      function(n) seq(4, n-1, by = 3)
    )
    
    idx <- idx + 1
  }
}

counts_mat <- do.call(rbind, counts_list)



# ---------- DFCRP wRad ----------
rad_counts_list <- vector("list", 125 * 4)
idx <- 1

for (i in 1:125){
  for (j in 1:4){
    
    name <- paste0(
      "DFCRP wRad Partitions/DFCRP_Rad_Rand_job",
      i,
      "_run_",
      j,
      ".csv"
    )
    
    rad_counts_list[[idx]] <- process_partition_file(
      name,
      function(n) seq(4, n-1, by = 3)
    )
    
    idx <- idx + 1
  }
}

rad_counts_mat <- do.call(rbind, rad_counts_list)



# ---------- CRP ----------
crp_counts_list <- vector("list", 125 * 4)
idx <- 1

for (i in 1:125){
  for (j in 1:4){
    
    name <- paste0(
      "CRP Partitions/CRP_Rand_job",
      i,
      "_run_",
      j,
      ".csv"
    )
    
    crp_counts_list[[idx]] <- process_partition_file(
      name,
      function(n) seq(4001, n-1, by = 2)
    )
    
    idx <- idx + 1
  }
}

crp_counts_mat <- do.call(rbind, crp_counts_list)



# ---------- Final summary ----------
full_mat <- rbind(
  colMeans(counts_mat),
  colMeans(rad_counts_mat),
  colMeans(crp_counts_mat)
)

write.csv(full_mat, "Counting_Summary.csv", row.names = FALSE)