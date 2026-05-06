library(microbenchmark)
Rcpp::sourceCpp("DFCRP_Prior_Optimized.cpp")

#for i in {1..200}; do
#R CMD BATCH --no-save --no-restore "--args $i" PriorTiming.R &
#done

args=(commandArgs(TRUE))
args<-as.numeric(args)
set.seed(3132026 + args)

crater_data<-read.table("Crater_Meas_data.txt", header = T, sep = " ")
which_con <- which(crater_data$Observer == "Concensus")
crater_data <- crater_data[-which_con, ]
too_small <- which(crater_data$Diameter < 18)
crater_data <- crater_data[-too_small, ]
crater_data$Observer <- match(crater_data$Observer, c("Antonenko1", "Antonenko2", "Antonenko3", "Chapman", "Fassett", "Herrick", "Kirchoff", "Robbins1", "Robbins2", "Singer", "Zanetti"))

append_mean_time <- function(mean_time,
                             file = "PriorTiming.csv",
                             job_id = NA,
                             extra = NULL) {
  
  # Build row
  row <- data.frame(
    job_id = job_id,
    mean_time_sec = mean_time,
    extra = if (is.null(extra)) NA else extra
  )
  
  # If file doesn't exist, write with header
  if (!file.exists(file)) {
    write.table(
      row,
      file = file,
      sep = ",",
      row.names = FALSE,
      col.names = TRUE,
      append = FALSE
    )
  } else {
    # Append without header
    write.table(
      row,
      file = file,
      sep = ",",
      row.names = FALSE,
      col.names = FALSE,
      append = TRUE
    )
  }
}

draws <- read.csv(paste0("Final_Draws_cont", args, ".csv"))
idx <- seq(1, 750, by = 3)

for (i in seq_along(idx)){
  
  part <- as.numeric(draws[idx[i], 3:ncol(draws)])
  perm <- as.numeric(draws[idx[i]+1, 3:ncol(draws)])
  alpha <- draws[idx[i]+2, 3]
  
  mb <- microbenchmark(
    log_dfcrp_pmf_cond_cpp(crater_data[,4], cluster_vec = part, alpha = alpha, permutation = perm),
    times = 100
  )
  
  append_mean_time(mean_time = mean(mb$time)/1e9, job_id = args)
}


quit(save = "no")

times <- read.csv("PriorTiming.csv")

# 0.005880282
mean(times[,2])





