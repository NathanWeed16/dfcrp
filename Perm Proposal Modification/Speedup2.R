source("Old_DFCRP_Gibbs.R")
library(coda)

#for i in {1..30}; do
#R CMD BATCH --no-save --no-restore "--args $i" Speedup2.R &
#done

args=(commandArgs(TRUE))
args<-as.numeric(args)
set.seed(3142026 + args)

# Read in the data
crater_data<-read.table("Crater_Meas_data.txt", header = T, sep = " ")
which_con <- which(crater_data$Observer == "Concensus")
crater_data <- crater_data[-which_con, ]
too_small <- which(crater_data$Diameter < 18)
crater_data <- crater_data[-too_small, ]
crater_data$Observer <- match(crater_data$Observer, c("Antonenko1", "Antonenko2", "Antonenko3", "Chapman", "Fassett", "Herrick", "Kirchoff", "Robbins1", "Robbins2", "Singer", "Zanetti"))

# Define a function that will run the sampler, record the runtime, compute the effective sample size, and return the timing, ESS,
# and states to initialize the next run if necessary.
compute_ess <- function(niter, job_id, run_id, partition, permutation, alpha, init_mus, init_sigmas) {
  out_file <- paste0("ESS2_job_", job_id, "_run_", run_id, ".csv")
  
  t0 <- Sys.time()
  dfcrp_sampler(
    crater_data,
    family = "Observer",
    features = c("X","Y","Diameter"),
    location_vars = c("X","Y"),
    radius = 75,
    niter = niter,
    output_filename = paste0("ESS2_job_", job_id, "_run_", run_id),
    print_vec = c(1:niter),
    starting_assignment = partition,
    permutation = permutation,
    alpha = alpha, 
    init_mus = init_mus,
    init_sigmas = init_sigmas
  )
  t1 <- Sys.time()
  
  runtime <- as.numeric(t1 - t0, units="secs")
  
  ESS <- read.csv(out_file)
  
  partition <- as.numeric(unlist(ESS[nrow(ESS)-3, 3:ncol(ESS)]))
  permutation_mat <- as.matrix(ESS[seq(2, nrow(ESS), by = 3), c(3:ncol(ESS))])
  alpha <- as.numeric(ESS[nrow(ESS)-1, 3])
  
  mu_file <- paste0("mus ", out_file)
  mus <- read.csv(mu_file, header = F)
  na_rows <- which(apply(mus, 1, function(x) all(is.na(x))))
  
  if (length(na_rows) == 1){
    latest_mu_draws <- mus[1:(na_rows - 1), ]
  } else {
    last_na <- tail(na_rows, 2)[1]
    latest_mu_draws <- mus[(last_na + 1):(nrow(mus) - 1), ]
  }
  
  init_mus <- lapply(1:nrow(latest_mu_draws), function(i) as.numeric(latest_mu_draws[i, ]))
  
  sigma_file <- paste0("sigmas ", out_file)
  sigmas <- read.csv(sigma_file)
  na_rows <- which(apply(sigmas, 1, function(x) all(is.na(x))))
  
  if (length(na_rows) == 0){
    latest_sigma_draws <- sigmas[1:(na_rows - 1), ]
  } else {
    last_na <- tail(na_rows, 2)[1]
    latest_sigma_draws <- sigmas[(last_na + 1):(nrow(sigmas) - 1), ]
  }
  
  init_sigmas <- lapply(1:nrow(latest_sigma_draws), function(i) {
    matrix(as.numeric(latest_sigma_draws[i, ]), nrow = 3, byrow = TRUE)
  })
  
  return(list(
    time = runtime,
    partition = partition,
    permutation_mat = permutation_mat,
    init_mus = init_mus,
    init_sigmas = init_sigmas,
    alpha = alpha
  ))
}

# Initialize values
target_ess <- 5
ess <- 0
niter <- 5
total_time <- 0
run_counter <- 1
partition <- NULL
permutation <- NULL
init_mus <- NULL
init_sigmas <- NULL
alpha <- 40
total_vals <- c()

while (ess < target_ess) {
  
  # Run the sampler, then compute ESS and runtime
  res <- compute_ess(niter, args, run_id = run_counter, 
                     partition, permutation, alpha, init_mus, init_sigmas)
  for (i in 1:nrow(res$permutation_mat)){
    total_vals <- c(total_vals, mean(res$permutation_mat[i, 1:(floor(ncol(res$permutation_mat)/2))]))
  }

  ess <- effectiveSize(total_vals)
  
  total_time <- total_time + res$time
  
  # Store initialization values
  partition <- as.numeric(res$partition)
  permutation <- as.numeric(res$permutation_mat[nrow(res$permutation_mat), ])
  init_mus <- res$init_mus
  init_sigmas <- res$init_sigmas
  alpha <- as.numeric(res$alpha)
  
  run_counter <- run_counter + 1
}

# Write both the ESS and the time to a CSV
write.csv(c(total_time, ess), paste0("Speedup2Res_", args, ".csv"), append = T)

quit(save = "no")
