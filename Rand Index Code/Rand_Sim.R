library(flexclust)
source("DFCRP_Gibbs.R") # Contains the dfcrp_sampler function
source("SimData.R") # Includes the simulate_family_data function
source("CRP_Gibbs.R") # Contains the crp_sampler function with the alpha prior included

#for i in {1..125}; do
#R CMD BATCH --no-save --no-restore "--args $i" Rand_Sim.R &
#done

args=(commandArgs(TRUE))

args<-as.numeric(args)

set.seed(3102026 + args)

################################
#' Computes the rand index between simulated data and posterior samples from the crp or dfcrp.
#' 
#' @param table_assign Posterior draws from the function dfcrp_sampler(). 
#' @param simdat The output of the function simulate_family_data().
#' 
#' @return Vector of the rand index for each sample from the posterior.
#' 
#' 
rand_index_sampler <- function(table_assign,simdata){
  randIndex.df <- NULL
  for(i in 1:nrow(table_assign)){
    randIndex.df <- c(randIndex.df,as.numeric(randIndex(table_assign[i,],simdata$Cluster,correct=T)))
  }
  return(randIndex.df)
} #---rand_index_sampler----#

################################
#' Computes the average percentage of clusters over multiple draws from the posterior that have more than one observational crater from a single expert assigned to it-thus forming an illegal cluster.
#' @param crp_df The output of the function dfcrp_sampler with DFCRP=TRUE
#' @param simdat  The output of the function simulate_family_data(). 
#' 
#' @return The average number of percentages of illegal assignments from a sample of cluster assignments.
#' 
#' 
illegal_assignments <- function(crp_df,simdat){
  
  # Create a vector to store the number of illegal assignments
  illegal <- NULL
  
  # Loop over each sample from the posterior
  for(i in 1:nrow(crp_df)){
    tempassign <- crp_df[i,]
    tempsim <- data.frame(cluster=tempassign,family=simdat$Family)
    unqcluster <- unique(tempassign)
    tempncluster <- length(unqcluster)
    tempillegal <- 0
    for(j in 1:tempncluster){
      tempcluster <- unqcluster[j]
      tempdat     <- subset(tempsim,cluster==tempcluster)
      
      if(nrow(tempdat) != length(unique(tempdat$family))){
        tempillegal <- tempillegal + 1
      }
    }
    illegal <- c(illegal, tempillegal/tempncluster)
  }
  return(mean(illegal))
}


################################
#' Simulates data and computes the rand index.
#'
#' @param  simnumber An identifier for the simulated data set. 
#' @param K  The number of clusters.
#' @param  X_limits Vector of the max and min X values.
#' @param Y_limitsVector of the max and min Y values.
#' @param  a_d  Shape parameter for the gamma distribution that samples diameter.
#' @param  b_d Rate parameter for the gamma distribution that samples diameter.
#' @param J  (integer) The number of families.
#' @param p (vector) Vector of probabilities that each observation will belong to each family. Must be of length J. If entries do not sum to 1, values will be adjusted so that they do.
#' @param niter The number of Gibbs iterations
#' 
#' @return List of [[1]] randoutput, where randoutput includes (in this order): simulation number, the mean and standard deviation of the distribution rand indexes between the DFCRP samples and the truth, the mean, median, and standard deviation of the distribution rand indexes between the CRP samples and the truth, the mean, median, and standard deviation of the distribution of differences between DFCRP and CRP between each methods samples and the truth, the average percentage of "illegal" craters from each sample from the CRP, where illegal is defined to have more than one observational crater from a given expert within the same cluster, and [[2]] simdata, which is a data.frame with columns of the three simulated features, the family and cluster assignment, and the simulation number (names = c(X, Y, Diameter, Family, Cluster, simdat)).
#' 
#' 
simrand <- function(arg, simnumber, K, J, X_limits, Y_limits, a_diameter, b_diameter, true_rate, error_rate, niter, mu0, sigma0){
  
  ## simulate the data
  simdat <- simulate_data(K, J, X_limits, Y_limits, a_diameter, b_diameter, true_rate, error_rate)
  # Determine the size of the largest table
  km     <- max(table(simdat$Family))
  
  ## define burnin
  burnin <- floor(niter*0.2)
  
  filename <- paste0("DFCRP_Rand_job", arg, "_run", simnumber)
  
  ## Fit DFCRP
  res<-dfcrp_sampler(subset(simdat,select=c("X","Y","Diameter","Family")),
                          family="Family",
                          features=c("X","Y","Diameter"),
                          niter=niter,
                          mu0 = mu0,
                          sigma0 = sigma0,
                          print_vec = seq(burnin, niter, by = 30),
                          output_filename = filename,
                          starting_assignment = c(1:nrow(simdat)))
  
  dfcrp_df <- read.csv(paste0(filename, ".csv"), header = T)
  dfcrp_df <- dfcrp_df[seq(1, nrow(dfcrp_df)-1, by = 3), 3:ncol(dfcrp_df)]
  
  filename <- paste0("DFCRP_Rad_Rand_job", arg, "_run_", simnumber)
  
  ## Fit DFCRP w radius
  dfcrp_rad_df<-dfcrp_sampler(subset(simdat,select=c("X","Y","Diameter","Family")),
                          family="Family",
                          features=c("X","Y","Diameter"),
                          niter=niter,
                          mu0 = mu0,
                          sigma0 = sigma0,
                          location_vars = c("X", "Y"),
                          radius = 75,
                          print_vec = seq(burnin, niter, by = 30),
                          output_filename = filename,
                          starting_assignment = c(1:nrow(simdat)))
  
  dfcrp_rad_df <- read.csv(paste0(filename, ".csv"), header = T)
  dfcrp_rad_df <- dfcrp_rad_df[seq(1, nrow(dfcrp_rad_df)-1, by = 3), 3:ncol(dfcrp_rad_df)]
  
  filename <- paste0("CRP_Rand_job", arg, "_run_", simnumber)
  
  # ## Fit CRP
  crp_df<-crp_sampler(subset(simdat,select=c("X","Y","Diameter","Family")),
                      family="Family",
                      features=c("X","Y","Diameter"),
                      niter=(niter/30), 
                      print_status=F, 
                      mu0 = mu0, 
                      sigma0 = sigma0,
                      starting_assignment = c(1:nrow(simdat)),
                      output_filename = filename,
                      print_vec = c(1:(niter/30)))
  
  # Number of iterations that ran
  niter2 <- nrow(crp_df)
  
  ## Illegal cluster assignments for CRP
  illassign <- illegal_assignments(crp_df[(0.2*niter2):niter2,], simdat)
  
  ## Get Rand Index DFCRP
  dfcrp_rand <- rand_index_sampler(as.matrix(dfcrp_df[1:nrow(dfcrp_df), ]), simdat)
  
  ## Get Rand Index DFCRP (with rad)
  dfcrp_rad_rand <- rand_index_sampler(as.matrix(dfcrp_rad_df[1:nrow(dfcrp_rad_df), ]), simdat)

  ## Get Rand Index CRP
  crp_rand <- rand_index_sampler(crp_df[(0.2*niter2):niter2,], simdat)
  
  ## mean difference
  output <- data.frame(post_mean_dfcrp = mean(dfcrp_rand),
                       post_mean_rad_dfcrp = mean(dfcrp_rad_rand),
                       post_mean_crp = mean(crp_rand),
                       post_mean_diff = mean(dfcrp_rand-crp_rand),
                       perc_crp_illegal = illassign)
  
  return(list(randoutput = output,
              simdata    = data.frame(simdat,simdat = simnumber), 
              count_vec  = avg_count_vec,
              count_rad_vec = avg_rad_count_vec,
              count_crp_vec = avg_crp_count_vec))
  
}

output_mat<-matrix(nrow = 4, ncol = 5)
for (i in 1:4){
   output<-simrand(args, i, K=30, J=6, X_limits = c(0,700), Y_limits = c(0,500), a_diameter = 64, b_diameter = 16, true_rate=c(0.98, 0.96, 0.94, 0.92, 0.9, 0.88), error_rate=c(0.12, 0.1, 0.08, 0.06, 0.04, 0.02), niter=10000*30, mu0 = c(350, 250, 3.9), sigma0 = matrix(c(300^2, 0, 0, 0, 225^2, 0, 0, 0, 0.45^2), nrow = 3))
   output_mat[i, 1:5]<-as.numeric(output[[1]][1, 1:5])
}

write.table(output_mat, file='My_Rand.csv', sep = ",", append=T, row.names = F, quote = F)

quit(save="no")

res<-read.csv('Rand Index Code/My_Rand.csv')
# Keep only rows where all entries are numeric (you can ignore the warnings)
numeric_rows <- apply(res, 1, function(row) all(!is.na(as.numeric(row))))

# Subset the matrix to those rows
res <- res[numeric_rows, , drop = FALSE]

res <- apply(res, 2, as.numeric)

# Initialize a matrix for our results to be stored
sum_mat<-matrix(nrow=4, ncol=5)
rownames(sum_mat)<-c("DFCRP", "DFCRP w/Rad", "CRP", "Diff")
colnames(sum_mat)<-c("Min", "25%", "Mean", "75%", "Max")

# Fill in the matrix with the summary results
sum_mat[1,1:5]<-c(min(res[, 1]), quantile(res[, 1], 0.25), mean(res[, 1]), quantile(res[, 1], 0.75), max(res[, 1]))
sum_mat[2,1:5]<-c(min(res[, 2]), quantile(res[, 2], 0.25), mean(res[, 2]), quantile(res[, 2], 0.75), max(res[, 2]))
sum_mat[3,1:5]<-c(min(res[, 3]), quantile(res[, 3], 0.25), mean(res[, 3]), quantile(res[, 3], 0.75), max(res[, 3]))
sum_mat[4,1:5]<-c(min(res[, 4]), quantile(res[, 4], 0.25), mean(res[, 4]), quantile(res[, 4], 0.75), max(res[, 4]))

# In the paper we also reference this figure
mean(res[, 5])

# Write the summary to a .csv
write.table(sum_mat, file="My_Rand_Summary.csv", sep = ",", quote = F)

# See Rand_Processing.R for processing of partitions to produce the cluster counts table