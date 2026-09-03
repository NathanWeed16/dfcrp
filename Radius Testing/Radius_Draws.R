source("DFCRP_Gibbs.R")

#for i in {1..200}; do
#R CMD BATCH --no-save --no-restore "--args $i" Radius_Draws.R &
#done

# Read in the data and subset it for testing the radius
crater_data<-read.table("/Users/nathanweed/Research/DFCRP/Code/Crater_Meas_data.txt", header = T, sep = " ")
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

args=(commandArgs(TRUE))

args<-as.numeric(args)

set.seed(3102026 + args)

name <- paste0("Draws", args)

dfcrp_sampler(data_to_cluster, "Observer", c("X", "Y", "Diameter"), niter = 359*10000, mu0 = c(2050, -150, 3.2), sigma0 = matrix(c(200^2, 0, 0, 0, 110^2, 0, 0, 0, 0.5^2), nrow = 3), print_vec = 359*seq(2000, 10000, by = 10), output_filename = name, starting_assignment = c(1:359))
