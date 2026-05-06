library(coda)
library(ggplot2)
library(dplyr)
library(patchwork)
#for i in {1..200}; do
#R CMD BATCH --no-save --no-restore "--args $i" TheBigOne.R &
#done

Rprof("TheFinal")
args=(commandArgs(TRUE))
args<-as.numeric(args)
set.seed(3132026 + args)

crater_data<-read.table("Crater_Meas_data.txt", header = T, sep = " ")
which_con <- which(crater_data$Observer == "Concensus")
crater_data <- crater_data[-which_con, ]
too_small <- which(crater_data$Diameter < 18)
crater_data <- crater_data[-too_small, ]
crater_data$Observer <- match(crater_data$Observer, c("Antonenko1", "Antonenko2", "Antonenko3", "Chapman", "Fassett", "Herrick", "Kirchoff", "Robbins1", "Robbins2", "Singer", "Zanetti"))

source("DFCRP_Gibbs.R")

filename <- paste0("Final_Draws", args)

results <- dfcrp_sampler(crater_data, family = "Observer", features = c("X", "Y", "Diameter"), location_vars = c("X", "Y"), radius = 75, niter=9517*2500, print_vec = 9517*seq(10, 2500, by = 10), output_filename = filename, alpha_prior_shape = 3, alpha_prior_rate = 0.04)

Rprof(NULL)
out <- summaryRprof("TheFinal")

capture.output(out, file = paste0("FinalTime_", args, ".txt"))

quit(save = "no")







