library(ggplot2)
source("DFCRP_Gibbs.R")

set.seed(3162026)

#for i in {1..61}; do
#R CMD BATCH --no-save --no-restore "--args $i" Timing.R &
#done

# Read in the data and format the subset
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

# Read in the correct radius values
args=(commandArgs(TRUE))
radius_index<-as.numeric(args[1])
radii <- seq(0, 600, by = 10)
radius <- radii[radius_index]

# Test out the timing
val<-system.time(dfcrp_sampler(data_to_cluster, "Observer", c("X", "Y", "Diameter"), niter = 359*10000, mu0 = c(2050, -150, 3.2), sigma0 = matrix(c(200^2, 0, 0, 0, 110^2, 0, 0, 0, 0.5^2), nrow = 3), location_vars = c("X", "Y"), radius = radius))
res<-c("Radius:", radius, "Time:", val[[3]])

# Write the results to a .csv
write.table(res, file="Timing316.csv", sep = ",", append=T, row.names = F, col.names = F)

quit(save="no")

# Read in the data
timing_res<-read.table("Timing/Timing316.csv", header=F)

# Take out the text columns
text<-seq(1, 324, by=2)
timing_res<-timing_res[-text, ]

# Add another row for the time values
timing_res<-as.matrix(timing_res)
blank_row<-matrix(nrow = nrow(timing_res))
timing_res<-cbind(timing_res, blank_row)

# Move the times to the second column, to the row that corresponds to the appropriate radius
times<-seq(2, nrow(timing_res), by=2)
temp<-vector(length=35)
temp<-timing_res[times, 1]
timing_res<-timing_res[-times, ]
timing_res[(times/2), 2]<-temp

# Change the times to hours rather than seconds
timing_res[,2]<-(as.numeric(timing_res[,2])/3600)

# Truncate the radius values and order them
timing_res[,1]<-trunc(as.numeric(timing_res[,1]))
timing_res <- timing_res[order(as.numeric(timing_res[,1])), ]
timing_res<-as.data.frame(timing_res)

#Rename and reformat the matrix
names(timing_res) <- c("Radius", "Time")
timing_res$Radius <- as.numeric(timing_res$Radius)
timing_res$Time <- as.numeric(timing_res$Time)

radius_mat <- read.csv("Radius Testing/Radius_Results.csv", row.names = 1)
averages <- numeric(length = nrow(timing_res))
for (i in 1:ncol(radius_mat)){
  averages[i] <- mean(radius_mat[, i])
}

df <- data.frame(
  x = seq(10, 600, by = 10),
  averages = averages[2:61],
  timing = timing_res[2:61, 2]
)

# Scaling factor
scale_factor <- max(df$averages) / max(df$timing)

df$timing_scaled <- df$timing * scale_factor

# Plot the results
timing_plot <- ggplot(df, aes(x = x)) +
  geom_line(aes(y = averages), linetype = "dashed", linewidth = 1) +
  annotate("point",
           x = 75,
           y = averages[df$x == 75],
           shape = 18,
           size = 5,
           color = "grey") +
  
  #geom_point(aes(y = timing_scaled), size = 2, alpha = 0.6) +
  geom_smooth(
    aes(y = timing_scaled),
    method = "loess",
    se = FALSE,
    linewidth = 1, 
    span = 0.57, 
    color = "black"
  ) +
  geom_point(data = df, aes(y = timing_scaled), size = 2) +
  scale_y_continuous(
    name = "Proportion of Violating Clusters",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Time (hours)"
    )
  ) +
  labs(x = expression("Radius (" * rho * ")")) +
  theme_minimal() +
  geom_vline(aes(xintercept = 75), linetype = "dotted", linewidth = 0.9)

timing_plot

timing_plot <- timing_plot +
  theme(
    axis.text = element_text(size = 12),   # tick labels
    axis.title = element_text(size = 18)   # axis titles
  )

ggsave("Timing_Plot.pdf", plot = timing_plot, width = 9, height = 4.5)
