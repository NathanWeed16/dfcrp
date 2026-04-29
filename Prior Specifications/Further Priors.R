min(data_to_cluster[,1])
max(data_to_cluster[,1])
hist(data_to_cluster[,1])
(1700.16+2400)/2

sum(data_to_cluster[,1] < 1800)
475/47
sum(rnorm(10, 2050, 200) < 1800)

sum(data_to_cluster[,1] > 2350)
475/36
sum(rnorm(13, 2050, 200) > 2350)

min(data_to_cluster[,2])
max(data_to_cluster[,2])
hist(data_to_cluster[,2])

sum(data_to_cluster[,2] < -280)
475/44
sum(rnorm(11, -150, 110) < -280)

sum(data_to_cluster[,2] > -20)
475/33
sum(rnorm(14, -150, 110) > -20)


min(log(data_to_cluster[,3]))
max(log(data_to_cluster[,3]))
hist(log(data_to_cluster[,3]))
(2+5.269)/2
mean(log(data_to_cluster[,3]))

sum(log(data_to_cluster[,3]) < 2.5)
475/31
sum(rnorm(15, 3.15, 0.5) < 2.5)
sum(log(data_to_cluster[,3]) > 3.6)
475/60
sum(rnorm(8, 3.15, 0.5) > 3.6)

# So the analyses above would suggest a MVN like:
test_mat <- matrix(nrow = 50, ncol = 3)
for (i in 1:50){
  test_mat[i, ] <- rmvnorm(1, c(2050, -150, 3.15), matrix(c(200^2, 0, 0, 0, 110^2, 0, 0, 0, 0.5^2), nrow = 3))
}

p <- plot_craters(data_to_cluster[,1], data_to_cluster[,2], data_to_cluster[,3])

p + add_craters(test_mat[,1], test_mat[,2], exp(test_mat[,3]))





# For Rand Index calculations
source("Rand Index Code/SimData.R")
simdat <- simulate_data(K=30, J=6, X_limits = c(0,700), Y_limits = c(0,500), a_diameter = 64,b_diameter = 16,
                        true_rate=c(0.98, 0.96, 0.94, 0.92, 0.9, 0.88), error_rate=c(0.12, 0.1, 0.08, 0.06, 0.04, 0.02))



(min(simdat[,1])+max(simdat[,1]))/2
(min(simdat[,2])+max(simdat[,2]))/2
(min(simdat[,3])+max(simdat[,3]))/2

sum(simdat[,1] > 600)
180/29
sum(rnorm(6, 350, 300) > 600)

sum(simdat[,1] < 100)
180/18
sum(rnorm(10, 350, 300) < 100)


sum(simdat[,2] > 450)
180/22
sum(rnorm(8, 250, 225) > 450)

sum(simdat[,2] < 50)
180/35
sum(rnorm(5, 250, 225) < 50)


sum(simdat[,3] < 3.5)
180/26
sum(rnorm(7, 3.9, 0.45) < 3.5)

sum(simdat[,3] > 4.5)
180/13
sum(rnorm(14, 3.9, 0.45) > 4.5)


# So the analyses above would suggest a MVN like:
test_mat <- matrix(nrow = 100, ncol = 3)
for (i in 1:100){
  test_mat[i, ] <- rmvnorm(1, c(350, 250, 3.9), matrix(c(300^2, 0, 0, 0, 225^2, 0, 0, 0, 0.45^2), nrow = 3))
}

p <- plot_craters(simdat[,1], simdat[,2], exp(simdat[,3]))

p + add_craters(test_mat[,1], test_mat[,2], exp(test_mat[,3]))
