min(crater_data[,1])
max(crater_data[,1])
hist(crater_data[,1])
(64.1624+4072.18)/2

sum(crater_data[,1] < 100)
14216/90
sum(rnorm(158, 2068.171, 920) < 100)

sum(crater_data[,1] > 4000)
14216/218
sum(rnorm(65, 2068.171, 920) > 4000)


min(crater_data[,2])
max(crater_data[,2])
hist(crater_data[,2])
(-2221.1+10.1056)/2

sum(crater_data[,2] < -2200)
14216/134
sum(rnorm(106, -1105.497, 600) < -2200)

sum(crater_data[,2] > -100)
14216/773
sum(rnorm(18, -1105.497, 600) > -100)


min(log(crater_data[,3]))
max(log(crater_data[,3]))
hist(log(crater_data[,3]))
(1.658459+6.227773)/2
mean(log(crater_data[,3]))
sd(log(crater_data[,3]))

sum(log(crater_data[,3]) < 2)
14216/14
sum(rnorm(1015, 3.8, 0.65) < 2)

sum(log(crater_data[,3]) > 5.5)
14216/61
sum(rnorm(233, 3.8, 0.65) > 5.5)


# So the analyses above would suggest a MVN like:
test_mat <- matrix(nrow = 200, ncol = 3)
for (i in 1:200){
test_mat[i, ] <- rmvnorm(1, c(2068.171, -1105.497, 3.8), matrix(c(920^2, 0, 0, 0, 600^2, 0, 0, 0, 0.65^2), nrow = 3))
}

p <- plot_craters(crater_data[,1], crater_data[,2], crater_data[,3])

p + add_craters(test_mat[,1], test_mat[,2], exp(test_mat[,3]))


# How likely are the big ones under the prior?
# Make a similar plot under the old prior?

## Ok now for the variance priors
draw <- rmvnorm(1, c(2068.171, -1105.497, 3.8), matrix(c(920^2, 0, 0, 0, 600^2, 0, 0, 0, 0.65^2), nrow = 3))

p + add_craters(c(2737.879, 2722), c(-1402.647, -1402.647), c(exp(3.575954), exp(3.575954)))

# The original is 2737.88. 2710 is too far. I think 2720 is too. I think 2722 is about as far apart as I'd allow

p + add_craters(c(2737.879, 2754), c(-1402.647, -1402.647), c(exp(3.575954), exp(3.575954)))

# That makes sense, because on the other hand, 2754 is about as far as I'd allow on the other side.

p + add_craters(c(2737.879, 2737.879), c(-1402.647, -1387), c(exp(3.575954), exp(3.575954)))
p + add_craters(c(2737.879, 2737.879), c(-1402.647, -1418.5), c(exp(3.575954), exp(3.575954)))

# The two plots above are just evidence that the same adjustment to the y axis is also acceptable

p + add_craters(c(2202.256, 2179), c(-1424.671, -1424.671), c(exp(4.6), exp(4.6)))

# I think 2179 is as far as I would allow. That's a max distance of 23.

# One equation that's pretty close is D^0.356 * 1.49. Let's test it

# exp(5)^0.356 * 1.49 = 8.8355. Max distance of 26.5.

p + add_craters(c(2202.256, 2165), c(-1424.671, -1424.671), c(exp(5), exp(5)))


# New rule: D^0.607 * 0.54 = SD. At D = exp(3), SD = 3.336.

p + add_craters(c(2202.256, 2193.256), c(-1424.671, -1424.671), c(exp(3), exp(3)))

# New rule: D^0.71 * 0.38 = SD. At D = exp(5.5), SD = 18.87.

p + add_craters(c(2202.256, 2145.65), c(-1424.671, -1424.671), c(exp(5.5), exp(5.5)))

# I think that's close enough
diameter <- exp(3)
c <- 1000

# This looks pretty good!
rgamma(1, shape = ((diameter^0.71 * 0.38)^2)*c, scale = 1/c)

p + add_craters(c(2202.256, 2202.256-38), c(-1424.671, -1424.671), c(exp(5), exp(5)))


# Ok, now for the diameter variance.

p + add_craters(c(2202.256, 2202.256-16), c(-1424.671+20, -1424.671), c(exp(5.2), exp(5)))

# Ok, 3 and 5 is too big a jump. Is 4 and 6 equally bad?
# Yeah. 
# Even 4 and 5 is a little much for me. I think 4.7 and 5 is as far out as I would go
# The problem is, I don't think 5.3 and 5 make sense. I think 0.25 units apart is a good balance.

p + add_craters(c(2202.256, 2202.256-6), c(-1424.671+10, -1424.671), c(exp(4.5), exp(4.1)))

# I think at a smaller diameter, I would allow a little more. Say 0.3 units different

# Ok the points I have (logD, SD): (5, 0.25/3); (4.7, 0.3/3); (4.1, 0.4/3)
# Rule: 0.3613 - 0.0556lD
# This would mean that at a logD of 6.2, SD = 0.014

p + add_craters2(c(2202.256, 2202.256-2), c(-1424.671+3, -1424.671), c(exp(5.2), exp(4.6)))


# This rule looks pretty good!
sqrt(rgamma(1, shape = ((0.352*(log(diameter))^-0.7)^2)*c, scale = 1/c))*3

# Now to combine them
lD <- 5
c1 <- 10
c2 <- 1000
dist1 <- rnorm(1)*sqrt(rgamma(1, shape = (lD^5.6 * 0.08)*c1, scale = 1/c1))
dist2 <- rnorm(1)*sqrt(rgamma(1, shape = (lD^5.6 * 0.08)*c1, scale = 1/c1))
dist3 <- rnorm(1)*sqrt(rgamma(1, shape = ((0.352*(lD)^-0.7)^2)*c2, scale = 1/c2))
p + add_craters2(c(2202.256-dist1, 2202.256), c(-1424.671+dist2, -1424.671), c(max(exp(lD), exp(lD-dist3)), min(exp(lD), exp(lD-dist3))))

# Done!
lD <- 4
Vd = rgamma(1, shape = (lD^-0.8 * 0.124)*c2, scale = 1/c2)
Vxy <- rgamma(1, shape = (lD^5.5 * 0.08)*c1, scale = 1/c1)
sig <- matrix(c(Vxy, 0, 0, 0, Vxy, 0, 0, 0, Vd), nrow = 3)

x <- runif(1, 0, 100)
y <- runif(1, 0, 100)
d <- runif(1)
test <- c(2000+x, -1000-y, lD+d)
iden <- dmvnorm(c(2000,-1000,lD), c(2000,-1000,lD), sig, log = T)
print(iden)
shift <- dmvnorm(test, c(2000,-1000,lD), sig, log = T)
print(shift)
exp(iden)/exp(shift)
p + add_craters2(c(2000, test[1]), c(-1000, test[2]), c(max(exp(lD), exp(test[3])), min(exp(lD), exp(test[3]))))

# 5.5, 1000, (400 times more likely)
# 5, 300, (400 times more likely)
# 4.5, 100, (457.1447 times more likely)
# 4, 52, (407.3244 times more likely)
# 3.5, 19, (372.7843 times more likely)
# 3, 8.5, (358.6332 times more likely)

# 5.5, 0.02, (404 times more likely)
# 5, 0.017, (400 times more likely)
# 4.5, 0.024, (300 times more likely)
# 4, 0.03, (407.3244 times more likely)
# 3.5, 19, (372.7843 times more likely)
# 3, 8.5, (358.6332 times more likely)


