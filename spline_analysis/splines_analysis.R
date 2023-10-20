library(splines2)

ages <- 0:100
ages_no0 <- ages[-1]

# Construction of the spline basis functions
S.temp <- bSpline(ages_no0, 
                  knots = c(seq(5, 40, by = 5), 50, 60, seq(70, 95, by = 5)), 
                  degree = 2, 
                  intercept = TRUE)
# The number of B-spline basis is equal to the number of internal nodes plus
#   the order, that is equal to degree plus one. Hence, in this case
#   the number is 16+3=19.
K <- ncol(S.temp)

# Normalization
for (k in 1:K){
  S.temp[,k] <- S.temp[,k]/max(S.temp[,k])
}

S <- S.temp

# Plot of B-splines except for the Dirac delta
matplot(x = 1:100, y = S, type = "l")
knots <- c(seq(5, 40, by = 5), 50, 60, seq(70, 95, by = 5))
abline(v = knots)


b1 <- rnorm(19)
b2 <- rnorm(19)
b2[10:12] <- b1[10:12]

f1 <- S %*% b1
f2 <- S %*% b2

par(mfrow = c(1, 2))
plot(f1, type = "l", col = 1)
lines(f2, col = 2)
abline(v = knots, 
       col = 3, lty = 2)
middle_points <- (c(0, knots[1:16]) + c(knots[1:16], 100))/2
# points(middle_points, rep(0, length(knots) + 1), pch = 20, 
       # col = ifelse(b1 == b2, "blue", "red"))
ages_max <- apply(S, 2, which.max)
text(x = ages_max, y = rep(0, length(knots) + 1), 
     round(b1 - b2, 2),
     cex = 0.5)

matplot(x = 1:100, y = S, type = "l")
abline(v = knots, 
       col = 3, lty = 2)
points(middle_points, rep(0, length(knots) + 1), pch = 20, 
       col = ifelse(b1 == b2, "blue", "red"))
