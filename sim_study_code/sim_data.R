rm(list = ls())
source("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/utils/funzioni.R")





### ### ### ### ## #
#### Parameters ####
### ### ### ### ## #
n <- 3
T_final <- 10
ages <- 0:100
knots <- c(20, 40)
randomseed <- 17
sigma_vec <- c(0.05, 0.1, 0.25); names(sigma_vec) <- paste0("unit", 1:n)




### ### ### ### #
#### Splines ####
### ### ### ### #
library(splines2)
set.seed(randomseed)
ages_no0 <- ages[-1]

# Construction of the spline basis functions
S.temp <- bSpline(ages_no0, 
                  knots = knots, 
                  degree = 2, 
                  intercept = TRUE)
# The number of B-spline basis is equal to the number of internal nodes plus
#   the order, that is equal to degree plus one.
K <- ncol(S.temp)

# Normalization
for (k in 1:K){
  S.temp[,k] <- S.temp[,k]/max(S.temp[,k])
}

# Add Dirac delta in 0 as first spline
S <- rbind(c(1, rep(0, k)),
           cbind(rep(0, length(ages_no0)), S.temp)
)





### ### ### ### ### ### ### ###
#### Var-cov matrix of GP ####
### ### ### ### ### ### ### ###
grid_yr <- expand.grid(1:T_final, 1:T_final)
val_ker <- apply(grid_yr, 1, function(x) sq_exp_ker(x[1] - x[2], 3/2))
Sigma <- matrix(val_ker, T_final, T_final, byrow = F)
rownames(Sigma) <- colnames(Sigma) <- paste0("t", 1:T_final)





### ### ### ### ### ### ### ### ###
#### Sample betas from the GPs ####
### ### ### ### ### ### ### ### ###
beta <- list(t(cbind(t(rmvnorm(1, -2 + (1:T_final) * (-0.02), 0.05^2*Sigma)),
                     t(rmvnorm(1, -3.2 + (1:T_final) * (-0.02), 0.05^2*Sigma)),
                     t(rmvnorm(1, -4.4 + (1:T_final) * (-0.02), 0.05^2*Sigma)))),
             
             t(cbind(t(rmvnorm(1, -4.5 + (1:T_final) * (-0.02), 0.1^2*Sigma)),
                     t(rmvnorm(1, -5.8 + (1:T_final) * (-0.02), 0.1^2*Sigma)),
                     t(rmvnorm(1, -7.1 + (1:T_final) * (-0.02), 0.1^2*Sigma)))),
             
             t(cbind(t(rmvnorm(1, -3 + (1:T_final) * (-0.02), 0.1^2*Sigma)),
                     t(rmvnorm(1, -4 + (1:T_final) * (-0.02), 0.1^2*Sigma)),
                     t(rmvnorm(1, -5 + (1:T_final) * (-0.02), 0.1^2*Sigma)))),
             
             t(cbind(t(rmvnorm(1, -4 + (1:T_final) * (-0.02), 0.1^2*Sigma)),
                     t(rmvnorm(1, -4.6 + (1:T_final) * (-0.02), 0.1^2*Sigma)),
                     t(rmvnorm(1, -5.2 + (1:T_final) * (-0.02), 0.1^2*Sigma)))),
             
             t(cbind(t(rmvnorm(1, -1 + (1:T_final) * (-0.02), 0.01^2*Sigma)),
                     t(rmvnorm(1, -1.7 + (1:T_final) * (-0.02), 0.01^2*Sigma)),
                     t(rmvnorm(1, -2.5 + (1:T_final) * (-0.02), 0.01^2*Sigma)))),
             
             t(cbind(t(rmvnorm(1, 0 + (1:T_final) * (-0.02), 0.02^2*Sigma)),
                     t(rmvnorm(1, -0.5 + (1:T_final) * (-0.02), 0.02^2*Sigma)),
                     t(rmvnorm(1, -1.5 + (1:T_final) * (-0.02), 0.02^2*Sigma)))))
names(beta) <- paste0("spline", 1:length(beta))
beta <- lapply(beta, 
               function(x) {dimnames(x) <- list(paste("cluster", 1:n, sep = ""),
                                                paste("t", 1:T_final, sep = "")); x})




### ### ### ### ### ###
#### Sample labels ####
### ### ### ### ### ###
# Fix the labels to pre-determined patterns
labs <- list(# 1 and 3 together in the firs half
             # All in different groups in the second half
             matrix(c(rep(1, T_final),
                      rep(3, T_final), 
                      rep(1, T_final/2), rep(2, T_final/2)),
                    nrow = n,
                    byrow = T),
             # All in different in groups in the first half
             # 2 and 3 together in the second half
             matrix(c(rep(1, T_final),
                      rep(3, T_final), 
                      rep(2, T_final/2), rep(3, T_final/2)),
                    nrow = n,
                    byrow = T),
             # 1 and 3 together in the firs half
             # 2 and 3 together in the second half
             # It's a mix of the first two settings
             matrix(c(rep(1, T_final),
                      rep(3, T_final/2), rep(2, T_final/2),
                      rep(1, T_final/2), rep(2, T_final/2)),
                    nrow = n,
                    byrow = T),
             # 2 and 3 together in the first half
             # 1 and 2 together in the second half
             matrix(c(rep(3, T_final/2), rep(2, T_final/2),
                      rep(2, T_final),
                      rep(2, T_final/2), rep(1, T_final/2)),
                    nrow = n,
                    byrow = T),
             # Random
             t(replicate(n, 
                         rep(sample(n, T_final/2, replace = T), each = 2))),
             # All in different groups all the time, but in a different order
             #  in the first and second half
             matrix(c(rep(1, T_final/2), rep(2, T_final/2),
                      rep(2, T_final/2), rep(3, T_final/2),
                      rep(3, T_final/2), rep(1, T_final/2)),
                    nrow = n,
                    byrow = T))
names(labs) <- paste0("spline", 1:length(labs))
labs <- lapply(labs,
               function(x) {dimnames(x) <- list(paste("unit", 1:n, sep = ""),
                                                paste("t", 1:T_final, sep = "")); x})




### ### ### ### ### ### ##
#### Compute the data ####
### ### ### ### ### ### ##
# Beta divided for countries
beta_c1 <- sapply(1:ncol(S), function(j) beta[[j]][cbind(labs[[j]][1, ], 1:T_final)])
beta_c2 <- sapply(1:ncol(S), function(j) beta[[j]][cbind(labs[[j]][2, ], 1:T_final)])
beta_c3 <- sapply(1:ncol(S), function(j) beta[[j]][cbind(labs[[j]][3, ], 1:T_final)])

# Compute simulated data
data1 <- beta_c1 %*% t(S) + rnorm(length(ages) * T_final, 
                                  mean = 0, 
                                  sd = sigma_vec[1])
data2 <- beta_c2 %*% t(S) + rnorm(length(ages) * T_final, 
                                  mean = 0, 
                                  sd = sigma_vec[2])
data3 <- beta_c3 %*% t(S) + rnorm(length(ages) * T_final,
                                   mean = 0, 
                                   sd = sigma_vec[3])
data <- list(t(data1), t(data2), t(data3))
names(data) <- paste("unit", 1:n, sep = "")
data <- lapply(data, 
               function(x) {dimnames(x) <- list(paste("age", ages, sep = ""),
                                                paste("t", 1:T_final, sep = "")); x})





### ### ### ### ### ### ### ## # 
#### Plots to check results ####
### ### ### ### ### ### ### ## #
# Plot of labels assignment
windows();par(mfrow = c(2, 3))
for (j in 1:6){matplot(t(labs[[j]]), type = "p", 
                       col = c(1, 2, 4), pch = c("-", "|", "/"),
                       cex = 2, ylim = c(3, 1))}

# Plot of betas for each country
windows(); par(mfrow = c(2, 3)); for (j in 1:6){
  tmp <- cbind(beta[[j]][cbind(labs[[j]][1, ], 1:T_final)],
               beta[[j]][cbind(labs[[j]][2, ], 1:T_final)],
               beta[[j]][cbind(labs[[j]][3, ], 1:T_final)])
  matplot(tmp, lty = 1, type = "b", col = c(1, 2, 4), pch = c("-", "|", "/"), cex = 2)
}


# Plot of simulated data (as function of the age)
windows(); par(mfrow = c(1, 3)); for (i in 1:n){
  matplot(data[[i]], type = "l",
          xlab = "age", ylab = "log-mort. rate",
          ylim = c(-8, 1),
          main = paste("Unit", i))
}


# Plot of simulated data (as function of the time)
par(mfrow = c(1, 3))
for (i in 1:3){
  matplot(t(data[[i]]), type = "l")
}





### ### ### ### ### ### ### ### ### ### ##
#### Simulation of the two covariates ####
### ### ### ### ### ### ### ### ### ### ##


# Var-cov matrix of covariates
val_ker_x <- apply(grid_yr, 1, function(x) sq_exp_ker(x[1] - x[2], 1))
Sigma_x <- matrix(val_ker, T_final, T_final, byrow = F)
rownames(Sigma_x) <- colnames(Sigma_x) <- paste0("t", 1:T_final)

# I don't know if it's better to simulate covariates from GP or from 
#   independent Gaussians. I leave GP for now, w/ other options commented.

# x1 <- matrix(c(rnorm(T_final, mean = -5, sd = 0.5),
#                rnorm(T_final, mean = 3, sd = 0.5),
#                rnorm(T_final/2, mean = -5, sd = 0.5), rnorm(T_final/2, mean = 0, sd = 0.5)),
#              nrow = T_final, ncol = n,
#              byrow = F)
par(mfrow = c(1, 1))
x1 <- t(cbind(t(rmvnorm(1, rep(-5, T_final), 0.5^2*Sigma_x)),
              t(rmvnorm(1, rep(3, T_final), 0.5^2*Sigma_x)),
              t(rmvnorm(1, rep(c(-5, -1), each = T_final/2), 0.5^2*Sigma_x))
              )
        )
dimnames(x1) <- list(paste("unit", 1:n, sep = ""),
                     paste("t", 1:T_final, sep = ""))
matplot(t(x1), type = "b", lty = 1, col = c(1, 2, 4))


# What if one covariate is NOT Gaussian? Does the similarity function
#   still work?
# Function to simulate from a gamma distribution from mean and sd
myrgamma <- function(n, mean, sd){
  rgamma(n, shape = mean^2/sd^2, scale = sd^2/mean)
}

x2 <- t(matrix(c(myrgamma(T_final, mean = 1, sd = 0.1),
               myrgamma(T_final, mean = 7, sd = 0.1),
               myrgamma(T_final/2, mean = 4, sd = 0.1), myrgamma(T_final/2, mean = 7, sd = 0.1)),
             nrow = T_final, ncol = n,
             byrow = F)
        )
dimnames(x2) <- list(paste("unit", 1:n, sep = ""),
                     paste("t", 1:T_final, sep = ""))

matplot(t(x2), type = "b", lty = 1, col = c(1, 2, 4))





### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### Sim. of 2 cov. according to G-G sampling model ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### 
v2 <- 1
s2 <- 1

val_ker_x <- apply(grid_yr, 1, function(x) sq_exp_ker(x[1] - x[2], 3))
Sigma_x <- matrix(val_ker, T_final, T_final, byrow = F)
rownames(Sigma_x) <- colnames(Sigma_x) <- paste0("t", 1:T_final)

mu <- rbind(rmvnorm(1, rep(-5, T_final), v2*Sigma_x),
            rmvnorm(1, rep(0, T_final), v2*Sigma_x),
            rmvnorm(1, rep(5, T_final), v2*Sigma_x))
matplot(t(mu), type = "l", lty = 1, col = 1:3, xlab = "Time", ylab = "mu")

x1 <- sapply(1:T_final, function(t) mu[labs[[1]][, t], t]) + rnorm(n*T_final, 0, s2)
x2 <- sapply(1:T_final, function(t) mu[labs[[2]][, t], t]) + rnorm(n*T_final, 0, s2)

matplot(t(x1), type = "l", lty = 1, col = 1:3, xlab = "Time", ylab = "mu")
matplot(t(x2), type = "l", lty = 1, col = 1:3, xlab = "Time", ylab = "mu")





### ### ### ### ### ### ##
#### Object to return ####
### ### ### ### ### ### ##
out <- list(data = data,
            beta = beta,
            labs = labs,
            x1 = x1,
            x2 = x2,
            Sigma_x = Sigma_x,
            Sigma_beta = Sigma,
            T_final = T_final,
            n = n,
            sigma_vec = sigma_vec,
            knots = knots,
            randomseed = randomseed)

save(out, file = "C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/data/simdata_PPMx/simulated_data.RData")
