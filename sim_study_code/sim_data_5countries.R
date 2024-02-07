rm(list = ls())
source("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/utils/funzioni.R")




### ### ### ### ## #
#### Parameters ####
### ### ### ### ## #
n <- 5
T_final <- 10
ages <- 0:100
knots <- c(20, 40)
randomseed <- 290497
sigma_vec <- c(0.05, 0.05, 0.05, 0.05, 0.05); names(sigma_vec) <- paste0("unit", 1:n)




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
p <- ncol(S)




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
min_betas <- c(-4.4, -7.1, -5, -5.4, -2.5, -1.5)
max_betas <- c(-2, -4.5, -3, -3.8, -1, 0)
delta <- c(0.05, 0.1, 0.05, 0.05, 0.05, 0.05)
beta <- replicate(p,
                  matrix(NA, n, T_final),
                  simplify = F)
for (j in 1:p){
  for (i in 1:n){
    beta[[j]][i,] <- seq(min_betas[j], max_betas[j], length = n)[i] +
      drop(rmvnorm(1, 0 + (1:T_final) * (-0.02), delta[j]^2*Sigma))
  }
  dimnames(beta[[j]]) <- list(paste("cluster", 1:n, sep = ""),
                              paste("t", 1:T_final, sep = ""))
}

# for (j in 1:p){matplot(t(beta[[j]]), type = "l", lty = 1)}
names(beta) <- paste0("spline", 1:length(beta))
beta <- lapply(beta, 
               function(x) {dimnames(x) <- list(paste("cluster", 1:n, sep = ""),
                                                paste("t", 1:T_final, sep = "")); x})



### ### ### ### ### ###
#### Sample labels ####
### ### ### ### ### ###
# Fix the labels to pre-determined patterns
labs <- list(
  #1
  matrix(c(rep(1, T_final/2), rep(1, T_final/2),
           rep(1, T_final/2), rep(1, T_final/2), 
           rep(1, T_final/2), rep(3, T_final/2),
           rep(2, T_final/2), rep(4, T_final/2),
           rep(2, T_final/2), rep(5, T_final/2)),
         nrow = n,
         ncol = T_final,
         byrow = T),
  #2
  matrix(c(rep(1, T_final/2), rep(2, T_final/2),
           rep(2, T_final/2), rep(5, T_final/2),
           rep(4, T_final/2), rep(2, T_final/2),
           rep(3, T_final/2), rep(5, T_final/2),
           rep(4, T_final/2), rep(5, T_final/2)),
         nrow = n,
         ncol = T_final,
         byrow = T),
  #3
  matrix(c(rep(1, T_final/2), rep(1, T_final/2),
           rep(2, T_final/2), rep(2, T_final/2),
           rep(3, T_final/2), rep(3, T_final/2),
           rep(4, T_final/2), rep(4, T_final/2),
           rep(5, T_final/2), rep(5, T_final/2)),
         nrow = n,
         ncol = T_final,
         byrow = T),
  #4
  matrix(c(rep(3, T_final/2), rep(3, T_final/2),
           rep(3, T_final/2), rep(3, T_final/2),
           rep(3, T_final/2), rep(3, T_final/2),
           rep(3, T_final/2), rep(3, T_final/2),
           rep(3, T_final/2), rep(3, T_final/2)),
         nrow = n,
         ncol = T_final,
         byrow = T),
  #5
  matrix(c(rep(1, T_final/2), rep(3, T_final/2),
           rep(1, T_final/2), rep(4, T_final/2),
           rep(2, T_final/2), rep(5, T_final/2),
           rep(1, T_final/2), rep(5, T_final/2),
           rep(3, T_final/2), rep(5, T_final/2)),
         nrow = n,
         ncol = T_final,
         byrow = T),
  #6
  matrix(c(rep(sample(1:n, n, TRUE), each = 2),
           rep(sample(1:n, n, TRUE), each = 2),
           rep(sample(1:n, n, TRUE), each = 2),
           rep(sample(1:n, n, TRUE), each = 2),
           rep(sample(1:n, n, TRUE), each = 2)),
         nrow = n,
         ncol = T_final,
         byrow = T)
  )
names(labs) <- paste0("spline", 1:length(labs))
labs <- lapply(labs,
               function(x) {dimnames(x) <- list(paste("unit", 1:n, sep = ""),
                                                paste("t", 1:T_final, sep = "")); x})
# par(mfrow = c(2, 3))
# for (j in 1:p){
#   matplot(t(labs[[j]]), type = "p", ylim = c(0, n+1), cex = 1.2)
# }



### ### ### ### ### ### ##
#### Compute the data ####
### ### ### ### ### ### ##
# Beta divided for countries
beta_countries <- lapply(1:n, 
                         function(i) sapply(1:p, 
                                            function(j) beta[[j]][cbind(labs[[j]][i, ], 1:T_final)]))

# Compute simulated data
data <- lapply(1:n, 
               function(i) t(beta_countries[[i]] %*% t(S) + 
                               rnorm(length(ages) * T_final, 
                                     mean = 0, 
                                     sd = sigma_vec[i])))
names(data) <- paste("unit", 1:n, sep = "")
data <- lapply(data, 
               function(x) {dimnames(x) <- list(paste("age", ages, sep = ""),
                                                paste("t", 1:T_final, sep = "")); x})





### ### ### ### ### ### ### ## # 
#### Plots to check results ####
### ### ### ### ### ### ### ## #
# Plot of labels assignment
# windows()
pdf("simulated_labels.pdf", width = 16, height = 10)
par(mar = c(4, 4.5, 2, 1),
    oma = rep(0, 4))
layout(matrix(c(1, 1, 1,
                2:4,
                5:7,
                8, 8, 8), ncol = 3, byrow = T),
       heights=c(0.5, 1, 1, 0.6))
plot.new()
text(0.5, 0.3, 
     "Simulated cluster assignment.pdf", cex=1.7, font=2)
for (j in 1:p){matplot(t(labs[[j]]), type = "p", 
                       col = c(1:4, 7),
                       ylim = c(0.5, 5.5),
                       pch = c("_", "|", "/", "\\", "*"),
                       cex = 2,
                       xlab = "time",
                       ylab = "cluster assignment",
                       main = bquote(beta[.(j)]))}
plot(1, type = "n", axes = FALSE, xlab="", ylab="")
legend("center", 
       legend = paste("Unit", 1:n, sep = ""), 
       col = c(1:4, 7),
       cex = 1.2,
       pt.cex = 2,
       pch = c("_", "|", "/", "\\", "*"),
       inset = 0,
       horiz = T,
       title = "Legend")
dev.off()

# Plot of betas for each country
windows(); par(mfrow = c(2, 3)); for (j in 1:6){
  tmp <- sapply(beta_countries, 
                function(x) x[, j])
  matplot(tmp, lty = 1, type = "b", 
          col = c(1:4, 7), 
          pch = c("-", "|", "/", "\\", "*"),
          cex = 2)
}


# Plot of simulated data (as function of the age)
windows(); par(mfrow = c(2, 3)); for (i in 1:n){
  matplot(data[[i]], type = "l",
          xlab = "age", ylab = "log-mort. rate",
          ylim = c(-8, 1),
          main = paste("Unit", i))
}


# Plot of simulated data (as function of the time)
windows(); par(mfrow = c(2, 3))
for (i in 1:n){
  matplot(t(data[[i]]), type = "l")
}





### ### ### ### ### ### ### ### ### ### ##
#### Simulation of the two covariates ####
### ### ### ### ### ### ### ### ### ### ##


# Var-cov matrix of covariates
val_ker_x <- apply(grid_yr, 1, function(x) sq_exp_ker(x[1] - x[2], 1))
Sigma_x <- matrix(val_ker, T_final, T_final, byrow = F)
rownames(Sigma_x) <- colnames(Sigma_x) <- paste0("t", 1:T_final)


x1_mean <- seq(-2, 2, length = n)
x1 <- matrix(x1_mean[labs[[1]]],
             nrow = n,
             ncol = T_final,
             byrow = F) + 
  rmvnorm(n, rep(0, T_final), 0.1^2*Sigma_x)
dimnames(x1) <- list(paste("unit", 1:n, sep = ""),
                     paste("t", 1:T_final, sep = ""))
par(mfrow = c(1, 2))
matplot(t(x1), type = "b", lty = 1, col = c(1:4, 7), pch = c("-", "|", "/", "\\", "*"))
matplot(t(labs[[1]]), type = "p", 
        col = c(1:4, 7), 
        pch = c("-", "|", "/", "\\", "*"),
        cex = 2, ylim = c(1, 5))

x2_mean <- seq(-0.2, 0.2, length = n)
x2 <- matrix(x1_mean[labs[[2]]],
             nrow = n,
             ncol = T_final,
             byrow = F) + 
  rmvnorm(n, rep(0, T_final), 0.05^2*Sigma_x)
dimnames(x1) <- list(paste("unit", 1:n, sep = ""),
                     paste("t", 1:T_final, sep = ""))
par(mfrow = c(1, 2))
matplot(t(x2), type = "b", lty = 1, col = c(1:4, 7), pch = c("-", "|", "/", "\\", "*"))
matplot(t(labs[[2]]), type = "p", 
        col = c(1:4, 7), 
        pch = c("-", "|", "/", "\\", "*"),
        cex = 2, ylim = c(1, 5))
dimnames(x2) <- list(paste("unit", 1:n, sep = ""),
                     paste("t", 1:T_final, sep = ""))





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
            delta_vec = delta,
            knots = knots,
            randomseed = randomseed)

save(out, file = "C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/data/simdata_no_hierarchy/simulated_data_5countries.RData")
