rm(list = ls())

mysample <- function(m){
  min <- 1
  res <- c()
  for (j in 1:m){
    if (length(min:m) == 1){
      s <- m
    } else {
      s <- sample(min:m, 1, prob = (min:m)[order(min:m, decreasing = T)])
      min <- max(1, s)
    }
    
    res <- c(res, s)
  }
  return(res)
}


simulateData <- function(n, T_final, 
                         ages, knots,
                         sigma_vec,
                         randomseed # seed for set.seed()
                         ){
  library(splines2)
  set.seed(randomseed)
  ages_no0 <- ages[-1]
  
  # Construction of the spline basis functions
  S.temp <- bSpline(ages_no0, 
                    knots = knots, 
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
  
  # Add Dirac delta in 0 as first spline
  S <- rbind(c(1, rep(0, k)),
             cbind(rep(0, length(ages_no0)), S.temp)
  )
  
  
  # beta <- list(cbind(rnorm(T_final, -1 + (1:T_final) * (-0.05), 0.1),
  #                    rnorm(T_final, -1 + (1:T_final) * (-0.1), 0.1),
  #                    rnorm(T_final, -2 + (1:T_final) * 0, 0.1),
  #                    rnorm(T_final, -3 + (1:T_final) * (-0.05), 0.1),
  #                    rnorm(T_final, -4 + (1:T_final) * 0, 0.1)),
  #              
  #              cbind(rnorm(T_final, -2 + (1:T_final) * (-0.05), 0.1),
  #                    rnorm(T_final, -2 + (1:T_final) * (-0.1), 0.1),
  #                    rnorm(T_final, -3 + (1:T_final) * 0, 0.1),
  #                    rnorm(T_final, -4 + (1:T_final) * (-0.05), 0.1),
  #                    rnorm(T_final, -5 + (1:T_final) * 0, 0.1)),
  #              
  #              cbind(rnorm(T_final, -4 + (1:T_final) * 0, 0.1),
  #                    rnorm(T_final, -4.8 + (1:T_final) * 0, 0.1),
  #                    rnorm(T_final, -5.5 + (1:T_final) * 0, 0.1),
  #                    rnorm(T_final, -5 + (1:T_final) * 0.5, 0.1),
  #                    rnorm(T_final, -1 + (1:T_final) * (-0.3), 0.1)),
  #              
  #              cbind(rnorm(T_final, 0 + (1:T_final) * (-0.1), 0.1),
  #                    rnorm(T_final, -1 + (1:T_final) * 0, 0.1),
  #                    rnorm(T_final, -1.5 + (1:T_final) * (-0.05), 0.1),
  #                    rnorm(T_final, -1.8 + (1:T_final) * 0, 0.1),
  #                    rnorm(T_final, -2.5 + (1:T_final) * 0, 0.1)),
  #              
  #              cbind(rnorm(T_final, 0.6 + (1:T_final) * (-0.05), 0.05),
  #                    rnorm(T_final, -0.2 + (1:T_final) * 0, 0.05),
  #                    rnorm(T_final, -0.6 + (1:T_final) * 0.1, 0.05),
  #                    rnorm(T_final, -0.6 + (1:T_final) * 0, 0.05),
  #                    rnorm(T_final, -1 + (1:T_final) * (0.02), 0.05)))
  
  beta <- list(cbind(rnorm(T_final, -2 + (1:T_final) * (-0.02), 0.05),
                     rnorm(T_final, -2.6 + (1:T_final) * (-0.02), 0.05),
                     rnorm(T_final, -3.2 + (1:T_final) * (-0.02), 0.05),
                     rnorm(T_final, -3.8 + (1:T_final) * (-0.02), 0.05),
                     rnorm(T_final, -4.4 + (1:T_final) * (-0.02), 0.05)),
               
               cbind(rnorm(T_final, -4.5 + (1:T_final) * (-0.02), 0.1),
                     rnorm(T_final, -5.1 + (1:T_final) * (-0.02), 0.1),
                     rnorm(T_final, -5.8 + (1:T_final) * (-0.02), 0.1),
                     rnorm(T_final, -6.5 + (1:T_final) * (-0.02), 0.1),
                     rnorm(T_final, -7.1 + (1:T_final) * (-0.02), 0.1)),
               
               cbind(rnorm(T_final, -3 + (1:T_final) * (-0.02), 0.1),
                     rnorm(T_final, -3.5 + (1:T_final) * (-0.02), 0.1),
                     rnorm(T_final, -4 + (1:T_final) * (-0.02), 0.1),
                     rnorm(T_final, -4.5 + (1:T_final) * (-0.02), 0.1),
                     rnorm(T_final, -5 + (1:T_final) * (-0.02), 0.1)),
               
               cbind(rnorm(T_final, -4 + (1:T_final) * (-0.02), 0.1),
                     rnorm(T_final, -4.3 + (1:T_final) * (-0.02), 0.1),
                     rnorm(T_final, -4.6 + (1:T_final) * (-0.02), 0.1),
                     rnorm(T_final, -4.9 + (1:T_final) * (-0.02), 0.1),
                     rnorm(T_final, -5.2 + (1:T_final) * (-0.02), 0.1)),
               
               cbind(rnorm(T_final, -1 + (1:T_final) * (-0.02), 0.01),
                     rnorm(T_final, -1.3 + (1:T_final) * (-0.02), 0.01),
                     rnorm(T_final, -1.7 + (1:T_final) * (-0.02), 0.01),
                     rnorm(T_final, -2.1 + (1:T_final) * (-0.02), 0.01),
                     rnorm(T_final, -2.5 + (1:T_final) * (-0.02), 0.01)),
               
               cbind(rnorm(T_final, 0.5 + (1:T_final) * (-0.02), 0.02),
                     rnorm(T_final, 0 + (1:T_final) * (-0.02), 0.02),
                     rnorm(T_final, -0.5 + (1:T_final) * (-0.02), 0.02),
                     rnorm(T_final, -1 + (1:T_final) * (-0.02), 0.02),
                     rnorm(T_final, -1.5 + (1:T_final) * (-0.02), 0.02)))
                    
  
  labs <- list(matrix(c(rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6)),
                      nrow = n,
                      ncol = T_final,
                      byrow = T),
               matrix(c(rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6)),
                      nrow = n,
                      ncol = T_final,
                      byrow = T),
               matrix(c(rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6)),
                      nrow = n,
                      ncol = T_final,
                      byrow = T),
               matrix(c(rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6)),
                      nrow = n,
                      ncol = T_final,
                      byrow = T),
               matrix(c(rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6)),
                      nrow = n,
                      ncol = T_final,
                      byrow = T),
               matrix(c(rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6),
                        rep(sample(5, replace = T), each = 6)),
                      nrow = n,
                      ncol = T_final,
                      byrow = T)
  )

  
  
  # Beta divided for countries
  beta_c1 <- sapply(1:ncol(S), function(j) beta[[j]][cbind(1:30, labs[[j]][1, ])])
  beta_c2 <- sapply(1:ncol(S), function(j) beta[[j]][cbind(1:30, labs[[j]][2, ])])
  beta_c3 <- sapply(1:ncol(S), function(j) beta[[j]][cbind(1:30, labs[[j]][3, ])])
  beta_c4 <- sapply(1:ncol(S), function(j) beta[[j]][cbind(1:30, labs[[j]][4, ])])
  beta_c5 <- sapply(1:ncol(S), function(j) beta[[j]][cbind(1:30, labs[[j]][5, ])])
  
  # Data
  data1 <- sapply(1:T_final, function(t) (beta_c1 %*% t(S))[t, ]) + 
    rnorm(length(ages) * T_final, mean = 0, sd = sigma_vec[1])
  data2 <- sapply(1:T_final, function(t) (beta_c2 %*% t(S))[t, ]) + 
    rnorm(length(ages) * T_final, mean = 0, sd = sigma_vec[2])
  data3 <- sapply(1:T_final, function(t) (beta_c3 %*% t(S))[t, ]) + 
    rnorm(length(ages) * T_final, mean = 0, sd = sigma_vec[3])
  data4 <- sapply(1:T_final, function(t) (beta_c4 %*% t(S))[t, ]) + 
    rnorm(length(ages) * T_final, mean = 0, sd = sigma_vec[4])
  data5 <- sapply(1:T_final, function(t) (beta_c5 %*% t(S))[t, ]) + 
    rnorm(length(ages) * T_final, mean = 0, sd = sigma_vec[5])
  
  out <- list(data = list(data1, data2, data3, data4, data5),
              beta = beta,
              labs = labs)
}


n <- 5
T_final <- 30
ages <- 0:100
ages <- 0:100
knots <- c(20, 40)
randomseed <- 6
sigma_vec <- c(0.05, 0.1, 0.1, 0.15, 0.25)

sim2 <- simulateData(n, T_final, ages, knots, sigma_vec, randomseed)

par(mfrow = c(2, 3))
for (i in 1:5){
  matplot(sim2$data[[i]], type = "l", ylim = c(-8, -2))
}

# save.image("simdata/simdata4.RData")


# Verifica che i dati simulati siano sensati
# Devo caricare il dataset mortality
data_list_man <- list(ita_man = log(Y_ita_man / N_ita_man),
                      swe_man = log(Y_swe_man / N_swe_man),
                      uk_man = log(Y_uk_man / N_uk_man),
                      us_man = log(Y_us_man / N_us_man))
Y <- lapply(data_list_man, t)

windows(); par(mfrow = c(2, 3))
for (i in 1:5){
  matplot(sim2$data[[i]], type = "l",
          xlab = "age", ylab = "log-mort. rate",
          main = paste("Unit", i))
}


windows(); par(mfrow = c(2, 3))
for (i in 1:5){
  matplot(Y[[i]][, 1:30], type = "l")
}




windows(); par(mfrow = c(2, 3))
for (i in 1:5){
  matplot(t(sim2$data[[i]])[,20:40], type = "l")
}


windows(); par(mfrow = c(2, 3))
for (i in 1:5){
  matplot(t(Y[[i]])[1:30, 20:40], type = "l")
}


