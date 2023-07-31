rm(list=ls())

library(drpm)

load("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/BSP_Pavone/output/mortality.Rdata")

Y <- rbind(ita_man = log(Y_ita_man / N_ita_man)[1:87, 1],
           swe_man = log(Y_swe_man / N_swe_man)[1:87, 1],
           uk_man = log(Y_uk_man / N_uk_man)[1:87, 1],
           us_man = log(Y_us_man / N_us_man)[1:87, 1])
rownames(Y) <- substr(rownames(Y), 1, 2)

set.seed(4238)



### ### ### ### ### ### ### ### ### ### ###
#### Definizione parametri del modello ####
### ### ### ### ### ### ### ### ### ### ###

# Numero istanti temporali
T_final <- ncol(Y); T_final
# Numero stati considerati
n <- nrow(Y); n
# Iperparametri della prior dell'ultimo layer
m0 <- 0; s02 <- 100
# Uniforms' hyperparameters
#   Following Page's idea in paragraph 2.5 I choose these hyperparams
A_tau <- A_delta <- 10
A_sigma <- 5
# Hyperparams for the Beta prior on alpha
a_alpha <- b_alpha <- 1
# Parametro di concentrazione del CRP
M <- 2



### ### ### ### ### ### ### ### ### ###
### Definizione parametri del Gibbs ###
### ### ### ### ### ### ### ### ### ###

# Numero iterazioni
n_iter <- 5000
# RW step sizes
#   I've done some diagnostics and it seems that the best thing to do is to 
#   to pick a large stepsize, so I set it equal to the size of the domain
eps_tau <- A_tau
eps_delta <- A_delta
eps_sigma <- A_sigma





draws <- drpm_fit(y = Y,
                  s_coords=NULL,
                  M = M,
                  initial_partition = NULL,
                  starting_alpha = 0,
                  unit_specific_alpha = FALSE,
                  time_specific_alpha = FALSE,
                  alpha_0 = FALSE,
                  eta1_0 = FALSE,
                  phi1_0 = FALSE,
                  modelPriors = c(m0, s02, A_sigma, A_tau, A_delta, 1),
                  alphaPriors = rbind(c(a_alpha, b_alpha)),
                  simpleModel = 0,
                  mh = c(eps_sigma, eps_tau, eps_delta, 0.1, 0.1),
                  verbose = TRUE,
                  draws = 5000, burn = 1000, thin = 1)



# Avg number of clusters
windows()
plot(apply(draws$Si, 1, function(x) mean(apply(x, 2, function(y) max(y)))), 
     type = "l", 
     ylim = c(0, 4),
     xlab = "time",
     ylab = "avg # of cluster",
     main = "Coeff of infant")


# Estimated partition
library(salso)
system.time(est_clust <- t(apply(draws$Si, 2,
                               function(y)
                                 salso(x = t(y),
                                       loss = binder(a = NULL),
                                       maxNClusters = 4,
                                       maxZealousAttempts = 10,
                                       nRuns = 11,
                                       nCores = 11
                                 )
))
)

rownames(est_clust) <- dimnames(Y)[[1]]; colnames(est_clust) <- dimnames(Y)[[2]]
est_clust


# Prob of co-clustering
same_cluster0 <- matrix(NA, 6, 87)

up3 <- upper.tri(matrix(1, 4, 4))

same_cluster0 <- sapply(1:T_final, 
                        function(t) rowMeans(apply(draws$Si[t, , ], 2, 
                                                   function(x) outer(x, x, "==")[up3]))
)

rownames(same_cluster0) <- 
  outer(rownames(Y), rownames(Y), paste, sep = ":")[up3]
colnames(same_cluster0) <- 1933:2019

sc_melt0 <- melt(same_cluster0)
colnames(sc_melt0) <- c("Countries", "Year", "Value")

pl0 <- ggplot(sc_melt0,
              aes(x = Year, y = Value, 
                  group = Countries, 
                  colour = Countries, lty = Countries)) + 
  geom_line(lwd = 1) + 
  ylim(0, 1) +
  labs(title = "MonteCarlo prob. of co-clustering (C code - complete model)", 
       x = "Year", y = "")
windows()
pl0


windows()
par(mfrow = c(1, 2))
n_iter <- 4000
plot(1:87, 
     sapply(1:87, function(t) mean(draws$mu[t, , ][cbind(draws$Si[t, 1, ], 1:n_iter)])), 
     type = 'l',
     ylab = "", xlab = "",
     main = "Beta of each country (C code)",
     ylim = c(-6, -1.5))
lines(1:87, 
      sapply(1:87, function(t) mean(draws$mu[t, , ][cbind(draws$Si[t, 2, ], 1:n_iter)])), 
      type = 'l',
      ylab = "", xlab = "",
      col = 2)
lines(1:87, 
      sapply(1:87, function(t) mean(draws$mu[t, , ][cbind(draws$Si[t, 3, ], 1:n_iter)])), 
      type = 'l',
      ylab = "", xlab = "",
      col = 3)
lines(1:87, 
      sapply(1:87, function(t) mean(draws$mu[t, , ][cbind(draws$Si[t, 4, ], 1:n_iter)])), 
      type = 'l',
      ylab = "", xlab = "",
      col = 4)
legend("topright", col = 1:4, lty = 1, legend = rownames(Y))


plot(1:87, 
     Y[1, ], 
     type = 'l',
     ylab = "", xlab = "",
     main = "Data of each country",
     ylim = c(-6, -1.5))
lines(1:87, 
      Y[2, ],
      type = 'l',
      ylab = "", xlab = "",
      col = 2)
lines(1:87, 
      Y[3, ],
      type = 'l',
      ylab = "", xlab = "",
      col = 3)
lines(1:87, 
      Y[4, ],
      type = 'l',
      ylab = "", xlab = "",
      col = 4)
legend("topright", col = 1:4, lty = 1, legend = rownames(Y))


