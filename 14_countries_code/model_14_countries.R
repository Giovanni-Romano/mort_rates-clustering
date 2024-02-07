rm(list=ls())
library(splines2)

# setwd("/mnt/c/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering")
setwd("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering")
load("data/death_rates_1x1/14_countries_preproc.RData")
source("utils/funzioni_no_hier_dyn_on_theta.R")


set.seed(4238)

Y <- lapply(rates_male, log)


### ### ### ### ### ###
#### Basi B-Spline ####
### ### ### ### ### ###
# Maximum age
Z <- ncol(rates_male[[1]])
ages <- as.numeric(colnames(rates_male[[1]]))
ages_no0 <- ages[-1]

# Construction of the spline basis functions
S.temp <- bSpline(ages_no0, 
                 knots = c(seq(5, 40, by = 5), 50, 60, seq(70, 95, by = 5)), 
                 degree = 2, 
                 intercept = TRUE)
# Il numero di basi B-spline è pari al numero di nodi interni più l'ordine
#   scelto per le basi, che è pari al degree più uno. Quindi in questo caso
#   il numero di basi B-spline è 16+3=19.
K <- ncol(S.temp)

# Normalizzazione delle basi (perché Federico lo fa?)
for (k in 1:K){
  S.temp[,k] <- S.temp[,k]/max(S.temp[,k])
}
# In più aggiungo una cornice a L come prima riga e prima colonna fatta tutta
#   di 0 tranne l'elemento in posizione (1, 1) che è 1. In questo modo
#   aggiungo anche l'età 0 e la base attiva solo in 0.
S <- rbind(c(1, rep(0, K)),
           cbind(rep(0, length(ages_no0)), S.temp)
           )
colnames(S) <- NULL





### ### ### ### ### ### ### ### ### ### ###
#### Definizione parametri del modello ####
### ### ### ### ### ### ### ### ### ### ###

# Numero istanti temporali
T_final <- nrow(Y[[1]]); T_final
# Numero stati considerati
n <- length(Y); n
# Dimensione stato latente (-> numero coeff)
p <- ncol(S); p
# Hyperparameters prior on sigma_i^2 (observations' variance)
a_sigma <- 0.001; b_sigma <- 0.001
# Hyperparams for the Beta prior on alpha
a_alpha <- 1; b_alpha <- 1
# Parametro di concentrazione del CRP
a_M <- 0.002; b_M <- 0.001
# Sigma --> parte della matrice varcov dei vettori di beta
l <- 3/2 # With l = 2 the correlations at distance (1, 2, 3, 4, 5, 6, 7) are
# (0.88 0.61 0.32 0.14 0.04 0.01 0.00);
# with l = 2 they are (0.80 0.41 0.14 0.03 0.00 0.00 0.00)
SIGMA <- outer(1:T_final, 1:T_final, function(x1, x2) sq_exp_ker(x1-x2, l = l))
cholSIG <- chol(SIGMA)
# The Cholesky decomposition is useful because mvnfast::rmvn is way faster 
#   specifying the Cholesky decomp. of the varcov matrix if it has been already
#   computed. Otherwise, computing Chol plus applying mvnfast::rmvn(..., isChol = T)
#   is as long as applying mvnfast::rmvn(..., isChol = F)
SIGMAinv <- chol2inv(cholSIG)
# The inverse of SIGMA is used in the update of the thetas, because it 
#   contributes to the posterior variance/precision
cholSIGinv <- solve(cholSIG)

# REMEMBER: Chol. decomp. on R gives the upper triangular!!!!

### ### ### ### ### ### ### ### ### ###
### Definizione parametri del Gibbs ###
### ### ### ### ### ### ### ### ### ###

# Numero iterazioni
n_iter <- 5000




### ### ### ### ### ### ### ### ### ### ### ### ### ### #
#### Estimate of mean-level of GP through regression ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### #
# We want to center the GPs for the beta's around a mean level estimated as 
#   the regression of the data on the B-spline basis. We do this for each year 
#   and then we smooth the global trend through loess.
GP_means_tmp <- matrix(NA, nrow = p, ncol = T_final)
cat("Start estimate of GP means \n\n")
for (t in 1:T_final){
  
  cat(t, "")
  data_tmp <- sapply(Y, function(x) x[t, ])
  
  melt_tmp <- reshape2::melt(data_tmp)
  
  df <- data.frame(Y = melt_tmp$value, matrix(t(S), 
                                              nrow = nrow(melt_tmp), ncol = p, 
                                              byrow = T))
  mod <- lm(Y ~ -1 + ., data = df)
  
  GP_means_tmp[, t] <- mod$coefficients
}

GP_means <- list()
for (j in 1:p){
  df <- data.frame(x = 1:T_final, y = GP_means_tmp[j, ])
  GP_means[[j]] <- loess(y ~ x, data = df)$fitted
}
cat("End estimate of GP means \n\n")

# Hyperparameters delta_j^2
a_delta <- 0.001; b_delta <- 0.001
# Hyperparameters omega_j^2
a_omega <- 0.001; b_omega <- 0.001



### ### ### ### ### ### ### ### ### ### ###
### Creo gli oggetti in cui salverò gli ###
### output dei vari step del Gibbs      ###
### ### ### ### ### ### ### ### ### ### ###

# Ho provato ad inizializzare questi oggetti con 99999 sperando di velocizzare,
#   ma non cambia niente. Perciò, preferisco tenere NA, almeno se qualcosa va
#   storto me ne accorgo.
gamma_res <- labels_res <- replicate(p,
                                     array(NA,
                                           dim = c(n, # numero di osservazioni
                                                   T_final, # numero istanti temporali
                                                   n_iter # numero iterazioni
                                           )
                                     ),
                                     simplify = FALSE)

beta_res <- replicate(p,
                      array(NA,
                            dim = c(n, # numero di osservazioni
                                    T_final, # numero istanti temporali
                                    n_iter # numero iterazioni
                            )
                      ),
                      simplify = FALSE)

phi_res <- replicate(p,
                     array(NA,
                           dim = c(T_final, # numero istanti temporali
                                   n_iter # numero iterazioni
                           )
                     ),
                     simplify = FALSE)

delta_res <- omega_res <- replicate(p,
                                    array(NA,
                                          dim = c(n_iter) # numero iterazioni
                                    ),
                                    simplify = FALSE)


alpha_res <- replicate(p,
                       array(NA,
                             dim = c(n_iter # numero iterazioni
                             )
                       ),
                       simplify = FALSE)


sigma_res <- array(NA,
                   dim = c(n,
                           n_iter # numero iterazioni
                           )
                   )
M_res <- replicate(p, 
                   array(NA, dim = c(n_iter)),
                   simplify = F)

# Temp objects
gamma_temp <- labels_temp <-
  beta_temp <- replicate(p,
                         array(NA,
                                dim = c(n, # numero di osservazioni
                                        ncol = T_final # numero istanti temporali
                                        )
                         ),
                         simplify = FALSE)
phi_temp <- replicate(p,
                     array(NA,
                           dim = c(T_final # numero istanti temporali
                                   )
                     ),
                     simplify = FALSE)

delta_temp <- omega_temp <- replicate(p, list(NA))
alpha_temp <- replicate(p, list(NA))
sigma_temp <- array(NA, dim = c(n))
M_temp <- replicate(p, list(NA))





### ### ### ### ### ###
### Inizializzazioni ##
### ### ### ### ### ###

sigma_res[ , 1] <- sigma_temp <- sapply(Y, sd)#sqrt(1/rgamma(n, shape = a_sigma, rate = b_sigma))
cat("Start inizialization of values for first iteration of Gibbs Sampler \n\n")
for (j in 1:p){
  cat(j, "")
  
  omega_res[[j]][1] <- omega_temp[[j]] <- 1 #sqrt(1/rgamma(1, shape = a_delta, rate = b_delta))
  delta_res[[j]][1] <- delta_temp[[j]] <- 1 #sqrt(1/rgamma(1, shape = a_delta, rate = b_delta))
  
  alpha_res[[j]][1] <- alpha_temp[[j]] <- rbeta(1, a_alpha, b_alpha)
  
  # Per ora inizializzo tutti i gamma = 0, così da lasciare piena libertà di 
  #   movimento alla prima iterazione. Poi posso pensare di inizializzare anche
  #   loro dalla prior.
  gamma_res[[j]][ , , 1] <- gamma_temp[[j]][ , ] <- 0
  
  phi_res[[j]][ , 1] <- phi_temp[[j]] <-
    drop(rmvn(1,
              mu = GP_means[[j]],
              sigma = omega_temp[[j]] * cholSIG,
              isChol = TRUE))
  M_temp[[j]] <- 3 #rgamma(1, a_M, b_M)
  for (t in 1:T_final){
    lab <- rho2lab(rCRP(n, M_temp[[j]]))
    labels_temp[[j]][, t] <- labels_res[[j]][, t, 1] <- lab
    
    beta_res[[j]][1:max(lab) , t, 1] <- beta_temp[[j]][1:max(lab), t] <- 
      rnorm(max(lab), mean = phi_temp[[j]][t], sd = delta_temp[[j]])
  }
  
}
cat("End nizialization of values for first iteration of Gibbs Sampler \n\n")
rm(j); rm(t)




### ### ### ### ### ###
### GIBBS SAMPLING ####
### ### ### ### ### ###

inizio <- last_time <- Sys.time()

for (d in 1:n_iter){ # Ciclo sulle iterazioni
  
  for (j in 1:p){ # Ciclo sui coefficienti
    
    for (t in 1:T_final){ # Ciclo sugli istanti
      ### ### ### ### ### ### ####
      ### ### UPDATE GAMMA ### ###
      for (i in 1:n){ # Ciclo sulle osservazioni per gamma
        # Potrei direttamente saltare il caso t==1 perché tanto
        #   gamma_temp è inizializzato pieno di 0 e la colonna
        #   del tempo t=1 non verrebbe mai toccata.
        if (t == 1){ # Se t = 1, allora i gamma vengono fissati a 0
          # perché non c'è alcun movimento
          # Assegno la nuova gamma nel tmp
          gamma_temp[[j]][i, 1] <- 0
        } else if (t > 1){
          # Assegno la nuova gamma nel tmp
          gamma_temp[[j]][i, t] <- up_gamma_i(i = i, 
                                              gamma = gamma_temp[[j]][, t], 
                                              alpha_t = alpha_temp[[j]],
                                              lab_t = labels_temp[[j]][, t],
                                              lab_tm1 = labels_temp[[j]][, t-1],
                                              M = M_temp[[j]])
        }
      } # Fine ciclo sulle osservazioni per gamma
      
      
      ### ### ### ### ### ### ### #
      ### ### UPDATE LABELS ### ###
      # I can compute beta_units outside the loop in "i" because when I update 
      #   the i-th unit its labels to use to find beta_units are those at the 
      #   iteration d-1, that are those in labels_temp.
      beta_units <- vapply(1:p, 
                            function(x) beta_temp[[x]][labels_temp[[x]][ , t], t],
                            FUN.VALUE = double(n))
      can_move <- which(gamma_temp[[j]][, t] == 0)
      for (i in can_move){
        
        newclustervalue <- rnorm(1, phi_temp[[j]][t], delta_temp[[j]])
        
        # Non assegno qua anche a labels_res perché dovrò prima sistemarlo
        labels_temp[[j]][i, t] <-  up_label_i(i = i, 
                                              j = j,
                                              Y_it = Y[[i]][t, ],
                                              beta_i = beta_units[i, ],
                                              beta_cluster = beta_temp[[j]][, t],
                                              sigma_i = sigma_temp[i],
                                              lab_t = labels_temp[[j]][,t],
                                              lab_tp1 = labels_temp[[j]][,t+1], 
                                              gamma_tp1 = if (t == T_final) {'last time'} 
                                              else {gamma_temp[[j]][, t+1]},
                                              newclustervalue = newclustervalue,
                                              spline_basis = S,
                                              M = M_temp[[j]])
        
        
        # Sistemo le labels in modo tale che non ci siano buchi (tipo non 
        # voglio c(1, 2, 2, 4) ma c(1, 2, 2, 4)) e che siano ordinate
        # (tipo non voglio c(2, 1, 2, 3) ma c(1, 2, 1, 3)).
        labels_new <- c(1L)
        counter <- 2L
        for (ii in 2:n){
          if (!(labels_temp[[j]][ii, t] %in% labels_temp[[j]][1:(ii-1), t])){
            labels_new[ii] <- counter
            counter <- counter+1
          } else{
            labels_new[ii] <- labels_new[which(labels_temp[[j]][, t] == 
                                                 labels_temp[[j]][ii, t])[1]]
          }
        }
        
        # Voglio riordinare anche beta_temp in accordo con le labels:
        #   facendo unique(...) ottengo i valori unici delle labels in ordine
        #   di apparizione nel vettore (e non in ordine cresc./decresc.).
        #   Questo vettore dà la permutazione necessaria.
        
        if (all(labels_temp[[j]][i, t] > labels_temp[[j]][-i, t])){
          # Se la label dell'unità i è maggiore di tutte le altre, vuol dire
          #   che è andato in un nuovo cluster e non in uno di quelli esistenti,
          #   quindi aggiungo il beta del nuovo cluster al vettore dei beta
          #   dei cluster.
          if (labels_temp[[j]][i, t] > n){
            # If the new cluster is the (n+1)st, then I cannot modify directly
            #   beta_temp because I would exceed its size.
            beta_extended <- c(beta_temp[[j]][, t], newclustervalue)
            sort_beta_temp <- unique(labels_temp[[j]][, t])
            beta_new <- beta_extended[sort_beta_temp]
          } else { # If new cluster is <= n
            beta_temp[[j]][labels_temp[[j]][i, t], t] <- newclustervalue
            sort_beta_temp <- unique(labels_temp[[j]][, t])
            beta_new <- beta_temp[[j]][sort_beta_temp, t]
          }
        } else { # If not new cluster
          sort_beta_temp <- unique(labels_temp[[j]][, t])
          beta_new <- beta_temp[[j]][sort_beta_temp, t]
        }
        
        beta_temp[[j]][, t] <- NA
        beta_temp[[j]][1:length(sort_beta_temp), t] <- beta_new
        
        # Finally, I assign the new fixed label to the storing object
        labels_temp[[j]][, t] <- labels_new
        
      } # Fine ciclo sulle osservazioni per labels
        
      ### ### ### ### ### ### ###
      ### ### UPDATE BETA ### ###
      ### ### ### ### ### ### ###
      n_cluster <- max(labels_temp[[j]][ , t])
      
      Y_t <- t(vapply(Y, 
                      function(y) y[t, ], 
                      FUN.VALUE = double(length(ages))))
      
      beta_t <- vapply(beta_temp, 
                       function(b) b[, t], 
                       FUN.VALUE = double(n))
      
      labels_t <- vapply(labels_temp, 
                         function(lab) lab[, t], 
                         FUN.VALUE = double(n))
      
      for (k in 1:n_cluster){ # Inizio ciclo sui cluster per i beta
        beta_temp[[j]][k, t] <- up_beta(j = j, 
                                        k = k, 
                                        t = t,
                                        Y_t = Y_t,
                                        beta_t =  beta_t,
                                        phi_jt = phi_temp[[j]][t], 
                                        delta_j = delta_temp[[j]],
                                        sigma_vec = sigma_temp,
                                        labels_t = labels_t,
                                        spline_basis = S)
        beta_t[k] <- beta_temp[[j]][k, t]
      } # Fine ciclo sui cluster per i beta
          
    } # END OF FOR LOOP OVER TIME ISTANTS "t"
      
      
    
    # Aggiungo gli elementi anche negli oggetti non-temp
    labels_res[[j]][ , , d] <- as.integer(labels_temp[[j]])
    gamma_res[[j]][ , , d] <- gamma_temp[[j]]
    beta_res[[j]][ , , d] <- beta_temp[[j]]
    
    
    ### ### ### ### ####
    ### UPDATE DELTA ###
    to_up_deltaj <- (t(beta_temp[[j]]) - GP_means[[j]])[!is.na(t(beta_temp[[j]]))]
    delta_res[[j]][d] <- delta_temp[[j]] <-
      sqrt(1/rgamma(1, 
                    shape = a_delta + 0.5 * length(to_up_deltaj),
                    rate = b_delta + 
                      0.5 * sum(to_up_deltaj^2)))
    
    
    ### ### ### ### ##
    ### UPDATE PHI ###
    ### ### ### ### ##
    prec.lik_phi <- diag( apply(beta_temp[[j]], 2, function(x) sum(!is.na(x))) / 
                            delta_temp[[j]]^2,
                          nrow = T_final)
    var.post_phi <- mysolve( SIGMAinv / omega_temp[[j]]^2 + prec.lik_phi)
    phi_res[[j]][ ,d] <- phi_temp[[j]] <-
      drop(rmvn(1,
                mu = var.post_phi %*% ( (SIGMAinv / omega_temp[[j]]^2) %*% GP_means[[j]] + 
                                          colSums(beta_temp[[j]]/delta_temp[[j]]^2, na.rm = T) ),
                sigma = var.post_phi))
      
    ### ### ### ### ###
    ### UPDATE OMEGA ##
    ### ### ### ### ###
    omega_res[[j]][d] <- omega_temp[[j]] <-
      sqrt(1/rgamma(1, 
                    shape = a_omega + 0.5 * T_final,
                    rate = b_omega + 
                      0.5 * sum( (t(cholSIGinv) %*% 
                                    (phi_temp[[j]] - GP_means[[j]]))^2 ))
      )

    ### ### ### ### ####
    ### UPDATE ALPHA ###
    alpha_res[[j]][d] <- alpha_temp[[j]] <-
      up_alpha_j(sum(gamma_temp[[j]]),
                 n * T_final,
                 a_alpha,
                 b_alpha)
    
    
    ### ### ### ## #
    ### UPDATE M ###
    # Update Dirichlet Process concentration parameter with gamma prior
    # Escobar 1995
    eta_temp <- rbeta(1, M_temp[[j]] + 1, n)
    k_temp <- length(unique(labels_temp[[1]][, 1]))
    pi_temp <- (a_M + k_temp - 1) / ( a_M + + k_temp - 1 + n*(b_M - log(eta_temp)) ) 
    coin_temp <- (runif(1, 0, 1) < pi_temp)
    M_res[[j]][d] <- M_temp[[j]] <- 
      ifelse(coin_temp == 1, rgamma(1, a_M + k_temp, b_M - log(eta_temp)),
             rgamma(1, a_M + k_temp - 1, b_M - log(eta_temp)))
    
  } # END OF FOR LOOP OVER COEFFICIENTS "j"
  
  ### ### ### ### ### ### ### ##
  ### ### UPDATE SIGMA_i ### ###
  means <- vapply(1:T_final,
                  function(t) vapply(1:p, 
                                     function(x) beta_temp[[x]][labels_temp[[x]][ , t], t],
                                     FUN.VALUE = double(n)) %*% t(S),
                  FUN.VALUE = matrix(1, n, length(ages))
  )
  for (i in 1:n){
    sigma_res[i, d] <- sigma_temp[i] <- 
      sqrt(1/rgamma(1, shape = a_sigma + 0.5 * T_final * length(ages),
                    rate = b_sigma + 0.5 * sum((Y[[i]] - t(means[i, , ]))^2))
      )
  }
  
  if ((d %% 10) == 1) {
    now <- Sys.time()
    cat("iter. ", d-1, " fatta",
        "\n tempo passato da inizio: ", difftime(now, inizio, units = "mins"),
        "\n tempo passato da ultima iter: ", difftime(now, last_time, units = "mins"), " mins", 
        "\n ora print: ", format(Sys.time(), "%H:%M:%S"), 
        "\n\n",
        sep = "")
    last_time <- now
  }
  
} # END OF FOR LOOP OVER ITERATIONS "d"

fine <- Sys.time()

exec_time <- difftime(fine, inizio)

# save.image("/mnt/c/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/res/dynamic_on_beta/no_hierarchy/14_countries/full_bayes_dyn_on_theta/res_incomplete.RData")
save.image("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/res/dynamic_on_beta/no_hierarchy/14_countries/full_bayes_dyn_on_theta/res.RData")
