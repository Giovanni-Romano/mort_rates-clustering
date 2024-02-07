rm(list=ls())
library(splines2)

setwd("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering")
load("data/simdata_no_hierarchy/simulated_data_5countries.Rdata")
source("utils/funzioni.R")


set.seed(4238)


data <- out$data
str(data)





### ### ### ### ### ###
#### Basi B-Spline ####
### ### ### ### ### ###
# Maximum age
ages <- 0:100
ages_no0 <- ages[-1]

# Construction of the spline basis functions
S.temp <- bSpline(ages_no0, 
                  knots = out$knots, 
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
S <- rbind(c(1, rep(0, k)),
           cbind(rep(0, length(ages_no0)), S.temp)
)
colnames(S) <- NULL





### ### ### ### ### ### ### ### ### ### ###
#### Definizione parametri del modello ####
### ### ### ### ### ### ### ### ### ### ###

# Numero istanti temporali
T_final <- ncol(data[[1]]); T_final
# Numero stati considerati
n <- length(data); n
# Dimensione stato latente (-> numero coeff)
p <- ncol(S); p
# beta_sigma <- 100; alpha_sigma <- beta_sigma/est_sigma^2 + 1
alpha_sigma <- 1; beta_sigma <- 0.01
beta_sigma / (alpha_sigma - 1)
beta_sigma^2 / (alpha_sigma - 1)^2 / (alpha_sigma - 2)

# Hyperparams for the Beta prior on alpha
a_alpha <- 1; b_alpha <- 1
# Parametro di concentrazione del CRP
M <- 0.5
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
n_iter <- 10000





### ### ### ### ### ### ### ### ### ### ### ### ### ### #
#### Estimate of mean-level of GP through regression ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### #
# We want to center the GPs for the beta's around a mean level estimated as 
#   the regression of the data on the B-spline basis. We do this for each year 
#   and then we smooth the global trend through loess.
GP_means_tmp <- matrix(NA, nrow = p, ncol = T_final)
for (t in 1:T_final){
  
  cat(t, "")
  data_tmp <- sapply(data, function(x) x[, t])
  
  melt_tmp <- reshape2::melt(data_tmp)
  
  df <- data.frame(Y = melt_tmp$value, matrix(t(S), 
                                              nrow = nrow(melt_tmp), ncol = p, 
                                              byrow = T))
  mod <- lm(Y ~ -1 + ., data = df)
  # summary(mod)
  
  GP_means_tmp[, t] <- mod$coefficients
  
  # if (t %% 20 == 1) {
  #   windows(); par(mfrow = c(4, 5))
  # }
  # matplot(data_tmp, type = "p", pch = 1)
  # lines(1:99, mod$coefficients %*% t(S), col = "gold", lwd = 2)
}

GP_means <- list()
for (j in 1:p){
  df <- data.frame(x = 1:T_final, y = GP_means_tmp[j, ])
  GP_means[[j]] <- loess(y ~ x, data = df)$fitted
}


# delta_j
alpha_delta <- 1; beta_delta <- 0.01
beta_delta / (alpha_delta - 1)
beta_delta^2 / (alpha_delta - 1)^2 / (alpha_delta - 2)



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


delta_res <- replicate(p,
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

# Temp objects
gamma_temp <- labels_temp <-
  beta_temp <- replicate(p,
                         array(NA,
                               dim = c(n, # numero di osservazioni
                                       ncol = T_final # numero istanti temporali
                               )
                         ),
                         simplify = FALSE)

delta_temp <- replicate(p, list(NA))
alpha_temp <- replicate(p, list(NA))
sigma_temp <- array(NA, dim = c(n))





### ### ### ### ### ###
### Inizializzazioni ##
### ### ### ### ### ###
sigma_res[ , 1] <- sigma_temp <- sqrt(1/rgamma(n, shape = alpha_sigma, rate = beta_sigma))
for (j in 1:p){
  
  cat(j, "")
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### #
  ### Here I initialize delta_j and phi_j as the constant value fixed above. ##
  ### Then in the for-loop I will NOT update them.                          ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### #
  delta_res[[j]][1] <- delta_temp[[j]] <- sqrt(1/rgamma(1, shape = alpha_delta, rate = beta_delta))
  
  alpha_res[[j]][1] <- alpha_temp[[j]] <- rbeta(1, a_alpha, b_alpha)
  
  # Per ora inizializzo tutti i gamma = 0, così da lasciare piena libertà di 
  #   movimento alla prima iterazione. Poi posso pensare di inizializzare anche
  #   loro dalla prior.
  gamma_res[[j]][ , , 1] <- gamma_temp[[j]][ , ] <- 0
  
  for (t in 1:T_final){
    lab <- rho2lab(rCRP(n, M))
    labels_temp[[j]][, t] <- labels_res[[j]][, t, 1] <- lab
  }
  beta_res[[j]][ , , 1] <- beta_temp[[j]] <- 
    rmvn(n,
         mu = GP_means[[j]],
         sigma = delta_temp[[j]] * cholSIG,
         isChol = TRUE)
}

rm(j); rm(t)




### ### ### ### ### ###
### GIBBS SAMPLING ####
### ### ### ### ### ###
Y <- lapply(data, function(y) t(y))

inizio <- last_time <- Sys.time()

for (d in 2:n_iter){ # Ciclo sulle iterazioni
  
  for (j in 1:p){ # Ciclo sui coefficienti
    
    b_j <- beta_temp[[j]]
    
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
                                              lab_tm1 = labels_temp[[j]][, t-1])
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
      for (i in 1:n){
        
        # Non assegno qua anche a labels_res perché dovrò prima sistemarlo
        labels_temp[[j]][i, t] <-  up_label_i(i = i,
                                              j = j,
                                              Y_it = Y[[i]][t, ],
                                              beta_i = beta_units[i, ],
                                              beta_cluster = b_j[, t], 
                                              sigma_i = sigma_temp[i],
                                              lab_t = labels_temp[[j]][,t],
                                              lab_tp1 = labels_temp[[j]][,t+1],
                                              gamma_tp1 = if (t == T_final) {'last time'}
                                              else {gamma_temp[[j]][, t+1]},
                                              spline_basis = S)
        
        
        # In this version I must not reorder labels to avoid "gaps" anymore 
        #   because now the value of the label is important because it 
        #   corresponds to a particular curve of (beta_kj1, ..., beta_kjT).
        
      } # Fine ciclo sulle osservazioni per labels
      
      
    } # END OF FOR LOOP OVER TIME ISTANTS "t"
    
    
    # Update on beta's is now joint, so I have to move also the updates of
    #   theta and tau because they must happen after the one of the beta's
    
    ### ### ### ### ### ### ###
    ### ### UPDATE BETA ### ###
    ### ### ### ### ### ### ###
    # I get the beta for every coeff, for every country and for every time
    # It's an array n*T_final*p
    beta_units <- vapply(1:p, 
                         function(h) vapply(1:T_final, 
                                            function(s) beta_temp[[h]][labels_temp[[h]][ , s], s],
                                            FUN.VALUE = double(n)),
                         FUN.VALUE = matrix(1, n, T_final))
    
    # I obtain an array n*X*T to center the Y with respect to all the splines
    #   but the j-th
    centering <- simplify2array(apply(beta_units, 1, function(B) B[, -j]%*%t(S[, -j]), simplify = FALSE))
    # I change the str() of Y to center it simply
    Y_simpl <- simplify2array(Y)
    Y.tilde <- Y_simpl - centering
    
    # Sum_k g_j(x)^2
    sum_sq_g_j <- sum(S[, j]^2)
    
    # I get the labels for j-th spline because I'll use them multiple times
    lab_j <- labels_temp[[j]]
    
    for (k in 1:n){ # Inizio ciclo sui cluster per i beta
      
      # Sum of the precision that contributes to the varcov matrix of the
      #   of the posterior in the part that comes from the likelihood
      sum_prec_k <- vapply(1:T_final, 
                           function(s) sum(1/sigma_temp[which(lab_j[, s] == k)]^2), 
                           FUN.VALUE = double(1))
      # Sum_{i : c_ijt = k} \tilde{y_ixt} / sigma_i^2 --> goes in the mean of 
      #   the posterior in the part that comes from the likelihood
      sum_y_sig_k <- vapply(1:T_final, 
                            function(s) colSums( t(Y.tilde[s, , which(lab_j[, s] == k)])/
                                                   sigma_temp[which(lab_j[, s] == k)]^2 ), 
                            FUN.VALUE = double(length(ages)))
      # Sum_x g_j(x)*"term just computed" --> contribution of the likelihood 
      #   to the posterior mean
      lik2mean <- S[, j] %*% sum_y_sig_k
      # I can't compute this outside the "k" loop because sum_y_sig depends on k
      
      # Contribution of the likelihood to the varcov (precision) matrix;
      #   it's the argument of the diagonal matrix
      lik2var <- sum_prec_k * sum_sq_g_j
      
      # What observations are in the k-th cluster?
      #   It's a list of 87 elements, each with the observations in the k-th
      #   cluster. I cannot simplify to a matrix because the elements have
      #   different lengths.
      obs_k <- apply(lab_j, 2, function(l) which(l == k))
      
      
      beta_temp[[j]][k, ] <- up_beta(Y.t = Y.tilde,
                                     l2v = lik2var,
                                     l2m = lik2mean,
                                     obs = obs_k,
                                     SIGinv = SIGMAinv,
                                     d_j = delta_temp[[j]],
                                     phi_j = GP_means[[j]])
    } # Fine ciclo sui cluster per i beta
    
    
    # Aggiungo gli elementi anche negli oggetti non-temp
    labels_res[[j]][ , , d] <- as.integer(labels_temp[[j]])
    gamma_res[[j]][ , , d] <- gamma_temp[[j]]
    beta_res[[j]][ , , d] <- beta_temp[[j]]
    
    
    
    ### ### ### ### ####
    ### UPDATE DELTA ###
    delta_res[[j]][d] <- delta_temp[[j]] <-
      sqrt(1/rgamma(1, 
                    shape = alpha_delta + 0.5 * T_final * n,
                    rate = beta_delta + 
                      0.5 * sum( (t(cholSIGinv) %*% 
                                    (t(beta_temp[[j]]) - GP_means[[j]]))^2 ))
      )
    # REMEMBER THAT R GIVES THE UPPER TRIANGULAR PART OF THE CHOLESKY DECOMPOSITION
    
    
    
    ### ### ### ### ####
    ### UPDATE ALPHA ###
    alpha_res[[j]][d] <- alpha_temp[[j]] <-
      up_alpha_j(sum(gamma_temp[[j]]),
                 n * T_final,
                 a_alpha,
                 b_alpha)
    
    
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
        sqrt(1/rgamma(1, shape = alpha_sigma + 0.5 * T_final * length(ages),
                      rate = beta_sigma + 0.5 * sum((Y[[i]] - t(means[i, , ]))^2))
        )
    }
  } # END OF FOR LOOP OVER COEFFICIENTS "j"
  
  
  if ((d %% 50) == 1) {
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

# save.image("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/res/dynamic_on_beta/no_hierarchy/sim_study/sim_study_full_bayes/res.RData")
