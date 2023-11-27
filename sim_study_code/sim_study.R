rm(list=ls())
library(splines2)

setwd("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering")
load("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/data/simdata_PPMx/simulated_data.Rdata")
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





### ### ### ### ### ### ### ### ### ### ###
#### Definizione parametri del modello ####
### ### ### ### ### ### ### ### ### ### ###

# Numero istanti temporali
T_final <- out$T_final; T_final
# Numero stati considerati
n <- out$n; n
# Numero coeff
p <- ncol(S); p
# Varianza delle osservazioni fissata
sigma <- out$sigma_vec; sigma
# Iperparametri della prior dell'ultimo layer
m0 <- -4; s02 <- 1^2
# Uniforms' hyperparameters
#   Following Page's idea in paragraph 2.5 I choose these hyperparams
A_delta <- 4 # It was 2, but from traceplots I notice it's too small for v02 = 0.5
A_xi <- 5
# Hyperparams for the Beta prior on alpha
a_alpha <- 2; b_alpha <- 3
# Parametro di concentrazione del CRP
M <- 2
# Sigma --> parte della matrice varcov dei vettori di beta
l <- 3/2
SIGMA <- outer(1:T_final, 1:T_final, function(x1, x2) sq_exp_ker(x1-x2, l = l))
cholSIG <- chol(SIGMA)
SIGMAinv <- chol2inv(cholSIG)

# Matrix of covariates
# X <- abind::abind(appl(y(out$x1, 2, scale), apply(out$x2, 2, scale), 
#                   rev.along = 0)
X <- abind::abind((out$x1 - mean(out$x1))/sd(out$x1), 
                  (out$x2 - mean(out$x2))/sd(out$x2),
                  rev.along = 0)
dimnames(X)[[3]] <- c("x1", "x2")

# g3() - similarity function - parameters
# For the moment I use the same value for all the covariates, but maybe 
#   it could be better to use different values for different covariates.

# Mean of the marginal likelihood
mu0 <- 0 # Because we standardize the covariates.

# Variance of the distribution that we use to marginalize the likelihood
v02 <- 1 # In Muller&Quintana&Rosner2011 they set v02 as 10 times the 
          #   empirical variance of the covariates, that we standardize.
          #   So, in our case, v02 = 10*1 = 10.
# Hyperpar. of the prior on the variance of the likelihood, that is the 
#   right extreme of a Uniform distribution.
A_s2 <- 10 # I don't know how to set it; I start with a value (10) that should
            # be enough (probably too much).




### ### ### ### ### ### ### ### ### ###
### Definizione parametri del Gibbs ###
### ### ### ### ### ### ### ### ### ###

# Numero iterazioni
n_iter <- 5000

# RW step sizes
#   I've done some diagnostics and it seems that the best thing to do is to 
#   to pick a large stepsize, so I set it equal to the size of the domain
eps_xi <- A_xi
# I noticed that the beahviour of delta in the model w/ time-dependence in the
#   beta's needed a different treatment with respect to the other variance's
#   parameters. A Uniform RW-MH with step-size equal to the length of the 
#   domain (eps = A_delta) was not okay because the likelihood is quite peaked
#   on an interval of length 0.2/0.3 and so the majority of proposed values 
#   were rejected. The solution that I am using is to reduce the eps and to 
#   use a Gaussian RW to have more proposed values close to the current point
#   but w/o completely avoid longer steps (as it would happen with a Unif RW
#   with eps = 0.1).
eps_delta <- 0.1

eps_s2 <- A_s2/4 # W/ A_s2 = 10, I think that eps_s2 = A_s2 could be too large





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


phi_res <- delta_res <- replicate(p,
                                  array(NA,
                                        dim = c(n_iter # numero iterazioni
                                        )
                                  ),
                                  simplify = FALSE)

alpha_res <- replicate(p,
                       array(NA,
                             dim = c(n_iter # numero iterazioni
                             )
                       ),
                       simplify = FALSE)

lambda_res <- xi_res <- rep(NA, n_iter)
s2_res <- replicate(p,
                    array(NA, dim = c(n_iter, 2)),
                    simplify = FALSE)


# Temp objects
gamma_temp <- labels_temp <-
  beta_temp <- replicate(p,
                         array(NA,
                               dim = c(n, # numero di osservazioni
                                       ncol = T_final # numero istanti temporali
                               )
                         ),
                         simplify = FALSE)

phi_temp <- delta_temp <- replicate(p, list(NA))
alpha_temp <- replicate(p, list(NA))
lambda_temp <- xi_temp <- NA
s2_temp <- replicate(p, 
                     rep(NA, 2),
                     simplify = FALSE)


### ### ### ### ### ###
### Inizializzazioni ##
### ### ### ### ### ###
xi_res[1] <- xi_temp <- 
  runif(1, 0, A_xi)
lambda_res[1] <- lambda_temp <- 
  rnorm(1, mean = m0, sd = sqrt(s02))

for (j in 1:p){
  s2_temp[[j]] <- s2_res[[j]][1, ] <- runif(2, min = 0, max = A_s2)
  
  delta_res[[j]][1] <- delta_temp[[j]] <- 
    runif(1, 0.1, A_delta)
  phi_res[[j]][1] <- phi_temp[[j]] <- 
    rnorm(1, lambda_temp, xi_temp)
  
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
         mu = rep(phi_temp[[j]], T_final), 
         sigma = cholSIG,
         isChol = TRUE)
}

rm(j); rm(t)




### ### ### ### ### ###
### GIBBS SAMPLING ####
### ### ### ### ### ###
Y <- lapply(data, t)
# debug(up_s2.RWM)
inizio <- last_time <- Sys.time()

for (d in 2:n_iter){ # Ciclo sulle iterazioni
  
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
                                              X_t = X[ , t, ],
                                              s2_vec = s2_temp[[j]],
                                              mu0_vec = rep(mu0, 2),
                                              v02_vec = rep(v02, 2))
        }
      } # Fine ciclo sulle osservazioni per gamma
      
      
      
      ### ### ### ### ### ### ### #
      ### ### UPDATE LABELS ### ###
      for (i in 1:n){
        beta_actual <- vapply(1:p, 
                              function(x) beta_temp[[x]][labels_temp[[x]][i, t], t],
                              FUN.VALUE = double(1))
        
        
        # Non assegno qua anche a labels_res perché dovrò prima sistemarlo
        labels_temp[[j]][i, t] <-  up_label_i(i = i, 
                                              j = j,
                                              Y_it = Y[[i]][t, ],
                                              beta_i = beta_actual,
                                              beta_cluster = beta_temp[[j]][, t],
                                              sigma_i = sigma[i],
                                              lab_t = labels_temp[[j]][,t],
                                              lab_tp1 = labels_temp[[j]][,t+1], 
                                              gamma_tp1 = if (t == T_final) {'last time'} 
                                              else {gamma_temp[[j]][, t+1]},
                                              spline_basis = S,
                                              X_t = X[ , t, ],
                                              s2_vec = s2_temp[[j]],
                                              mu0_vec = rep(mu0, 2),
                                              v02_vec = rep(v02, 2))
        
        
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
    # It's an array 3*10*6
    beta_actual <- vapply(1:p, 
                          function(h) vapply(1:T_final, 
                                             function(s) beta_temp[[h]][labels_temp[[h]][ , s], s],
                                             FUN.VALUE = double(n)),
                          FUN.VALUE = matrix(1, n, T_final))
    
    # I obtain an array 10*101*3 to center the Y with respect to all the splines
    #   but the j-th
    centering <- simplify2array(apply(beta_actual, 1, function(B) B[, -j]%*%t(S[, -j]), simplify = FALSE))
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
                           function(s) sum(1/sigma[which(lab_j[, s] == k)]^2), 
                           FUN.VALUE = double(1))
      # Sum_{y : c_ijt = k} \tilde{y_ixt} / sigma_i^2 --> goes in the mean of 
      #   the posterior in the part that comes from the likelihood
      sum_y_sig <- vapply(1:T_final, 
                          function(s) colSums( t(Y.tilde[s, , which(lab_j[, s] == k)])/
                                                 sigma[which(lab_j[, s] == k)]^2 ), 
                          FUN.VALUE = double(101))
      # Sum_x g_j(x)*"term just computed" --> contribution of the likelihood 
      #   to the posterior mean
      lik2mean <- S[, j] %*% sum_y_sig
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
                                     phi_j = phi_temp[[j]])
    } # Fine ciclo sui cluster per i beta
    
    
    # Aggiungo gli elementi anche negli oggetti non-temp
    labels_res[[j]][ , , d] <- labels_temp[[j]]
    gamma_res[[j]][ , , d] <- gamma_temp[[j]]
    beta_res[[j]][ , , d] <- beta_temp[[j]]
    
    
    
    
    ### ### ### ### ##
    ### UPDATE PHI ###
    phi_res[[j]][d] <- phi_temp[[j]] <- 
      up_phi_j(beta_j = beta_temp[[j]],
               delta_j = delta_temp[[j]],
               SIGinv = SIGMAinv,
               lambda = lambda_temp,
               xi = xi_temp)
    
    
    ### ### ### ### ####
    ### UPDATE DELTA ###
    delta_res[[j]][d] <- delta_temp[[j]] <- up_delta_j(val_now = delta_temp[[j]],
                                                       beta_j = beta_temp[[j]],
                                                       phi_j = phi_temp[[j]],
                                                       SIGchol = cholSIG,
                                                       A_delta = A_delta,
                                                       eps = eps_delta)
    
    
    ### ### ### ### ####
    ### UPDATE ALPHA ###
    alpha_res[[j]][d] <- alpha_temp[[j]] <-
      up_alpha_j(sum(gamma_temp[[j]]),
                 n * T_final,
                 a_alpha,
                 b_alpha)
    
    
    ### ### ### ### #
    ### UPDATE S2 ###
    s2_temp[[j]][1] <- s2_res[[j]][d, 1] <-  up_s2.RWM(s2_vec = s2_temp[[j]],
                                                       idx_cov = 1,
                                                       eps = eps_s2,
                                                       A_s2 = A_s2,
                                                       labels = labels_temp[[j]],
                                                       M = M,
                                                       X = X,
                                                       mu0_vec = rep(mu0, 2),
                                                       v02_vec = rep(v02, 2))
    s2_temp[[j]][2] <- s2_res[[j]][d, 2] <- up_s2.RWM(s2_vec = s2_temp[[j]],
                                                      idx_cov = 2,
                                                      eps = eps_s2,
                                                      A_s2 = A_s2,
                                                      labels = labels_temp[[j]],
                                                      M = M,
                                                      X = X,
                                                      mu0_vec = rep(mu0, 2),
                                                      v02_vec = rep(v02, 2))
  } # END OF FOR LOOP OVER COEFFICIENTS "j"
  
  ### ### ### ### ### #
  ### UPDATE LAMBDA ###
  lambda_res[d] <- lambda_temp <- 
    up_lambda(phi = unlist(phi_temp),
              xi = xi_temp,
              m0 = m0,
              s02 = s02)
  
  
  ### ### ### ### #
  ### UPDATE XI ###
  xi_res[d] <- xi_temp <- 
    up_sd.RWM(val_now = xi_temp,
              eps = eps_xi,
              data = unlist(phi_temp),
              mean = lambda_temp,
              hyppar = A_xi)
  
  
  # Print the time passed since the beginning of the simulation and since the
  #   last tot iteration.
  if ((d %% 100) == 1) {
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

# save.image("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/res/dynamic_on_beta/PPMx/sim_study/v02_1/res_v02_1.RData")
