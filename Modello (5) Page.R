rm(list=ls())

library(tidyverse)
library(splines2)
library(tidyverse)
library(invgamma)

load("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/BSP_Pavone/output/mortality.Rdata")
load("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/fit_indep.RData")
source("funzioni.R")

data_list_man <- list(ita_man = log(Y_ita_man / N_ita_man),
                      swe_man = log(Y_swe_man / N_swe_man),
                      uk_man = log(Y_uk_man / N_uk_man),
                      us_man = log(Y_us_man / N_us_man))

data_list_woman <- list(ita_woman = log(Y_ita_woman / N_ita_woman),
                        swe_woman = log(Y_swe_woman / N_swe_woman),
                        uk_woman = log(Y_uk_woman / N_uk_woman),
                        us_woman = log(Y_us_woman / N_us_woman))

set.seed(4238)





### ### ### ### ### ###
#### Basi B-Spline ####
### ### ### ### ### ###
# Maximum age
Z <- 87
ages <- 0:100
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
S <- rbind(c(1, rep(0, k)),
           cbind(rep(0, length(ages_no0)), S.temp)
           )





### ### ### ### ### ### ### ### ### ### ###
#### Definizione parametri del modello ####
### ### ### ### ### ### ### ### ### ### ###

# Numero istanti temporali
T_final <- min(sapply(c(data_list_man, data_list_woman), 
                      NROW)) # Italy ha un anno in meno degli altri,
# quindi mi fermo al suo ultimo anno
T_final
# Numero stati considerati
n <- length(data_list_man); n
# Dimensione stato latente (-> numero coeff)
p <- ncol(S)
# Varianza delle osservazioni fissata
sigma <- rep(0.1, n)
# Iperparametri della prior dell'ultimo layer
m0 <- -4; s02 <- 1^2
# Uniforms' hyperparameters
#   Following Page's idea in paragraph 2.5 I choose these hyperparams
A_tau <- 1
A_delta <- 1
A_xi <- 1
# Hyperparams for the Beta prior on alpha
a_alpha <- 2; b_alpha <- 3
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
eps_xi <- A_xi





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

theta_res <- tau_res <- replicate(p,
                                  array(NA,
                                        dim = c(T_final, # numero istanti temporali
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


# Temp objects
gamma_temp <- labels_temp <-
  beta_temp <- replicate(p,
                         array(NA,
                                dim = c(n, # numero di osservazioni
                                        ncol = T_final # numero istanti temporali
                                        )
                         ),
                         simplify = FALSE)

theta_temp <- tau_temp <- replicate(p,
                                  array(NA,
                                        dim = c(T_final # numero istanti temporali
                                                )
                                        ),
                                  simplify = FALSE)

phi_temp <- delta_temp <- replicate(p, list(NA))
alpha_temp <- replicate(p, list(NA))
lambda_temp <- xi_temp <- NA



### ### ### ### ### ###
### Inizializzazioni ##
### ### ### ### ### ###
xi_res[1] <- xi_temp <- 
  runif(1, 0, A_xi)
lambda_res[1] <- lambda_temp <- 
  rnorm(1, mean = m0, sd = sqrt(s02))

for (j in 1:p){
  delta_res[[j]][1] <- delta_temp[[j]] <- 
    runif(1, 0.1, A_delta)
  phi_res[[j]][1] <- phi_temp[[j]] <- 
    rnorm(1, lambda_temp, xi_temp)
  
  alpha_res[[j]][1] <- alpha_temp[[j]] <- rbeta(1, a_alpha, b_alpha)
  
  theta_res[[j]][ , 1] <- theta_temp[[j]] <- 
    rnorm(T_final, phi_temp[[j]], delta_temp[[j]])
  tau_res[[j]][ , 1] <- tau_temp[[j]] <- 
    runif(T_final, 0.1, A_tau)
  
  
  # Per ora inizializzo tutti i gamma = 0, così da lasciare piena libertà di 
  #   movimento alla prima iterazione. Poi posso pensare di inizializzare anche
  #   loro dalla prior.
  gamma_res[[j]][ , , 1] <- gamma_temp[[j]][ , ] <- 0
  
  for (t in 1:T_final){
    lab <- rho2lab(rCRP(4, 1))
    labels_temp[[j]][, t] <- labels_res[[j]][, t, 1] <- lab
    beta_res[[j]][1:max(lab), t, 1] <- beta_temp[[j]][1:max(lab), t] <- 
      rnorm(max(lab),
            mean = theta_temp[[j]][t], sd = tau_temp[[j]][t])
  }
}

rm(j); rm(t)




### ### ### ### ### ###
### GIBBS SAMPLING ####
### ### ### ### ### ###
# Ipotizzo di volerlo fare per gli uomini
Y <- lapply(data_list_man, function(y) y[1:87, ])
names(Y) <- substr(names(data_list_man), 1, 2)

inizio <- Sys.time()

for (d in 2:n_iter){ # Ciclo sulle iterazioni
  
  if ((d %% 100) == 0) {cat(d, difftime(Sys.time(), inizio), 
                            as.character(Sys.time()), "\n")}
  
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
                                              lab_tm1 = labels_temp[[j]][, t-1])
        }
      } # Fine ciclo sulle osservazioni per gamma
      
      
      
      ### ### ### ### ### ### ### #
      ### ### UPDATE LABELS ### ###
      for (i in 1:n){
        beta_actual <- vapply(1:p, 
                              function(x) beta_temp[[x]][labels_temp[[x]][i, t], t],
                              FUN.VALUE = double(1))
        
        
        
        
        # Campiono possibile nuovo valore per beta
        newclustervalue <- rnorm(1, theta_temp[[j]][t], tau_temp[[j]][t])
        
        
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
                                              newclustervalue = newclustervalue,
                                              spline_basis = S)
        
        
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
      n_cluster <- max(labels_temp[[j]][ , t])
      
      for (k in 1:n_cluster){ # Inizio ciclo sui cluster per i beta
        beta_temp[[j]][k, t] <- up_beta(j = j, 
                                        k = k, 
                                        t = t,
                                        Y_t = t(vapply(Y, 
                                                       function(y) y[t, ], 
                                                       FUN.VALUE = double(101))),
                                        beta_t =  vapply(beta_temp, 
                                                         function(b) b[, t], 
                                                         FUN.VALUE = double(4)),
                                        theta_jt = theta_temp[[j]][t], 
                                        tau_jt = tau_temp[[j]][t],
                                        sigma_vec = sigma,
                                        labels_t = vapply(labels_temp, 
                                                          function(lab) lab[, t], 
                                                          FUN.VALUE = double(4)),
                                        spline_basis = S)
      } # Fine ciclo sui cluster per i beta
      
      
      
      ### ### ### ### ### 
      ### UPDATE THETA ###
      theta_temp[[j]][t] <- up_theta_jt(beta_jt = beta_temp[[j]][, t],
                                        tau_jt = tau_temp[[j]][t],
                                        phi_j = phi_temp[[j]],
                                        delta_j = delta_temp[[j]])
      
      
      
      ### ### ### ### ###
      ### UPDATE TAU ###
      tau_temp[[j]][t] <- up_sd.RWM(val_now = tau_temp[[j]][t],
                                     eps = eps_tau,
                                     data = beta_temp[[j]][, t],
                                     mean = theta_temp[[j]][t],
                                     hyppar = A_tau)
      
      
    } # END OF FOR LOOP OVER TIME ISTANTS "t"
    
    # Aggiungo gli elementi anche negli oggetti non-temp
    labels_res[[j]][ , , d] <- as.integer(labels_temp[[j]])
    gamma_res[[j]][ , , d] <- gamma_temp[[j]]
    beta_res[[j]][ , , d] <- beta_temp[[j]]
    theta_res[[j]][ , d] <- theta_temp[[j]]
    tau_res[[j]][ , d] <- tau_temp[[j]]
    
    
    ### ### ### ### ##
    ### UPDATE PHI ###
    phi_res[[j]][d] <- phi_temp[[j]] <- 
      up_phi_j(theta_j = theta_temp[[j]],
               delta_j = delta_temp[[j]],
               lambda = lambda_temp,
               xi = xi_temp)
    
    
    ### ### ### ### ####
    ### UPDATE DELTA ###
    delta_res[[j]][d] <- delta_temp[[j]] <- 
      up_sd.RWM(val_now = delta_temp[[j]],
                eps = eps_delta,
                data = theta_temp[[j]],
                mean = phi_temp[[j]],
                hyppar = A_delta)
    
    ### ### ### ### ####
    ### UPDATE ALPHA ###
    alpha_res[[j]][d] <- alpha_temp[[j]] <-
      up_alpha_j(sum(gamma_temp[[j]]),
                 n * T_final,
                 a_alpha,
                 b_alpha)
  } # END OF FOR LOOP OVER COEFFICIENTS "j"
  
  ### ### ### ### ### #
  ### UPDATE LAMBDA ###
  lambda_res[d] <- lambda_temp <- 
    up_lambda(phi = unlist(phi_temp),
              xi = xi_temp,
              m0 = m0,
              s02 = s02)
  
  
  ### ### ### ### ####
  ### UPDATE XI ###
  xi_res[d] <- xi_temp <- 
    up_sd.RWM(val_now = xi_temp,
              eps = eps_xi,
              data = unlist(phi_temp),
              mean = lambda_temp,
              hyppar = A_xi)
  
  
} # END OF FOR LOOP OVER ITERATIONS "d"

fine <- Sys.time()

exec_time <- difftime(fine, inizio)

save.image("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/res/uniform_prior_on_sd_2.RData")
