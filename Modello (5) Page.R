rm(list = ls())
library(doParallel)
library(tidyverse)
library(splines2)
library(tidyverse)
library(invgamma)



# Valori da testare
sigma2_i_val <- c(0.1^2, 5^2)
A_tau_val <- c(1, 50)
A_delta_val <- c(1, 50)
m0_val <- c(0, -4, -10)
s02_val <- c(0.1^2, 10^2)
a_b_val <- rbind(c(1, 1), c(2, 3))
# sigma2_i_val <- 0.1^2
# A_tau_val <- 1
# A_delta_val <- 1
# m0_val <- -4
# s02_val <- 0.1^2
# a_b_val <- rbind(c(2, 3))
grid <- expand_grid(sigma2_i_val, A_tau_val, A_delta_val, 
                    m0_val, s02_val, a_b_val)
grid <- as.matrix(grid)

cl <- makeCluster(8, outfile = "res/sim_only_infant/unif_prior_on_sd/Progress.txt")
# cl <- makeCluster(1, outfile = "Progress.txt")
registerDoParallel(cl)


ng <- nrow(grid)
out <- foreach(row=1:ng) %dopar% {
  
  load("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/BSP_Pavone/output/mortality.Rdata")
  source("funzioni.R")
  
  Y <- rbind(ita_man = log(Y_ita_man / N_ita_man)[1:87, 1],
             swe_man = log(Y_swe_man / N_swe_man)[1:87, 1],
             uk_man = log(Y_uk_man / N_uk_man)[1:87, 1],
             us_man = log(Y_us_man / N_us_man)[1:87, 1])
  rownames(Y) <- substr(rownames(Y), 1, 2)
  
  
  set.seed(1)
  
  
  ### ### ### ### ### ### ### ### ### ### ###
  #### Definizione parametri del modello ####
  ### ### ### ### ### ### ### ### ### ### ###
  
  # Numero istanti temporali
  T_final <- ncol(Y); T_final
  # Numero stati considerati
  n <- nrow(Y); n
  # Iperparametri della prior dell'ultimo layer
  # m0 <- -10; s02 <- 0.1^2
  m0 <- grid[row, "m0_val"]; s02 <- grid[row, "s02_val"]
  # Uniforms' hyperparameters
  #   Following Page's idea in paragraph 2.5 I choose these hyperparams
  A_tau <- grid[row, "A_tau_val"] # 5
  A_delta <- grid[row, "A_delta_val"] # 10
  # A_sigma <- 10
  sigma2_i <- rep(grid[row, "sigma2_i_val"], 4)
  # Hyperparams for the Beta prior on alpha
  a_alpha <- grid[row, "a_b_val.1"] # 2 
  b_alpha <- grid[row, "a_b_val.2"] # 3
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
  eps_tau <- A_tau; eps_tau
  eps_delta <- A_delta; eps_delta
  # eps_sigma <- A_sigma
  
  
  
  
  
  ### ### ### ### ### ### ### ### ### ### ###
  ### Creo gli oggetti in cui salverò gli ###
  ### output dei vari step del Gibbs      ###
  ### ### ### ### ### ### ### ### ### ### ###
  
  # Ho provato ad inizializzare questi oggetti con 99999 sperando di velocizzare,
  #   ma non cambia niente. Perciò, preferisco tenere NA, almeno se qualcosa va
  #   storto me ne accorgo.
  gamma_res <- labels_res <- array(NA,
                                   dim = c(n, # numero di osservazioni
                                           T_final, # numero istanti temporali
                                           n_iter # numero iterazioni
                                   ),
                                   dimnames = list(rownames(Y), colnames(Y), 
                                                   paste("iter", 1:n_iter, sep = "_"))
  )
  
  
  beta_res <- array(NA,
                    dim = c(n, # numero di osservazioni
                            T_final, # numero istanti temporali
                            n_iter # numero iterazioni
                    ),
                    dimnames = list(paste("cl", 1:4, sep = ""), colnames(Y), 
                                    paste("iter", 1:n_iter, sep = "_"))
  )
  
  
  theta_res <- tau_res <- array(NA,
                                dim = c(T_final, # numero istanti temporali
                                        n_iter # numero iterazioni
                                ),
                                dimnames = list(colnames(Y), 
                                                paste("iter", 1:n_iter, sep = "_"))
  )
  
  
  phi_res <- delta_res <- array(NA,
                                dim = c(n_iter # numero iterazioni
                                ),
                                dimnames = list(paste("iter", 1:n_iter, sep = "_"))
  )
  
  
  alpha_res <- array(NA,
                     dim = c(n_iter # numero iterazioni
                     ),
                     dimnames = list(paste("iter", 1:n_iter, sep = "_"))
  )
  
  
  sigma_res <- array(NA,
                     dim = c(n, # numero osservazioni
                             n_iter # numero iterazioni
                     ),
                     dimnames = list(rownames(Y), paste("iter", 1:n_iter, sep = "_"))
  )
  
  
  # Temp objects
  gamma_temp <- labels_temp <-
    array(NA,
          dim = c(n, # numero di osservazioni
                  ncol = T_final # numero istanti temporali
          ),
          dimnames = list(rownames(Y), colnames(Y))
    )
  
  beta_temp <- 
    array(NA,
          dim = c(n, # numero di osservazioni
                  ncol = T_final # numero istanti temporali
          ),
          dimnames = list(paste("cl", 1:4, sep = ""), colnames(Y))
    )
  
  theta_temp <- tau_temp <- array(NA,
                                  dim = c(T_final # numero istanti temporali
                                  ),
                                  dimnames = list(colnames(Y))
  )
  
  
  phi_temp <- delta_temp <- NA
  alpha_temp <- NA
  # alpha_temp <- 0
  # sigma_temp <- rep(NA, n)
  sigma_temp <- sigma2_i
  names(sigma_temp) <- rownames(Y)
  
  
  
  ### ### ### ### ### ###
  ### Inizializzazioni ##
  ### ### ### ### ### ###
  # sigma_res[, 1] <- sigma_temp <- runif(n,
  # 0, A_sigma)
  sigma_res[, 1] <- sigma_temp <- sigma2_i
  
  alpha_res[1] <- alpha_temp <- rbeta(1, 
                                      a_alpha, b_alpha)
  
  delta_res[1] <- delta_temp <-runif(1, 
                                     0.1, A_delta)
  phi_res[1] <- phi_temp <- rnorm(1, 
                                  mean = m0, sd = sqrt(s02))
  
  
  theta_res[ , 1] <- theta_temp <- rnorm(T_final, 
                                         mean = phi_temp, sd = delta_temp)
  tau_res[ , 1] <- tau_temp <- runif(T_final, 
                                     0.1, A_tau)
  
  
  # Per ora inizializzo tutti i gamma = 0, così da lasciare piena libertà di 
  #   movimento alla prima iterazione. Poi posso pensare di inizializzare anche
  #   loro dalla prior.
  gamma_res[ , , 1] <- gamma_temp[ , ] <- 0
  
  for (t in 1:T_final){
    # lab <- rho2lab(rCRP(4, M))
    # Per ora inizializzo le unità in 4 cluster per forzarli a stare separati
    lab <- c(1, 2, 3, 4)
    labels_temp[, t] <- labels_res[, t, 1] <- lab
    
    beta_res[1:max(lab), t, 1] <- beta_temp[1:max(lab), t] <- 
      rnorm(max(lab),
            mean = theta_temp[t], sd = tau_temp[t])
  }
  rm(t)
  
  
  # debug(up_beta)
  # debug(up_theta_t)
  
  ### ### ### ### ### ###
  ### GIBBS SAMPLING ####
  ### ### ### ### ### ###
  
  inizio <- Sys.time()
  
  for (d in 2:n_iter){ # Ciclo sulle iterazioni
    
    if ((d %% 500) == 0) {cat("Job", row, "iter", d, Sys.time() - inizio, "\n")}
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
          gamma_temp[i, 1] <- 0
        } else if (t > 1){
          # Assegno la nuova gamma nel tmp
          gamma_temp[i, t] <- up_gamma_i(i = i, 
                                         gamma = gamma_temp[, t], 
                                         alpha_t = alpha_temp,
                                         lab_t = labels_temp[, t],
                                         lab_tm1 = labels_temp[, t-1],
                                         M = M)
        }
      } # Fine ciclo sulle osservazioni per gamma
      
      
      
      ### ### ### ### ### ### ### #
      ### ### UPDATE LABELS ### ###
      for (i in 1:n){
        
        # Campiono possibile nuovo valore per beta
        newclustervalue <- rnorm(1, 
                                 mean = theta_temp[t], 
                                 sd = tau_temp[t])
        
        
        # Non assegno qua anche a labels_res perché dovrò prima sistemarlo
        labels_temp[i, t] <-  up_label_i(i = i, 
                                         Y_it = Y[i, t],
                                         beta_cluster = beta_temp[, t],
                                         sigma2_i = sigma_temp[i],
                                         lab_t = labels_temp[,t],
                                         lab_tp1 = labels_temp[,t+1], 
                                         gamma_tp1 = if (t == T_final) {'last time'} 
                                         else {gamma_temp[, t+1]},
                                         newclustervalue = newclustervalue,
                                         M = M)
        
        
        # Sistemo le labels in modo tale che non ci siano buchi (tipo non 
        # voglio c(1, 2, 2, 4) ma c(1, 2, 2, 4)) e che siano ordinate
        # (tipo non voglio c(2, 1, 2, 3) ma c(1, 2, 1, 3)).
        labels_new <- c(1L)
        counter <- 2L
        for (ii in 2:n){
          if (!(labels_temp[ii, t] %in% labels_temp[1:(ii-1), t])){
            labels_new[ii] <- counter
            counter <- counter+1
          } else{
            labels_new[ii] <- labels_new[which(labels_temp[, t] == 
                                                 labels_temp[ii, t])[1]]
          }
        }
        
        # Voglio riordinare anche beta_temp in accordo con le labels:
        #   facendo unique(...) ottengo i valori unici delle labels in ordine
        #   di apparizione nel vettore (e non in ordine cresc./decresc.).
        #   Questo vettore dà la permutazione necessaria.
        
        if (all(labels_temp[i, t] > labels_temp[-i, t])){
          # Se la label dell'unità i è maggiore di tutte le altre, vuol dire
          #   che è andato in un nuovo cluster e non in uno di quelli esistenti,
          #   quindi aggiungo il beta del nuovo cluster al vettore dei beta
          #   dei cluster.
          if (labels_temp[i, t] > n){
            # If the new cluster is the (n+1)st, then I cannot modify directly
            #   beta_temp because I would exceed its size.
            beta_extended <- c(beta_temp[, t], newclustervalue)
            sort_beta_temp <- unique(labels_temp[, t])
            beta_new <- beta_extended[sort_beta_temp]
          } else { # If new cluster is <= n
            beta_temp[labels_temp[i, t], t] <- newclustervalue
            sort_beta_temp <- unique(labels_temp[, t])
            beta_new <- beta_temp[sort_beta_temp, t]
          }
        } else { # If not new cluster
          sort_beta_temp <- unique(labels_temp[, t])
          beta_new <- beta_temp[sort_beta_temp, t]
        }
        
        beta_temp[, t] <- NA
        beta_temp[1:length(sort_beta_temp), t] <- beta_new
        
        # Finally, I assign the new fixed label to the storing object
        labels_temp[, t] <- labels_new
        
      } # Fine ciclo sulle osservazioni per labels
      
      
      
      ### ### ### ### ### ### ###
      ### ### UPDATE BETA ### ###
      n_cluster <- max(labels_temp[ , t])
      
      for (k in 1:n_cluster){ # Inizio ciclo sui cluster per i beta
        beta_temp[k, t] <- up_beta(k = k, 
                                   t = t,
                                   Y_t = Y[, t],
                                   theta_t = theta_temp[t], 
                                   tau_t = tau_temp[t],
                                   sigma2_vec = sigma_temp,
                                   labels_t = labels_temp[, t])
      } # Fine ciclo sui cluster per i beta
      
      
      
      ### ### ### ### ### 
      ### UPDATE THETA ###
      theta_temp[t] <- up_theta_t(beta_t = beta_temp[, t],
                                  tau_t = tau_temp[t],
                                  phi = phi_temp,
                                  delta = delta_temp)
      
      
      
      ### ### ### ### ###
      ### UPDATE TAU ###
      tau_temp[t] <- up_sd.RWM(val_now = tau_temp[t],
                                eps = eps_tau,
                                data = beta_temp[, t][!is.na(beta_temp[, t])],
                                mean = theta_temp[t],
                                hyppar = A_tau)
      
      
    } # END OF FOR LOOP OVER TIME ISTANTS "t"
    
    # Aggiungo gli elementi anche negli oggetti non-temp
    labels_res[ , , d] <- as.integer(labels_temp)
    gamma_res[ , , d] <- gamma_temp
    beta_res[ , , d] <- beta_temp
    theta_res[ , d] <- theta_temp
    tau_res[ , d] <- tau_temp
    
    
    ### ### ### ### ##
    ### UPDATE PHI ###
    phi_res[d] <- phi_temp <-
      up_phi(theta = theta_temp,
             delta = delta_temp,
             m0 = m0,
             s02 = s02)
    
    
    ### ### ### ### ####
    ### UPDATE DELTA ###
    delta_res[d] <- delta_temp <-
      up_sd.RWM(val_now = delta_temp,
                 eps = eps_delta,
                 data = theta_temp,
                 mean = phi_temp,
                 hyppar = A_delta)
    
    
    
    ### ### ### ### ####
    ### UPDATE ALPHA ###
    alpha_res[d] <- alpha_temp <-
      up_alpha(sum(gamma_temp),
               n * T_final,
               a_alpha,
               b_alpha)
    
    
    
    ### ### ### ### ### ##
    ### UPDATE SIGMA_i ###
    sigma_res[, d] <- sigma_temp
    # for (i in 1:n){
    #   # T times 1 vector with the coefficients of the splines for the i-th obs
    #   beta <- beta_temp[cbind(labels_temp[i, ], 1:87)]
    #   sigma_res[i, d] <- sigma_temp[i] <-
    #     up_var.RWM(val_now = sigma_temp[i],
    #                eps = eps_sigma,
    #                data = Y[i, ],
    #                mean = beta,
    #                hyppar = A_sigma)
    # }
    
  } # END OF FOR LOOP OVER ITERATIONS "d"
  
  print("finito")
  
  fine <- Sys.time()
  
  exec_time <- fine - inizio; exec_time
  
  # Organize items
  params <- list(T_final = T_final,
                 n = n,
                 n_iter = n_iter,
                 eps_tau = eps_tau,
                 eps_delta = eps_delta,
                 sigma2_i = sigma2_i,
                 hyperpar = list(A_tau = A_tau,
                                 A_delta = A_delta,
                                 m0 = m0, 
                                 s02 = s02,
                                 a_alpha = a_alpha, 
                                 b_alpha = b_alpha,
                                 M = M))
  
  res <- list(labels = labels_res,
              gamma = gamma_res,
              beta = beta_res,
              theta = theta_res, 
              tau = tau_res,
              phi = phi_res,
              delta = delta_res,
              alpha = alpha_res,
              sigma = sigma_res,
              exec_time = exec_time)
  
  functions <- list(GGup_iid = GaussGaussUpdate_iid,
                    up_label = up_label_i,
                    up_gamma = up_gamma_i,
                    up_beta = up_beta,
                    up_theta = up_theta_t,
                    up_phi = up_phi,
                    up_alpha = up_alpha,
                    up_sd = up_sd.RWM)
  
  data <- list(Y = Y,
               gender = "man")
  
  print("finito")
  
  fine <- Sys.time()
  
  exec_time <- fine - inizio; exec_time
  
  # Organize items
  params <- list(T_final = T_final,
                 n = n,
                 n_iter = n_iter,
                 eps_tau = eps_tau,
                 eps_delta = eps_delta,
                 sigma2_i = sigma2_i,
                 hyperpar = list(A_tau = A_tau,
                                 A_delta = A_delta,
                                 m0 = m0, 
                                 s02 = s02,
                                 a_alpha = a_alpha, 
                                 b_alpha = b_alpha,
                                 M = M))
  
  res <- list(labels = labels_res,
              gamma = gamma_res,
              beta = beta_res,
              theta = theta_res, 
              tau = tau_res,
              phi = phi_res,
              delta = delta_res,
              alpha = alpha_res,
              sigma = sigma_res,
              exec_time = exec_time)
  
  functions <- list(GGup_iid = GaussGaussUpdate_iid,
                    up_label = up_label_i,
                    up_gamma = up_gamma_i,
                    up_beta = up_beta,
                    up_theta = up_theta_t,
                    up_phi = up_phi,
                    up_alpha = up_alpha,
                    up_sd = up_sd.RWM)
  
  data <- list(Y = Y,
               gender = "man")
  
  cat("Job", row, "fatto", "\n")
  # rm(list = setdiff(ls(), c("data", "params", "res", "functions")))
  list(params, res, functions, data)
}


stopCluster(cl)
save.image("res/sim_only_infant/unif_prior_on_sd/comparisons_hyppar.RData")
