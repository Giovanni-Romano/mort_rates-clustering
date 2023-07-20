library(here)
library(tidyverse)
library(KFAS)
library(mvtnorm)
library(mvnfast)
library(mggd)


#### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### Funzione per invertire matrici pd, come varcov ####
#### ### ### ### ### ### ### ### ### ### ### ### ### ###
mysolve <- function(x){
  chol2inv(chol(x))
}

#### ### ### ### ### ### ### ### ### ### ###
#### Da labels a partizione e viceversa ####
#### ### ### ### ### ### ### ### ### ### ###
lab2rho <- function(labels){
  nTbl <- max(labels)
  listTableOccupants <- lapply(1:nTbl, function(t) which(labels == t)) 
  return(listTableOccupants)
}

rho2lab <- function(rho){
  n <- max(unlist(rho))
  # lab <- sapply(1:n, function(x) which(sapply(rho, function(y) x %in% y)))
  lab <- vapply(1:n, 
                function(x) which(vapply(rho, is.element, el = x, FUN.VALUE = FALSE)), FUN.VALUE = integer(1))
  return(lab)
}


### ### ### ### ### ### ###
#### DISTRIBUZIONE CRP ####
### ### ### ### ### ### ###
dCRP <- function(size, M){
  # size è il vettore con le cluster size e M è il parametro di concentrazione
  size <- size[size>0]
  M^NROW(size) * prod(factorial(size - 1)) / pochhammer(M, sum(size))
}

rCRP <- function(n, M){
  labels <- rep(0, n)
  
  labels[1] <- 1
  for (dnr in 2:n) {
    
    # compute occupation probabilities for current diner
    vOcc <- table(labels[1:(dnr-1)])
    vProb <- c(vOcc, M) / (dnr - 1 + M)
    
    # add table label to diner
    nTbl <- as.numeric(names(vOcc)[length(vOcc)])  # avoid overhead of finding max of possibly large vector
    labels[dnr] <- sample.int(nTbl+1, size=1, prob=vProb)
  }
  
  nTbl <- max(c(nTbl, labels[n]))
  listTableOccupants <- lapply(1:nTbl, function(t) which(labels == t))
  
  return(listTableOccupants)
}


### ### ### ### ### ### ### ### ### ##
#### CONTROLLO PER COMPATIBILITA' ####
### ### ### ### ### ### ### ### ### ##
are.partitions.equal <- function(part1, part2){
  # part1_sort <- part1 %>%
  #   sapply(., '[[', 1) %>% 
  #   order() %>% 
  #   part1[.] %>% 
  #   lapply(., sort)
  # 
  # part2_sort <- part2 %>%
  #   sapply(., '[[', 1) %>% 
  #   order() %>% 
  #   part2[.] %>% 
  #   lapply(., sort)
  # 
  # # Se la setdiff tra gli insiemi delle due liste ha lunghezza 0, vuol dire
  # # che tutti gli insiemi sono uguali
  # check <- (length(setdiff(part1_sort, part2_sort)) == 0)
  
  
  
  ### ### ### ### ### ### ### ### ### ### ### ### #
  ### POSSIBILE ALTERNATIVA --> CONTROLLARE!!!  ###
  ### SE FA LA COSA GIUSTA, E' MOLTO PIU' VELOCE ##
  ### ### ### ### ### ### ### ### ### ### ### ### #
  check <- all(
    vapply(part1, 
           function(x) any(vapply(part2, function(y) setequal(x, y), FALSE)),
           FALSE)
    )
  
  return(check)
}





### ### ### ### ### ##
#### UPDATE GAMMA ####
### ### ### ### ### ##
up_gamma_i <- function(i, gamma, alpha_t, rho_t, rho_tm1){
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### - gamma è il vettore gamma_{t} più recente, quindi alcune componenti
  ###     saranno già aggiornate all'iterazione (d) e altre saranno ancora alla (d-1)
  ### - rho_tm1 è rho_{t-1}
  ### - nota che rho_{t-1} sarà già aggiornato all'iter. (d) mentre rho_t è ancora
  ###     alla vecchia iter (d-1)
  
  # R_t^{(+i)}: insieme di unità che rimangono fisse da t-1 a t
  R_t <- which(gamma == 1)
  R_tmi <- R_t[R_t != i]
  R_tpi <- c(R_tmi, i)
  
  # uso il ".R" per indicare le partizioni ridotte
  rho_tm1.R <- list()
  rho_t.R <- list()
  
  for (j in 1:length(rho_tm1)){
    if (sum(rho_tm1[[j]] %in% R_tpi) > 0){
      rho_tm1.R[[length(rho_tm1.R) + 1]] <- rho_tm1[[j]][rho_tm1[[j]] %in% R_tpi]
    }
  }
  
  for (j in 1:length(rho_t)){
    if (sum(rho_t[[j]] %in% R_tpi) > 0){
      rho_t.R[[length(rho_t.R) + 1]] <- rho_t[[j]][rho_t[[j]] %in% R_tpi]
    }
  }
  
  
  ############################################################################
  ### Voglio ordinare le "partizioni ridotte" al tempo t-1 e t per poterle ###
  ### confrontare facilmente. Ordino ciascun vettore all'interno delle     ###
  ### liste e poi ordino tra di loro i vettori in base al primo elemento   ###
  ############################################################################
  
  check <- are.partitions.equal(rho_tm1.R, rho_t.R)
  
  
  # Vedi formula S.9 del materiale supplementare di Page
  # Devo calcolare la predictive per l'i-esima osservazione data la 
  #   partizione ridotta a R_{t-1}^{(-i)}
  
  # Trova insieme della partizione con dentro i
  j <- which(vapply(rho_t.R, is.element, el = i, FUN.VALUE = FALSE))
  n <- sum(vapply(rho_t.R, length, FUN.VALUE = integer(1))) - 1
  
  if (n == 0) { # se la partizione a cui ci condizioniamo è vuota
    ratio <- 1
  } else {
    if (length(rho_t.R[[j]]) == 1) {
      # se l'unità "i" è in un nuovo cluster
      ratio <- M / (n + M)
    } else { # se l'unità "i" è in un cluster già esistente
      ratio <- (length(rho_t.R[[j]]) - 1) / (n + M)
    }
  }
  
  prob <- alpha_t / (alpha_t + (1 - alpha_t) * ratio) * check
  
  gamma_it <- as.integer(runif(1) < prob)
  
  return(gamma_it)
}





### ### ### ### ### ### ### ###
#### UPDATE CLUSTER LABELS ####
### ### ### ### ### ### ### ###
up_label_i <- function(i, j,
                       Y_it,
                       beta_i,
                       beta_cluster,
                       sigma2_i,
                       rho_t, rho_tp1, gamma_tp1,
                       newclustervalue,
                       spline_basis) {
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### - i: indice dell'unità che stiamo aggiornando
  ### - j: indice del coefficiente che stiamo aggiornando
  ### - Y_it: vettore con Y_{ixt} per ogni x
  ### - beta_i: vettore con i beta dell'osservazione i-esima al tempo t
  ### - beta_cluster: vettore con i beta dei cluster del j-esimo coeff. al tempo t
  ### - sigma2_i: varianza dell'i-esima osservazione
  ### - rho_t: partizione al tempo t PER IL COEFF j
  ### - rho_tp1: partizione al tempo t+1 PER IL COEFF j
  ### - gamma_tp1: gamma_{t+1}, per trovare la partizione ridotta;
  ###               all'istante finale passare 'last time'
  ### - newclustevalue: valore per il nuovo cluster
  ### - spline_basis: matrice con i valori delle basi spline
  
  # Inizializzazione vettore probabilità della multinomiale
  logp <- c()
  
  # Cluster in cui è l'i-esima osservazione
  k <- which(vapply(rho_t, is.element, el = i, FUN.VALUE = FALSE))
  
  # Partizione senza i-esima osservazione
  rho_tmi <- rho_t
  rho_tmi[[k]] <- rho_tmi[[k]][rho_tmi[[k]] != i]
  
  whichclusters <- which(vapply(rho_tmi, function(x) length(x)>0, FUN.VALUE = FALSE)) # cluster esistenti
  whichclusters <- c(whichclusters, max(whichclusters)+1L) # aggiungo il possibile cluster nuovo
  # Uso "1L" per farlo venire integer, altrimenti il "+1" lo rende un float
  
  for (h in whichclusters){
    # Devo costruire la partizione con l'i-esima osservazione inserita
    #   nell'h-esimo cluster
    rho_t.h <- rho_tmi
    if (h < max(whichclusters)){ # se assegno ad un cluster già esistente
      rho_t.h[[h]] <- c(rho_t.h[[h]], i) 
    } else { # se assegno al nuovo cluster
      rho_t.h[[h]] <- i 
    }
    
    # Nel paper di Page: Pr(c_it = h).
    prob_cit <- dCRP(vapply(rho_t.h, length, FUN.VALUE = integer(1)), M)
    
    # Nel paper di Page: N(Y_it | bla bla).we
    # Per me ci sarà il prodotto di Pr(Y_{ixt} |  beta, ecc) per tutti gli x.
    if (h < max(whichclusters)){ # se il cluster è già esistente, devo usare il beta di quel cluster
      beta <- beta_cluster[h]
    } else { # se il cluster è nuovo, devo simulare la media dalla prior
      beta <- newclustervalue
    }
    # I beta che devo usare sono il vettore beta_i, sostituendo
    #   per il j-esimo coeff. il beta che ho appena calcolato.
    #
    beta_i[j] <- beta
    
    # Nel paper di Page: indicatrice sulle partizioni ridotte
    #   per prima cosa devo trovare la partizione ridotta
    if (identical(gamma_tp1, 'last time')){
      check <- TRUE
    } else {
      R_tp1 <- which(gamma_tp1 == 1)
      rho_tp1.R <- list()
      rho_t.h.R <- list()
      
      for (kk in 1:length(rho_t.h)){
        if (sum(rho_t.h[[kk]] %in% R_tp1) > 0){
          rho_t.h.R[[length(rho_t.h.R) + 1]] <- rho_t.h[[kk]][rho_t.h[[kk]] %in% R_tp1]
        }
      }
      
      for (kkk in 1:length(rho_tp1)){
        # Ho aggiunto il seguente "if" per non mettere vettori 
        #   vuoti in rho_tp1.R
        if (sum(rho_tp1[[kkk]] %in% R_tp1) > 0){
          rho_tp1.R[[length(rho_tp1.R) + 1]] <- rho_tp1[[kkk]][rho_tp1[[kkk]] %in% R_tp1]
        }
      }
      check <- are.partitions.equal(rho_t.h.R, rho_tp1.R)
    }
    
    means <- drop(spline_basis %*% beta_i)
    logpnorm <- sum(dnorm(Y_it, 
                          mean = means,
                          sd = sqrt(sigma2_i),
                          log = TRUE))
    
    logp <- c(logp, 
              log(prob_cit) + logpnorm + log(check))
  }
  
  
  # Almeno nelle prime iterazioni sembra che le prob siano tutte 0, 
  #   quindi devo lavorare in scala logaritmica. Per campionare dalle 
  #   log-prob uso il trick del max con la Gumbel.
  lll <- logp + min(abs(min(logp)), .Machine$double.xmax)
  gumbel <- -log(-log(runif(length(lll))))
  lll <- lll + gumbel
  
  c_it <- whichclusters[which.max(lll)]
  return(c_it)
}






## ### ### ### ### ##
#### UPDATE BETA ####
## ### ### ### ### ##
up_beta <- function(j, k, t,
                    Y_t,
                    beta_t,
                    theta_jt, tau_jt,
                    sigma2_vec,
                    labels_t,
                    spline_basis){
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### - j: indice del coefficiente di cui stiamo aggiornando la media
  ### - k: indice del cluster di cui stiamo aggiornando la media
  ### - t: istante temporale attuale
  ### - Y_t: matrix with Y_{ixt} for all obs. i and all obs. x
  ### - beta_t: matrice con i beta per tutti i cluster per tutti i coeff. 
  ###           al tempo t
  ### - theta_jt, tau_jt: media e varianza della prior dei beta
  ### - sigma2_vec: vettore con le varianze delle osservazioni
  ### - labels_t: matrice delle label al tempo t per tutti i coeff.
  ### - spline_basis: matrice con i valori delle basi spline
  
  # Devo trovare quali osservazioni appartengono al cluster k per il 
  #   e coeff. j e poi uso conjugacy Gauss-Gauss.
  obs <- lab2rho(labels_t[ , j])[[k]]
  
  
  # Lavoro su \tilde{Y}_ixt, cioè le osservazioni Y_{ixt} meno le spline tranne
  #   la j-esima, così che beta_jkt*g_j(x) sia la media di questa
  #   nuova variabile. Inoltre tengo solo le osservazioni che servono, cioè 
  #   quelle che per il coefficiente j appartengono al cluster k
  
  # Devo recuperare i beta giusti per ogni osservazione per ogni coeff. != j
  p <- ncol(spline_basis)
  beta_actual <- vapply(1:p, 
                        function(x) beta_t[labels_t[ , x], x],
                        FUN.VALUE = double(4))[obs, -j] # tolgo la j-esima colonna
                                                        # e le oss. non nel cluster
  # Costruisco effettivamente \tilde{Y}_t
  Y_t.tilde <- Y_t[obs, ] - beta_actual %*% t(spline_basis[ , -j])
  # Se Y_t.tilde ha una riga sola, lo trasformo in vettore; se ha più di 
  #   una riga, non succede niente.
  Y_t.tilde <- drop(Y_t.tilde)
  
  
  # Devo togliere le sigma2_i delle osservazioni che non usiamo
  sigma2_vec <- sigma2_vec[obs]
  
  
  # Varianza a posteriori
  var.post <- ( 1 / tau_jt + sum(1/sigma2_vec)*sum(spline_basis[, j]^2) )^(-1)
  # Media a posteriori
  mean.post <- var.post * 
    ( sum( t(t(Y_t.tilde)*spline_basis[ , j])/sigma2_vec ) + 
        theta_jt / tau_jt )
  # Devo fare il doppio trasposto per sfruttare bene il prodotto element-wise
  
  beta_updated <- rnorm(1, mean = mean.post, sd = sqrt(var.post))
  
  return(beta_updated)
}





### ### ### ### ### ### ### ### ##
#### GAUSS-GAUSS CONJ. UPDATE ####
### ### ### ### ### ### ### ### ##
GaussGaussUpdate_iid <- function (xbar, n, datavar, priormean, priorvar) 
{
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### - xbar: sample mean of data
  ### - datavar: theoretical variance of the distribution of the data
  ### - priormean: mean of the Gaussian prior
  ### - priorvar: variance of the Gaussian prior
  
  postvar <- (n/datavar + 1/priorvar)^(-1)
  postmean <- postvar * (n * xbar/datavar + priormean/priorvar)
  
  res <- c(postmean = postmean, postvar = postvar)
  return(res)
}





### ### ### ### ### ##
#### UPDATE THETA ####
### ### ### ### ### ##
up_theta_jt <- function(beta_jt,
                        tau_jt,
                        phi_j,
                        delta_j){
  ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### - beta_jt: vector of beta_kjt for all clusters k
  ### - the other parameters' names are clear
  
  x <- beta_jt[!is.na(beta_jt)]
  xbar <- mean(x)
  n <- length(x)
  
  postpar <- GaussGaussUpdate_iid(xbar = xbar,
                                  n = n,
                                  datavar = tau_jt,
                                  priormean = phi_j,
                                  priorvar = delta_j)
  
  out <- rnorm(1, postpar[1], sqrt(postpar[2]))
  
  return(out)
}





### ### #### ### ### 
#### UPDATE PHI ####
### ### #### ### ### 
up_phi_j <- function(theta_j,
                     delta_j,
                     lambda,
                     xi){
  ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### - theta_j: vector with theta_jt for all t
  ### - the other parameters' names are clear
  
  x <- theta_j
  xbar <- mean(x)
  n <- length(x)
  
  postpar <- GaussGaussUpdate_iid(xbar = xbar,
                                  n = n,
                                  datavar = delta_j,
                                  priormean = lambda,
                                  priorvar = xi)
  
  out <- rnorm(1, postpar[1], sqrt(postpar[2]))
  
  return(out)
}





### ### #### ### ### ## 
#### UPDATE LAMBDA ####
### ### #### ### ### ##
up_lambda <- function(phi,
                     xi,
                     m0,
                     s02){
  ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### - phi: vector with phi_j for all j
  ### - the other parameters' names are clear
  
  x <- phi
  xbar <- mean(x)
  n <- length(x)
  
  postpar <- GaussGaussUpdate_iid(xbar = xbar,
                                  n = n,
                                  datavar = xi,
                                  priormean = m0,
                                  priorvar = s02)
  
  out <- rnorm(1, postpar[1], sqrt(postpar[2]))
  
  return(out)
}





# ## ### ### ### ### #
# #### UPDATE TAU ####
# ## ### ### ### ### #
# up_tau <- function(tau_now,
#                    epsilon,
#                    beta_t,
#                    theta_t,
#                    A_tau){
#   ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#   ### - tau_now: vector w/ current values of tau_jt for all j
#   ### - eps: RW step size
#   ### - beta_t: n*p matrix of beta_{kjt} for all clusters k and for all coeff. j
#   ### - theta_t: p-vector of theta_{jt} for all coeff.'s j
#   ### - A_tau: hyperparam of prior on tau's
#   
#   # Sample proposed value
#   tau_prop <- runif(ncol(beta_t), tau_now - eps, tau_now + eps)
#   
#   
#   # Vector with prob's of acceptance
#   
#   # Indicator to check if proposed values are in the domain
#   indicator <- log(as.integer(tau_prop > 0 & tau_prop < A_tau))
#   
#   # Likelihood of current values and proposed ones
#   #   I have to transpose beta_t because R fills matrices by columns and beta_t
#   #   is n*p and theta_t is p*1.
#   likel_now <- rowSums(dnorm(t(beta_t), 
#                              mean = theta_t, 
#                              sd = sqrt(tau_now),
#                              log = TRUE),
#                        na.rm = TRUE) 
#   likel_prop <- rowSums(dnorm(t(beta_t), 
#                               mean = theta_t,
#                               sd = sqrt(tau_prop),
#                               log = TRUE),
#                         na.rm = TRUE)
#   logprob <- indicator + likel_now - likel_prop
#   
#   # Create output
#   out <- tau_now
#   # If prob. of accept. is >= than 1 (so log>=0), accept proposed value
#   out[logprob >= 0] <- tau_prop[logprob >= 0]
#   
#   return(out)
# }





### ### ### ### ### ### ##
#### UPDATE VARIANCES ####
### ### ### ### ### ### ##
# Function for the RW Metropolis (RWM) for tau, delta and xi
# It is written to update one param. at a time
up_var.RWM <- function(val_now,
                   eps,
                   data,
                   mean,
                   hyppar){
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### - val_now: vector w/ current value of tau_jt/delta_t/xi
  ### - eps: RW step size
  ### - data: object that as the role of data in the posterior update, so
  ###         beta_kjt(for all k)/theta_jt(for all j)/phi_t(for all t) respectively
  ### - mean: the mean of the distrib. of the "data", so respectively
  ###         theta_jt/phi_t/lambda
  ### - hyppar: hyperparam of prior, so respectively A_tau/A_delta/A_xi
  
  # Sample proposed value
  val_prop <- runif(1, val_now - eps, val_now + eps)
  
  
  # Vector with prob's of acceptance
  
  # Indicator to check if proposed values are in the domain
  indicator <- (val_prop > 0 & val_prop < hyppar)
  
  if (indicator){
    # Likelihood of current value and proposed one
    likel_now <- sum(dnorm(data, 
                           mean = mean, 
                           sd = sqrt(val_now),
                           log = TRUE),
                     na.rm = TRUE) 
    likel_prop <- sum(dnorm(data, 
                            mean = mean,
                            sd = sqrt(val_prop),
                            log = TRUE),
                      na.rm = TRUE)
    logprob <- indicator + likel_prop - likel_now
  } else {
    logprob <- -Inf
  }
  
  # Create output
  acc <- (log(runif(1, 0, 1)) < logprob) 
  out <- ifelse(acc, val_prop, val_now)
  
  return(out)
}
