library(tidyverse)
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
  
  
  # Unit: microseconds
  # expr        min    lq     mean      median    uq      max       neval
  # my.func     3.4   3.601   3.985279  3.701     3.802   61.400    1000
  # nimble      7.0   7.301   8.135255  7.401     7.601   206.801   1000
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
  
  # Unit: microseconds (1e-6 seconds)
  # expr    min     lq      mean median      uq      max neval
  # old 87.600 90.301 109.87789 92.601 100.801 9225.001 10000
  # new 47.702 49.301  59.84257 50.301  54.201 7147.502 10000
  
  check <- all(
    vapply(part1, 
           function(x) any(vapply(part2, function(y) setequal(x, y), FALSE)),
           FALSE)
    )
  
  return(check)
}





### ### ### ### ### ### ### ### ###
#### SQUARED EXPONENTIAL KERNEL ####
### ### ### ### ### ### ### ### ###
sq_exp_ker <- function(d, l){
  ### ### ### ### ### ### ### ### ### ### ###
  ### - compute the squared exponential kernel in d = x1 - x2: exp{-(d^2)/(2l^2)}
  ### - l control the smoothness of the kernel (the larger l and the greater the smoothness)
  ### - the term to control the scale (which is included for example in BDA3 is 
  ###     not considered here because I take it outside the kernel and i put a 
  ###     prior on it)
  exp( - d^2 / (2*l^2) )
}





### ### ### ### ### ##
#### UPDATE GAMMA ####
### ### ### ### ### ##
up_gamma_i <- function(i, gamma, alpha_t, lab_t, lab_tm1){
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### - gamma è il vettore gamma_{t} più recente, quindi alcune componenti
  ###     saranno già aggiornate all'iterazione (d) e altre saranno ancora alla (d-1)
  ### - rho_tm1 è rho_{t-1}
  ### - nota che rho_{t-1} sarà già aggiornato all'iter. (d) mentre rho_t è ancora
  ###     alla vecchia iter (d-1)
  
  # R_t^{(+i)}: insieme di unità che rimangono fisse da t-1 a t unito a i
  R_t <- which(gamma == 1)
  R_tmi <- R_t[R_t != i]
  R_tpi <- c(R_tmi, i)
  
  # uso il ".R" per indicare le partizioni ridotte
  lab_tm1.R <- rep(-1, length(lab_t))
  lab_tm1.R[R_tpi] <- lab_tm1[R_tpi]
  lab_t.R <- rep(-1, length(lab_t))
  lab_t.R[R_tpi] <- lab_t[R_tpi]
  
  check <- all(lab_tm1.R == lab_t.R)
  
  
  # Vedi formula S.9 del materiale supplementare di Page
  # Devo calcolare la predictive per l'i-esima osservazione data la 
  #   partizione ridotta a R_{t}^{(-i)}
  
  
  # Trova insieme della partizione con dentro i
  j <- lab_t.R[i]
  n.R <- sum(lab_t.R > 0) - 1
  n_j <- sum(lab_t.R == j)
  
  if (n.R == 0) { # se la partizione a cui ci condizioniamo è vuota
    ratio <- 1
  } else {
    if (n_j == 1) {
      # se l'unità "i" è in un nuovo cluster
      ratio <- M / (n.R + M)
    } else { # se l'unità "i" è in un cluster già esistente
      ratio <- (n_j - 1) / (n.R + M)
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
                       sigma_i,
                       lab_t, lab_tp1, gamma_tp1,
                       newclustervalue,
                       spline_basis) {
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### - i: indice dell'unità che stiamo aggiornando
  ### - j: indice del coefficiente che stiamo aggiornando
  ### - Y_it: vettore con Y_{ixt} per ogni x
  ### - beta_i: vettore con i beta dell'osservazione i-esima al tempo t
  ### - beta_cluster: vettore con i beta dei cluster del j-esimo coeff. al tempo t
  ### - sigma_i: sd dell'i-esima osservazione
  ### - rho_t: partizione al tempo t PER IL COEFF j
  ### - rho_tp1: partizione al tempo t+1 PER IL COEFF j
  ### - gamma_tp1: gamma_{t+1}, per trovare la partizione ridotta;
  ###               all'istante finale passare 'last time'
  ### - newclustevalue: valore per il nuovo cluster
  ### - spline_basis: matrice con i valori delle basi spline
  
  # Cluster in cui è l'i-esima osservazione
  # k <- which(vapply(rho_t, is.element, el = i, FUN.VALUE = FALSE))
  k <- lab_t[i]
  
  # Partizione senza i-esima osservazione
  # rho_tmi <- rho_t
  # rho_tmi[[k]] <- rho_tmi[[k]][rho_tmi[[k]] != i]
  lab_tmi <- lab_t
  lab_tmi[i] <- -1
  
  # Clusters with at least one observetion assigned to it; I need this to choose
  #   if its dCRP is cluster size or M
  non_empty_clust <- unique(lab_tmi[lab_tmi > 0])
  
  means <- drop(spline_basis[, -j] %*% beta_i[-j])
  means_list <- c()
  prob_cit_list <- c()
  check_list <- c()
  
  for (h in 1:nrow(beta_j)){
    # Devo costruire la partizione con l'i-esima osservazione inserita
    #   nell'h-esimo cluster
    lab_t.h <- lab_tmi
    lab_t.h[i] <- h
    
    # Since here we always have the beta's for all the clusters (even the empty 
    #   ones), if the unit i is assigned to an empty cluster, I don't sample a 
    #   new value for the beta because I already have it. I think this 
    #   procedure is more similar to MacEachern&Muller1994 than to Neal200.
    # I don't think this is a problem because at each iteration I update all
    #   betas, also the ones corresponding to empty clusters.
    beta <- beta_cluster[h]
    
    # Now the only difference between empty and non-empty clusters is the weight
    #   given by the CRP. If the cluster is not empty than its weight is the 
    #   size of the cluster, otherwise it's M divided by the number of empty 
    #   clusters at that moment
    if (h  %in% non_empty_clust){ # se il cluster è già esistente, devo usare il beta di quel cluster
      prob_cit <- sum(lab_tmi == h)
    } else { # If the cluster to which i is assigned is an empty one
      prob_cit <- M / (nrow(beta_j) - length(non_empty_clust))
      # The division is not actually necessary because it could go into the
      #   normalizing constant, but I'm explicitly including it so that I 
      #   keep track of the "true" procedure that I'm following
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
      
      lab_tp1.R <- rep(-1, length(lab_t))
      lab_t.h.R <- rep(-1, length(lab_t))
      
      lab_t.h.R[R_tp1] <- lab_t.h[R_tp1]
      lab_tp1.R[R_tp1] <- lab_tp1[R_tp1]
      
      check <- all(lab_t.h.R == lab_tp1.R)
    }
    
    means_list <- cbind(means_list, means + beta*spline_basis[, j])
    prob_cit_list <- c(prob_cit_list, prob_cit)
    check_list <- c(check_list, check)
    
  }
  
  
  logpnorm_list <- colSums(dnorm(Y_it, 
                                 mean = means_list,
                                 sd = sigma_i,
                                 log = TRUE))
  
  logp <- log(prob_cit_list) + logpnorm_list + log(check_list)
  
  # Almeno nelle prime iterazioni sembra che le prob siano tutte 0, 
  #   quindi devo lavorare in scala logaritmica. Per campionare dalle 
  #   log-prob uso il trick del max con la Gumbel.
  # lll <- logp + min(abs(min(logp)), .Machine$double.xmax)
  # gumbel <- -log(-log(runif(length(lll))))
  # lll <- lll + gumbel
  gumbel <- -log(-log(runif(length(logp))))
  lll <- logp + gumbel
  
  c_it <- whichclusters[which.max(lll)]
  return(c_it)
}






## ### ### ### ### ##
#### UPDATE BETA ####
## ### ### ### ### ##
up_beta <- function(Y.t = Y.tilde,
                    l2v = lik2var,
                    l2m = lik2mean,
                    obs = obs_k,
                    SIGinv = SIGMAinv,
                    d_j = delta_temp[[j]],
                    phi_j = phi_temp[[j]]){
  
  # Posterior precision matrix
  prec.post <- (SIGinv/d_j^2) + diag(l2v, ncol = ncol(SIGinv))
  # Posterior varcov matrix
  varcov.post <- mysolve(prec.post)
  # Media a posteriori
  mean.post <- varcov.post %*%  
    ( (SIGinv/d_j^2) %*% rep(phi_j, ncol(SIGinv)) + t(l2m)  )
  
  beta_updated <- drop(rmvn(1, 
                            mu = mean.post, 
                            sigma = varcov.post))
  
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





# ### ### ### ### ### ##
# #### UPDATE THETA ####
# ### ### ### ### ### ##
# up_theta_jt <- function(beta_jt,
#                         tau_jt,
#                         phi_j,
#                         delta_j){
#   ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#   ### - beta_jt: vector of beta_kjt for all clusters k
#   ### - the other parameters' names are clear
#   
#   x <- beta_jt[!is.na(beta_jt)]
#   xbar <- mean(x)
#   n <- length(x)
#   
#   postpar <- GaussGaussUpdate_iid(xbar = xbar,
#                                   n = n,
#                                   datavar = tau_jt^2,
#                                   priormean = phi_j,
#                                   priorvar = delta_j^2)
#   
#   out <- rnorm(1, 
#                mean = postpar[1], 
#                sd = sqrt(postpar[2]))
#   
#   return(out)
# }





### ### #### ### ### 
#### UPDATE PHI ####
### ### #### ### ### 
up_phi_j <- function(beta_j,
                     delta_j,
                     SIGinv,
                     lambda,
                     xi){
  ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### - beta_j: vector with beta_kjt for all k and all t
  ### - the other parameters' names are clear
  
  x <- beta_kjt
  K <- nrow(x)
  
  # Precision of the posterior
  prec.post <- (K/delta_j^2) * sum(SIGinv) + 1/xi^2
  # Variance of the posterior
  var.post <- 1/prec.post
  # Mean of the posterior
  mean.post <- (1/delta_j^2) * sum( colSums(SIGinv) * colSums(x) ) + 
    lambda/xi^2 
  
  out <- rnorm(1, 
               mean = mean.post, 
               sd = sqrt(var.post))
  
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
                                  datavar = xi^2,
                                  priormean = m0,
                                  priorvar = s02)
  
  out <- rnorm(1, postpar[1], sqrt(postpar[2]))
  
  return(out)
}





### ### ### ### ### ##
#### UPDATE DELTA ####
### ### ### ### ### ##
# I need a specific update for delta_j's because now they aren't the variances,
#   but the "scale" parameter of the varcov matrix
up_delta_j <- function(val_now,
                       beta_j,
                       phi_j,
                       SIGchol,
                       A_delta,
                       eps){
  ### ### ### ### ### ### ###\
  ### - beta_j: vector with beta_kjt for all clusters k and all times t
  # Sample proposed value
  val_prop <- runif(1, val_now - eps, val_now + eps)
  
  
  # Vector with prob's of acceptance
  
  # Indicator to check if proposed values are in the domain
  indicator <- (val_prop > 0 & val_prop < A_delta)
  
  if (indicator){
    # Likelihood of current value and proposed one
    likel_now <- sum(dmvn(beta_j, 
                      mu = rep(phi_j, ncol(beta_j)), 
                      sigma =  val_now * SIGchol,
                      log = TRUE,
                      isChol = TRUE))
    
    likel_prop <- sum(dmvn(beta_j, 
                           mu = rep(phi_j, ncol(beta_j)), 
                           sigma =  val_prop * SIGchol,
                           log = TRUE,
                           isChol = TRUE))
    logprob <- min(likel_prop - likel_now, 0)
  } else {
    logprob <- -Inf
  }
  
  # Create output
  acc <- (log(runif(1, 0, 1)) < logprob) 
  out <- ifelse(acc, val_prop, val_now)
  
  return(out)
}




### ### ### ### ### ### ##
#### UPDATE VARIANCES ####
### ### ### ### ### ### ##
# Function for the RW Metropolis (RWM) for tau, delta and xi
# It is written to update one param. at a time
up_sd.RWM <- function(val_now,
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
                           sd = val_now,
                           log = TRUE),
                     na.rm = TRUE) 
    likel_prop <- sum(dnorm(data, 
                            mean = mean,
                            sd = val_prop,
                            log = TRUE),
                      na.rm = TRUE)
    logprob <- min(likel_prop - likel_now, 0)
  } else {
    logprob <- -Inf
  }
  
  # Create output
  acc <- (log(runif(1, 0, 1)) < logprob) 
  out <- ifelse(acc, val_prop, val_now)
  
  return(out)
}





### ### ### ### ### ##
#### UPDATE ALPHA ####
### ### ### ### ### ##
up_alpha_j <- function(sum_data,
                       n,
                       priorshape1,
                       priorshape2){
  out <- rbeta(1, 
               priorshape1 + sum_data, 
               priorshape2 + n - sum_data)
  
  return(out)
}