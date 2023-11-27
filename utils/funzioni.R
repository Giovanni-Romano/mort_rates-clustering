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
  require(mggd)
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





### ### ### ### ### ### ### ### ###
##### Similarity function g_3() ####
### ### ### ### ### ### ### ### ###
# Firstly, I define the function for the univariate case (i.e. just one covariate)
g3_univ <- function(s2, x, v02, mu0){
  n <- length(x)
  out <- (1/2/pi)^(n/2) * (1/s2)^((n-1)/2) * sqrt(1/(s2+n*v02)) * 
    exp( -(1/2/s2/(s2+n*v02)) * 
           ( (s2+n*v02)*sum(x^2) - v02*(sum(x)^2) + 
               s2*n*mu0^2 - 2*s2*mu0*sum(x)) 
    )
  
  # It's faster written in this way than with any R function that I tried; 
  #   probably because in my formula there is no inversion, since I know 
  #   analytically found it.
  
  return(out)
}

# I define the function for the multivariate case (i.e. more than one covariate) 
#   as the product of the univariate ones.
g3_multiv <- function(s2_vec, X, v02_vec, mu0_vec){
  # X is a matrix whose columns are the covariates for all the n units
  
  # If X has just 1 row, then R does not consider him as a matrix, but as a 
  #   vector. So, I have to force it to be a matrix.
  if (is.null(dim(X))) {X <- matrix(X, nrow = 1)}
  tmp <- sapply(1:ncol(X), 
                function(j) g3_univ(s2_vec[j], X[,j], v02_vec[j], mu0_vec[j]))
  
  out <- prod(tmp)
  
  return(out)
}



### ### ### ### ### ##
#### UPDATE GAMMA ####
### ### ### ### ### ##
up_gamma_i <- function(i, gamma,
                       alpha_t, 
                       lab_t, lab_tm1,
                       X_t, # matrix of covariates at time t (for all units)
                       s2_vec, # vector of variances of covariates in g3
                       mu0_vec, # vector of means of the priors in g3
                       v02_vec){ # vector of variances of the priors in g3
  
  # R_t^{(-i)}: set of units that don't change cluster from t-1 to t with {i} removed --> Rtmi
  # R_t^{(+i)}: set of units that don't change cluster from t-1 to t united to {i} --> Rtpi
  R_t <- which(gamma == 1)
  R_tmi <- R_t[R_t != i]
  R_tpi <- c(R_tmi, i)
  
  # I use ".R" for the reduced partitions.
  # "pi" stands for "plus i" and "mi" for "minus i".
  lab_tm1.Rpi <- rep(-1, length(lab_t)) # it's used to compute the indicator in formula (S.8) of suppl. mat.
  lab_tm1.Rpi[R_tpi] <- lab_tm1[R_tpi]
  lab_t.Rpi <- rep(-1, length(lab_t)) # it's used to compute the indicator in formula (S.8) of suppl. mat.
  lab_t.Rpi[R_tpi] <- lab_t[R_tpi]
  # I need also the reduced partitions w.r.t. R_{t}^{(-i)} to compute 
  #   the unnormalized prob. relative to formula (S.9) in suppl. mat.
  lab_t.Rmi <- rep(-1, length(lab_t))
  lab_t.Rmi[R_tmi] <- lab_t[R_tmi]
  check <- all(lab_tm1.Rpi == lab_t.Rpi)
  
  
  # If the reduced partitions (and hence the corresponding indicator is 1), 
  #   then I have to compute the rest of the probability.
  # On the other hand if the reduced partitions are different, then the 
  #   corresponding indicator is 0 and hence all the prob. is 0 and I don't 
  #   have to compute anything else.
  if (check){ 
    
    # See formula S.9 Page's supplem. material.
    # I have to compute the predictive distrib. for the i-th observation given
    #   the reduced partition w.r.t. R_{t}^{(-i)}.
    # In this case (w/ the PPMx) I don't have the closed form expression for the
    #   predictive, hence I have to compute its kernel for all possible 
    #   assignments of the i-th unit to all clusters.
    
    if (length(R_tmi) == 0) { # if the partition to wich we condition is empty
      ratio <- 1
    } else {
      # Find the non-empty clusters in the reduced partition
      non_empty_clust <- unique(lab_t.Rmi[lab_t.Rmi > 0])
      
      # Find the subset j of the partition w/ i-th unit inside
      j <- lab_t[i]
      
      # Initialize the vector of unnormalized probabilities
      unnorm_prob_list <- rep(NA, length(gamma))
      # Fill the vector of unnormalized probabilities for the empty clusters.
      #   These probabilities are all equal and they are equal to what written 
      #   in formula (S.9) of Page's suppl. mat. divided by the nuyber of empty 
      #   clusters (since we always have all "tables" open).
      unnorm_prob_list[-non_empty_clust] <- 
        M / (length(gamma) - length(non_empty_clust)) * 
        g3_multiv(s2_vec, X_t[i, ], v02_vec, mu0_vec)
      
      for (h in non_empty_clust){
        idx_h <- which(lab_t.Rmi == h)
        n_h <- length(idx_h)
        unnorm_prob <- n_h * 
          g3_multiv(s2_vec, X_t[c(idx_h, i), ], v02_vec, mu0_vec) / 
          g3_multiv(s2_vec, X_t[idx_h, ], v02_vec, mu0_vec)
          
        unnorm_prob_list[h] <- unnorm_prob
      }
      
      norm_prob_list <- unnorm_prob_list / sum(unnorm_prob_list)
      ratio <- norm_prob_list[j]

    }
    
    prob <- alpha_t / (alpha_t + (1 - alpha_t) * ratio)
    gamma_it <- as.integer(runif(1) < prob)
  } else {
    prob <- 0 # not necessary, here just for conceptual clarity
    gamma_it <- 0
  }
  
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
                       spline_basis,
                       X_t, # matrix of covariates at time t (for all units)
                       s2_vec, # vector of variances of covariates in g3
                       mu0_vec, # vector of means of the priors in g3
                       v02_vec){ # vector of variances of the priors in g3)
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
  beta_list <- c()
  prob_cit_list <- c()
  check_list <- c()
  
  # I save this object because I will use it many times
  spl_bas_j <- spline_basis[, j]
  
  # I create these objects here, but I will use it inside the loop.
  #   Up to now I used to build them inside the loop, but it is not necessary
  #   since it these quantities do NOT depend on "h" in any way.
  if (!(identical(gamma_tp1, 'last time'))){
    lab_tp1.R <- rep(-1, length(lab_t))
    R_tp1 <- which(gamma_tp1 == 1)
    lab_tp1.R[R_tp1] <- lab_tp1[R_tp1]
  }
  
  
  # I compute gamma(|S^{-i}|) g(x_{kjt}^{-i \star}) for all clusters S, where 
  #   I set it to 1 if S is empty 
  gamma_times_g_list <- rep(1, nrow(X_t))
  for (h in non_empty_clust){
    gamma_times_g_list[h] <- gamma(sum(lab_tmi == h)) * 
                                    g3_multiv(s2_vec, 
                                              X_t[lab_tmi == h, ], 
                                              v02_vec, mu0_vec)
  }
  
  
  # I compute outside the probability of the i-th unit to be assigned to an 
  #   empty cluster because it is the same for all the empty clusters and it 
  #   does not depend on "h"
  prob_empty_clust <- M / (length(beta_cluster) - length(non_empty_clust)) * 
    g3_multiv(s2_vec, 
              X_t[i, ], 
              v02_vec, mu0_vec) *
    prod(gamma_times_g_list)
  

  for (h in 1:length(beta_cluster)){
    # Devo costruire la partizione con l'i-esima osservazione inserita
    #   nell'h-esimo cluster
    lab_t.h <- lab_tmi
    lab_t.h[i] <- h
    
    
    # Nel paper di Page: indicatrice sulle partizioni ridotte
    #   per prima cosa devo trovare la partizione ridotta
    if (identical(gamma_tp1, 'last time')){
      check <- TRUE
    } else {
      
      lab_t.h.R <- rep(-1, length(lab_t))
      lab_t.h.R[R_tp1] <- lab_t.h[R_tp1]
      
      check <- all(lab_t.h.R == lab_tp1.R)
    }
    
    
    # Since here we always have the beta's for all the clusters (even the empty 
    #   ones), if the unit i is assigned to an empty cluster, I don't sample a 
    #   new value for the beta because I already have it. I think this 
    #   procedure is more similar to MacEachern&Muller1994 than to Neal200.
    # I don't think this is a problem because at each iteration of the Gibbs 
    # then I update all betas, also the ones corresponding to empty clusters.
    beta <- beta_cluster[h]
    
    
    # I compute prob_cit only if check == TRUE, because otherwise it is useless,
    #   since it will be multuplied by the indicator corrisponding to check, 
    #   which is 0.
    if (check){ # If the indicator (check) is equal to 1..
      if (h  %in% non_empty_clust){ # se il cluster è già esistente, devo usare il beta di quel cluster
        prob_cit <- M * gamma(sum(lab_tmi == h) + 1) * 
          g3_multiv(s2_vec, 
                    X_t[c(which(lab_tmi == h), i), ], 
                    v02_vec, mu0_vec) *
          prod(gamma_times_g_list[-h])
      } else { # If the cluster to which i is assigned is an empty one
        prob_cit <- prob_empty_clust
        # Now I compute this prob just once before the for loop because it does
        #   NOT depend on "h".
      }
    } else { # .. else if the indicator (check) is equal to 0..
      # I put it to 0, but it does not matter the value because then it will be
      #   multiplied by as.integer(check) that is 0.
      prob_cit <- 0
    }
    
    # No need to build beta_list, it will coincide with beta_cluster..
    #   For the moment I leave it in the code for clarity.
    beta_list <- c(beta_list, beta)
    prob_cit_list <- c(prob_cit_list, prob_cit)
    check_list <- c(check_list, check)
    
  }
  
  
  logpnorm_list <- colSums(dnorm(Y_it,
                                 mean = means + outer(spl_bas_j, beta_list),
                                 sd = sigma_i,
                                 log = TRUE))
  
  logp <- log(prob_cit_list) + logpnorm_list + log(check_list)
  
  # Almeno nelle prime iterazioni sembra che le prob siano tutte 0, 
  #   quindi devo lavorare in scala logaritmica. Per campionare dalle 
  #   log-prob uso il trick del max con la Gumbel.
  gumbel <- -log(-log(runif(length(logp))))
  lll <- logp + gumbel
  
  c_it <- which.max(lll)
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
  
  x <- beta_j
  K <- nrow(x)
  
  # Precision of the posterior
  prec.post <- (K/delta_j^2) * sum(SIGinv) + 1/xi^2
  # Variance of the posterior
  var.post <- 1/prec.post
  # Mean of the posterior
  mean.post <- var.post * ( (1/delta_j^2) * sum( colSums(SIGinv) * colSums(x) ) + 
                              lambda/xi^2 )
  
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

# I noticed that the beahviour of delta in the model w/ time-dependence in the
#   beta's needed a different treatment with respect to the other variance's
#   parameters. A Uniform RW-MH with step-size equal to the length of the 
#   domain (eps = A_delta) was not okay because the likelihood is quite peaked
#   on an interval of length 0.2/0.3 and so the majority of proposed values 
#   were rejected. The solution that I am using is to reduce the eps and to 
#   use a Gaussian RW to have more proposed values close to the current point
#   but w/o completely avoid longer steps (as it would happen with a Unif RW
#   with eps = 0.1).

up_delta_j <- function(val_now,
                       beta_j,
                       phi_j,
                       SIGchol,
                       A_delta,
                       eps){
  ### ### ### ### ### ### ###
  ### - beta_j: vector with beta_kjt for all clusters k and all times t
  # Sample proposed value
  val_prop <- rnorm(1, mean = val_now, sd = eps)
  
  
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





### ### ### ### ### ### ### ###
#### UPDATE S2 (param. g3) ####
### ### ### ### ### ### ### ###
# To compute the acceptance probability I need the likelihood of the PPMx for
#   s2 (it is not sufficient the kernel of the distribution of the PPM, since 
#   s2 will probably appear also in the normalizing constant).
# I'm not able to normalize analitycally the likelihood, so I should do that 
#   numerically. It is usually unfeasible because the number of possible 
#   partitions increases as n^n, but in the simulation study n is just 3, so 
#   it is feasible.
unnorm_ppmx <- function(labs,
                        M,
                        # g3_multiv. parameters
                        s2_vec, 
                        X, 
                        v02_vec, 
                        mu0_vec){
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### - labs: vector with the labels of the clusters for which we want to 
  ###         compute the unnorm_ppmx
  ### - M: concentration parameter of the Dirichlet process
  ### - s2_vec: vector of variances of the marginal likelihood
  ### - X: matrix with the covariates for each unit i
  ### - mu0_vec: vector of means of the marginal likelihood
  ### - v02_vec: vector of variances of the distribution used to marginalize 
  ###             the likelihood
  
  rho <- lab2rho(labs)
  
  # I need labs to be w/o gaps, so I can't have a label vector like c(1, 3, 1),
  #   but it should be c(1, 2, 1).
  rho <- rho[lapply(rho, length) > 0]
  
  
  tmp <- sapply(rho, function(r) M * gamma(length(r)) * g3_multiv(s2_vec, 
                                                                  X[r, ], 
                                                                  v02_vec, 
                                                                  mu0_vec))
  
  out <- prod(tmp)
  
  return(out)
}

norm_ppmx <- function(labs,
                      M,
                      # g3_multiv. parameters
                      s2_vec, 
                      X, 
                      v02_vec, 
                      mu0_vec){
  require(partitions)
  
  n <- length(labs)
  
  parts <- cbind(labs, setparts(n))
  
  # Compute the likelihood for each partition
  out <- apply(parts, 2, function(p) unnorm_ppmx(p,
                                                 M,
                                                 s2_vec,
                                                 X,
                                                 v02_vec,
                                                 mu0_vec))
  
  return(out[1]/sum(out[-1]))
}

# Function for the RW Metropolis (RWM) for the parameter s2 of the similarity 
#   function g3(), that is the variance of the covariates.
# It is written to update one param. at a time.
up_s2.RWM <- function(s2_vec,
                      idx_cov,
                      eps,
                      A_s2,
                      labels,
                      M,
                      X,
                      mu0_vec,
                      v02_vec){
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### - s2_vec: vector of variances of the likelihood
  ### - idx_cov: index of the covariate for which we want to update s2
  ### - eps: RW step size
  ### - A_s2: hyperparam of prior on s2
  ### - labels: matrix with the labels of the clusters for each time t and unit i
  ### - M: concentration parameter of the Dirichlet process
  ### - X: 3D-array with the covariates for each time t and unit i
  ### - mu0_vec: vector of means of the marginal likelihood
  ### - v02_vec: vector of variances of the distribution used to marginalize 
  ###             the likelihood
  
  # The current value for the s2 to update is in the idx_cov-th position of s2_vec
  val_now <- s2_vec[idx_cov]
  
  # Sample proposed value
  val_prop <- runif(1, val_now - eps, val_now + eps)
  
  # Indicator to check if proposed values are in the domain
  indicator <- (val_prop > 0 & val_prop < A_s2)
  
  if (indicator){
    # Likelihood of current value
    s2_vec_now <- s2_vec # no need to create s2_vec_now; here just for clarity
    likel_now <- sum(vapply(1:10, 
                            function(yr) log(norm_ppmx(labs = labels[, yr],
                                                       M = M,
                                                       s2_vec = s2_vec_now, 
                                                       X = X[ , yr, ], 
                                                       v02_vec = v02_vec, 
                                                       mu0_vec = mu0_vec)),
                            FUN.VALUE = numeric(1))
    )
    
    # Likelihood of the proposed value
    s2_vec_prop <- s2_vec; s2_vec_prop[idx_cov] <- val_prop
    likel_prop <- sum(vapply(1:10, 
                             function(yr) log(norm_ppmx(labs = labels[, yr],
                                                        M = M,
                                                        s2_vec = s2_vec_prop,
                                                        X = X[, yr, ], 
                                                        v02_vec = v02_vec, 
                                                        mu0_vec = mu0_vec)),
                             FUN.VALUE = numeric(1))
    )
    logprob <- min(likel_prop - likel_now, 0)
  } else {
    logprob <- -Inf
  }
  
  # Create output
  acc <- (log(runif(1, 0, 1)) < logprob) 
  out <- ifelse(acc, val_prop, val_now)
  
  return(out)
}