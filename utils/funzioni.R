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
up_gamma_i <- function(i, gamma, alpha_t, lab_t, lab_tm1, M){
  
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
  
  # If check is false it's useless to do all computations, because the 
  #   probability will be always 0.
  if (check){
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
    
  } else {
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
                       M) {
  
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
  
  n <- length(lab_t)
  n_cl <- length(beta_cluster)
  
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
  # tmp_means_list <- c()
  prob_cit_list <- c()
  check_list <- c()
  
  # I save this object because I will use it many times
  spl_bas_j <- spline_basis[, j]
  
  # I create the following objects here, but I will use it inside the loop.
  #   Up to now I used to build them inside the loop, but it is not necessary
  #   since it these quantities do NOT depend on "h" in any way.
  prob_empty_clust <- M / (n_cl - length(non_empty_clust))
  
  
  # I need to compute the indicator for the equality of the reduced partitions
  #   at time t and t+1 with i assigned to the h-th cluster (see formula S.10 
  #   of the suppl. material).
  is.last_time <- identical(gamma_tp1, 'last time')
  if (!(is.last_time)){
    # \mathcal{R}_{t+1}
    R_tp1 <- which(gamma_tp1 == 1)
    
    # I remove the i-th unit because its check will be computed inside the loop
    R_tp1_mi <- R_tp1[R_tp1 != i]
    
    # Check for all units but the i-th one.
    check_tmp <- all(lab_tmi[R_tp1_mi] == lab_tp1[R_tp1_mi])
    
    
    # Here I compute the vector of the indicator's values for the i-th unit for 
    #   the assignment of "i" to the different clusters. 
    # I will use it inside the loop.
    # I check if the i-th unit is a free one (!gamma_tp1[i]) and, if it 
    #   is not, then I check if the proposal "i to the h-th cluster" is valid.
    check_i <- !gamma_tp1[i] | (lab_tp1[i] == (1:n_cl))
  }
  
  
  
  for (h in 1:n_cl){
    # Devo costruire la partizione con l'i-esima osservazione inserita
    #   nell'h-esimo cluster
    
    # Here I finish the computation of the indicator started before the loop.
    #   If we are at the last time, then I set check = TRUE because there aren't
    #   other partitions to compare with.
    #   If we are not at the last time, then I check if the i-th unit is a 
    #   constrained one and, if it is, then I check if the proposal "i to the 
    #   h-th cluster" is valid.
    if (is.last_time){
      check <- TRUE
    } else {
      if (check_tmp){
        check <- check_tmp & check_i[h]
      } else {
        check <- FALSE
      }
      
    }
    
    # Since here we always have the beta's for all the clusters (even the empty 
    #   ones), if the unit i is assigned to an empty cluster, I don't sample a 
    #   new value for the beta because I already have it. I think this 
    #   procedure is more similar to MacEachern&Muller1994 than to Neal200.
    # I don't think this is a problem because at each iteration I update all
    #   betas, also the ones corresponding to empty clusters.
    beta <- beta_cluster[h]
    
    
    # If check is false, then it's useless to compute prob_cit because the 
    #   total probability will however be 0. So I set prob_cit = 0 (but every 
    #   constant would be ok, because it will be multiplied by 0).
    if (check){
      # Now the only difference between empty and non-empty clusters is the weight
      #   given by the CRP. If the cluster is not empty than its weight is the 
      #   size of the cluster, otherwise it's M divided by the number of empty 
      #   clusters at that moment
      if (h  %in% non_empty_clust){
        prob_cit <- sum(lab_tmi == h)
      } else { # If the cluster to which i is assigned is an empty one
        prob_cit <- prob_empty_clust # M / (length(beta_cluster) - length(non_empty_clust))
        # Now I compute this prob just once before the for loop because it does
        #   NOT depend on "h".
        # The division is not actually necessary because it could go into the
        #   normalizing constant, but I'm explicitly including it so that I 
        #   keep track of the "true" procedure that I'm following
      }
    } else {
      prob_cit <- 0
    }
    
    # means_list <- cbind(means_list, means + beta*spline_basis[, j])
    # I change the part concerning the means so that I:
    #   - call spline_basis[, j] just once, before the for loop
    #   - I do the sum of means + beta... just once, after the for loop
    beta_list <- c(beta_list, beta)
    # tmp_means_list <- cbind(tmp_means_list, beta*spl_bas_j)
    prob_cit_list <- c(prob_cit_list, prob_cit)
    check_list <- c(check_list, check)
    
  }
  
  
  # logpnorm_list <- colSums(dnorm(Y_it, 
  #                                mean = means + outer(spl_bas_j, beta_list),
  #                                # "means" is sum_{h != j} beta_{h{c_iht}t} * g(x) for every x = 0, ..., 98
  #                                # outer(spl_bas_j, beta_list) is a matrix with 
  #                                #  r-th column equal to beta_r * g(x) for every x = 0, ..., 98.
  #                                # So the sum of the two is a matrix with number
  #                                #  of columns equal to the number of clusters
  #                                #  and number of rows equal to the number of ages,
  #                                #  where the r-th column is equal to the mean 
  #                                #  of the Gaussian for Y_it if "i" is assigned
  #                                #  to the r-th cluster.
  #                                # In this way the log-density of the normal is
  #                                #  computed for Y_it for every possible assignment
  #                                #  of "i" to a cluster. This exploits the fact
  #                                #  that R automatically recycles the vector Y_it
  #                                #  to have the same size of the matrix passed to
  #                                #  the argument "mean".
  #                                # To check this, try to replace "Y_it" with 
  #                                #  matrix(Y_it, nrow = length(Y_it), ncol = length(beta_cluster)).
  #                                sd = sigma_i,
  #                                log = TRUE))
  
  logpnorm_list <- colSums(-0.5*log(2*pi*sigma_i^2) - 0.5/sigma_i^2 * (Y_it - means - outer(spl_bas_j, beta_list))^2)
  
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
up_beta <- function(Y.t,
                    l2v,
                    l2m,
                    obs,
                    SIGinv,
                    d_j,
                    phi_j){
  
  # Posterior precision matrix
  prec.post <- (SIGinv/d_j^2) + diag(l2v, ncol = ncol(SIGinv))
  # Posterior varcov matrix
  varcov.post <- mysolve(prec.post)
  # Media a posteriori
  mean.post <- varcov.post %*%  
    ( (SIGinv/d_j^2) %*% phi_j + t(l2m)  )
  
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