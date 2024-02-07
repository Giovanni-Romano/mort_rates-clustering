# library(tidyverse)
# library(mvtnorm)
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
  lab_t.R <- rep(-1, length(lab_t))
  lab_t.R[R_tpi] <- lab_t[R_tpi]
  
  adj_tm1 <- outer(lab_tm1[R_tpi], lab_tm1[R_tpi], function(x, y) as.integer(x == y))
  adj_t <- outer(lab_t[R_tpi], lab_t[R_tpi], function(x, y) as.integer(x == y))
  check <- all(adj_tm1 == adj_t)
  
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
                       newclustervalue,
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
  
  whichclusters <- unique(lab_tmi[lab_tmi > 0])
  whichclusters <- c(whichclusters, max(whichclusters)+1L) # I add the possible new cluster
  
  means <- drop(spline_basis[, -j] %*% beta_i[-j])
  beta_list <- c()
  # tmp_means_list <- c()
  prob_cit_list <- c()
  check_list <- c()
  
  # I save this object because I will use it many times
  spl_bas_j <- spline_basis[, j]
  
  
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
    adj_t <- outer(lab_t[R_tp1_mi], lab_t[R_tp1_mi], function(x, y) as.integer(x == y))
    adj_tp1 <- outer(lab_tp1, lab_tp1, function(x, y) as.integer(x == y))
    check_tmp <- all(adj_t == adj_tp1[R_tp1_mi, R_tp1_mi])
  }
  
  
  
  for (h in whichclusters){
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
        lab_t.h <- lab_t
        lab_t.h[i] <- h
        adj_t.h <- outer(lab_t.h[R_tp1], lab_t.h[R_tp1], function(x, y) as.integer(x == y))
        check <- all(adj_t.h == adj_tp1[R_tp1, R_tp1])
      } else {
        cat("check_tmp FALSE \n")
        check <- FALSE
      }
      
    }
    
    
    
    # If check is false, then it's useless to compute prob_cit because the 
    #   total probability will however be 0. So I set prob_cit = 0 (but every 
    #   constant would be ok, because it will be multiplied by 0).
    if (check){
      if (h %in% non_empty_clust){ # se il cluster è già esistente, devo usare il beta di quel cluster
        beta <- beta_cluster[h]
        prob_cit <- sum(lab_tmi == h)
      } else { # If the cluster to which i is assigned is an empty one..
        # The procedure in Neal's paper is a bit different from the one used in Page: 
        #   - Page always sample a new value from the prior
        #   - Neal sample a new value only if the cluster is actually new, i.e.
        #       if obs. i was a singleton before the update and now it is assigned
        #       to an empty cluster, than that cluster is its old one, so the 
        #       value of beta is the old one.
        if (sum(lab_t == k) == 1){
          beta <- beta_cluster[k]
        } else{
          beta <- newclustervalue
        }
        prob_cit <- M
      }
      
      # I beta che devo usare sono il vettore beta_i, sostituendo
      #   per il j-esimo coeff. il beta che ho appena calcolato.
      #
    } else {
      prob_cit <- 0
      beta <- 0
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
  
  logpnorm_list <- colSums(-0.5*log(2*pi*sigma_i^2) - 0.5/sigma_i^2 * (Y_it - means - outer(spl_bas_j, beta_list))^2)
  
  logp <- log(prob_cit_list) + logpnorm_list + log(check_list)
  
  # Almeno nelle prime iterazioni sembra che le prob siano tutte 0, 
  #   quindi devo lavorare in scala logaritmica. Per campionare dalle 
  #   log-prob uso il trick del max con la Gumbel.
  gumbel <- -log(-log(runif(length(logp))))
  lll <- logp + gumbel
  
  c_it <- whichclusters[which.max(lll)]
  return(c_it)
}





## ### ### ### ### ##
#### UPDATE BETA ####
## ### ### ### ### ##
up_beta <- function(j, k, t,
                    Y_t,
                    beta_t,
                    phi_jt, delta_j,
                    sigma_vec,
                    labels_t,
                    spline_basis){
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### - j: indice del coefficiente di cui stiamo aggiornando la media
  ### - k: indice del cluster di cui stiamo aggiornando la media
  ### - t: istante temporale attuale
  ### - Y_t: matrix with Y_{ixt} for all obs. i and all obs. x
  ### - beta_t: matrice con i beta per tutti i cluster per tutti i coeff. 
  ###           al tempo t
  ### - theta_jt, delta_j: media e sd della prior dei beta
  ### - sigma_vec: vettore con le sd delle osservazioni
  ### - labels_t: matrice delle label al tempo t per tutti i coeff.
  ### - spline_basis: matrice con i valori delle basi spline
  
  n <- length(sigma_vec)
  
  # Devo trovare quali osservazioni appartengono al cluster k per il 
  #   e coeff. j e poi uso conjugacy Gauss-Gauss.
  # obs <- lab2rho(labels_t[ , j])[[k]]
  obs <- which(labels_t[ , j] == k)
  
  # Lavoro su \tilde{Y}_ixt, cioè le osservazioni Y_{ixt} meno le spline tranne
  #   la j-esima, così che beta_jkt*g_j(x) sia la media di questa
  #   nuova variabile. Inoltre tengo solo le osservazioni che servono, cioè 
  #   quelle che per il coefficiente j appartengono al cluster k
  
  # Devo recuperare i beta giusti per ogni osservazione per ogni coeff. != j
  p <- ncol(spline_basis)
  beta_minusj_units_clusterk <- vapply(1:p, 
                                       function(x) beta_t[labels_t[ , x], x],
                                       FUN.VALUE = double(n))[obs, -j] # tolgo la j-esima colonna
  # e le oss. non nel cluster
  # Costruisco effettivamente \tilde{Y}_t
  Y_t.tilde <- Y_t[obs, ] - beta_minusj_units_clusterk %*% t(spline_basis[ , -j])
  # Se Y_t.tilde ha una riga sola, lo trasformo in vettore; se ha più di 
  #   una riga, non succede niente.
  Y_t.tilde <- drop(Y_t.tilde)
  
  
  # Devo togliere le sigma_i delle osservazioni che non usiamo
  sigma_clusterk <- sigma_vec[obs]
  
  
  # Varianza a posteriori
  var.post <- ( 1 / delta_j^2 + sum(1/sigma_clusterk^2)*sum(spline_basis[, j]^2) )^(-1)
  # Media a posteriori
  mean.post <- var.post * 
    ( sum( t(t(Y_t.tilde)*spline_basis[ , j])/sigma_clusterk^2 ) + 
        phi_jt / delta_j^2 )
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