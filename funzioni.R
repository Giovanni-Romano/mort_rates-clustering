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
                       b_it, b_itm1,
                       mustarvec_jt, muvec_t,
                       T_, Q,
                       rho_t, rho_tp1, gamma_tp1,
                       newclustermean) {
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### - i: indice dell'unità che stiamo aggiornando
  ### - j: indice del coefficiente che stiamo aggiornando
  ### - b_it: valore stato latente per unità i al tempo t
  ### - b_itm1: valore stato latente per unità i t-1
  ### - mustarvec_jt: vettore con le medie dei cluster per il coeff. j al tempo t
  ### - muvec_t: vettore con le medie degli errori dello stato latente al tempo t;
  ###             N.B.: i valori che vanno passati dipendono dal valore delle
  ###                   label al tempo t all'iterazione corrente
  ### - T_: matrice del sistema di equazioni dello stato latente (l'underscore
  ###       c'è solo perché non posso chiamarlo T per via di TRUE)
  ### - Q: matrice var-cov degli errori dello stato latente; per ora è fissata,
  ###       perciò non ha nessuna prior
  ### - rho_t: partizione al tempo t PER IL COEFF j
  ### - rho_tp1: partizione al tempo t+1 PER IL COEFF j
  ### - gamma_tp1: gamma_{t+1}, per trovare la partizione ridotta;
  ###               all'istante finale passare 'last time'
  ### - newclustermean: media per il nuovo cluster
  
  # Inizializzazione vettore probabilità della multinomiale
  logp <- c()
  
  # Cluster in cui è l'i-esima osservazione
  k <- which(vapply(rho_t, is.element, el = i, FUN.VALUE = FALSE))
  
  # Partizione senza i-esima osservazione
  rho_tmi <- rho_t
  rho_tmi[[k]] <- rho_tmi[[k]][rho_tmi[[k]] != i]
  
  whichcluster <- c(which(vapply(rho_tmi, function(x) length(x)>0, FUN.VALUE = FALSE)), 
                    length(rho_t)+1L)
  # Uso "1L" per farlo venire integer, altrimenti il "+1" lo rende un float
  
  for (h in whichcluster){
    # Devo costruire la partizione con l'i-esima osservazione inserita
    #   nell'h-esimo cluster
    rho_t.h <- rho_tmi
    if (h < max(whichcluster)){ # se assegno ad un cluster già esistente
      rho_t.h[[h]] <- c(rho_t.h[[h]], i) 
    } else { # se assegno al nuovo cluster
      rho_t.h[[h]] <- i 
    }
    
    # Nel paper di Page: Pr(c_it = h)
    # prob_cit <- dCRP(sapply(rho_t.h, length), M)
    prob_cit <- dCRP(vapply(rho_t.h, length, FUN.VALUE = integer(1)), M)
    
    # Nel paper di Page: N(Y_it | bla bla);
    #   io devo usare la densità della normale multivariata
    if (h < max(whichcluster)){ # se il cluster è già esistente, devo usare la media di quel cluster
      mu <- mustarvec_jt[h]
    } else { # se il cluster è nuovo, devo simulare la media dalla prior
      mu <- newclustermean
    }
    # La media che devo usare è il vettore delle medie mu_vec, sostituendo in
    #   per il j-esima coeff. la mu che ho appena calcolato.
    #
    muvec_t[1 + 3 * (j-1)] <- mu
    
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
      
    meantopass <- muvec_t + T_ %*% b_itm1
    logpnorm <- mvnfast::dmvn(b_it,
                              mu = meantopass, sigma = Q,
                              log = TRUE)
    
    logp <- c(logp, 
              log(prob_cit) + logpnorm + log(check))
  }
  
  
  # Almeno nelle prime iterazioni sembra che le prob siano tutte 0, 
  #   quindi devo lavorare in scala logaritmica. Per campionare dalle 
  #   log-prob uso il trick del max con la Gumbel.
  lll <- logp + min(abs(min(logp)), .Machine$double.xmax)
  gumbel <- -log(-log(runif(length(lll))))
  lll <- lll + gumbel
  
  c_it <- whichcluster[which.max(lll)]
  return(c_it)
}



#@@@@@@@@@@@@@@@@@@@@@@@#
#### UPDATE MUSTAR_t ####
#@@@@@@@@@@@@@@@@@@@@@@@#
up_mustar <- function(j, k, t,
                      b_t,  b_tm1,
                      mu_t, T_list, Q_list,
                      rho_t,
                      priorm, priorv){
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### - j: indice del coefficiente di cui stiamo aggiornando la media
  ### - k: indice del cluster di cui stiamo aggiornando la media
  ### - t: istante temporale attuale
  ### - b_t: matrice con gli stati latenti al tempo t per tutti i
  ###         coefficienti e tutte le unità; ogni riga è relativa ad 
  ###         un'osservazione e ogni colonna ad un coeff.;
  ###         non basta solo lo stato latente relativo al coefficiente j
  ###         perché i coefficienti sono correlati
  ### - b_tm1: come b_t ma al tempo t-1
  ### - mu_t: matrice delle medie degli errori dello stato latente al tempo t;
  ###             ogni riga si riferisce un'osservazione e ogni colonna ad un
  ###             coeff.
  ###             N.B: i valori che vanno passati dipendono dal valore delle
  ###                   label al tempo t all'iterazione corrente
  ### - T_list: lista con le matrici T del sistema di equaz. dello stato latente
  ### - Q_list: lista con le matrici var-cov degli errori dello stato latente
  ### - rho_t: partizione al tempo t PER IL COEFF j
  ### - priorm: media a priori per la mu dei cluster
  ### - priorv: varianza a priori per la mu dei cluster
  
  # Devo trovare quali elementi di b_jt appartengono al cluster k e poi uso
  #   classica conjugacy Gauss-Gauss.
  obs <- rho_t[[k]]
  
  # Lavoro su b_it - T*b_itm1, così che mustar sia proprio la media di questa
  #   nuova variabile. Inoltre tengo solo le osservazioni che servono.
  b.tilde <- c()
  if (t == 1){
    b.tilde <- b_t[obs, ]
  } else{
    for (i in obs){
      b.tilde <- rbind(b.tilde, t(b_t[i, ] - T_list[[i]] %*% b_tm1[i, ]))
    }
    # Se b.tilde ha una riga sola, lo trasformo in vettore; se ha più di 
    #   una riga, non succede niente.
    b.tilde <- drop(b.tilde)
  }
  
  mu <- mu_t[obs, ]
  # Se mu ha una riga sola, lo trasformo in vettore; se ha più di 
  #   una riga, non succede niente.
  mu <- drop(mu)
  
  # Devo togliere le Q delle osservazioni che non usiamo
  Q_list.new <- list()
  for (iii in obs){
    Q_list.new[[length(Q_list.new) + 1]] <- Q_list[[iii]]
  }
  
  n_k <- length(Q_list.new)
  p <- nrow(Q_list.new[[1]])
  if (n_k > 1){
    list_cond <- lapply(1:n_k, 
                        function(x) condMVN(mu[x, ], Q_list.new[[x]], j, 
                                            (1:p)[-j], b.tilde[x, -j]))
    
    cvar <- vapply(list_cond, function(x) x$condVar, FUN.VALUE = double(1))
    cmean <- vapply(list_cond, function(x) x$condMean, FUN.VALUE = double(1))
    
    post.var <- 1 / (1/priorv + sum(1/cvar))
    post.mean <- post.var * (priorm/priorv + sum(b.tilde[,j]/cvar))
  } else {
    list_cond <- condMVN(mu, Q_list.new[[1]], j, (1:p)[-j], b.tilde[-j])
    
    cvar <- list_cond$condVar
    cmean <- list_cond$condMean
    
    post.var <- 1 / (1/priorv + sum(1/cvar))
    post.mean <- post.var * (priorm/priorv + b.tilde[j]/cvar)
  }
  
  rnorm(1, mean = post.mean, sd = sqrt(post.var))
}
