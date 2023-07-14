rm(list=ls())

library(tidyverse)
library(splines2)
library(tidyverse)
library(invgamma)

load("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/BSP_Pavone/output/mortality.Rdata")
source("funzioni.R")

data_list_man <- list(ita_man = Y_ita_man / N_ita_man,
                      swe_man = Y_swe_man / N_swe_man,
                      uk_man = Y_uk_man / N_uk_man,
                      us_man = Y_us_man / N_us_man)

data_list_woman <- list(ita_woman = Y_ita_woman / N_ita_woman,
                        swe_woman = Y_swe_woman / N_swe_woman,
                        uk_woman = Y_uk_woman / N_uk_woman,
                        us_woman = Y_us_woman / N_us_woman)

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
# Parametro Bernoulli per i gamma
alpha <- 0.5
# Iperparametri della prior dell'ultimo layer
m0 <- 0; s02 <- 100 
# Iperparametri Inverse-Gamma
#   Nel paper di Page le prior sulle varianze sono Uniformi tra 0 e 5 o 10;
#     impostando i seguenti iperparametri ottengo una prior molto piatta, con 
#     peak tra 0 e 20: la moda è in ~5 e la varianza è "infinita" (l'InvGamma
#     ha varianza indefinita per alfa <= 2).
a_tau <- a_delta <- a_xi <- 1
b_tau <- b_delta <- b_xi <- 5
# Parametro di concentrazione del CRP
M <- 1





### ### ### ### ### ### ### ### ### ###
### Definizione parametri del Gibbs ###
### ### ### ### ### ### ### ### ### ###

# Numero iterazioni
n_iter <- 10





### ### ### ### ### ### ### ### ### ### ###
### Creo gli oggetti in cui salverò gli ###
### output dei vari step del Gibbs      ###
### ### ### ### ### ### ### ### ### ### ###

# Ho provato ad inizializzare questi oggetti con 99999 sperando di velocizzare,
#   ma non cambia niente. Perciò, preferisco tenere NA, almeno se qualcosa va
#   storto me ne accorgo.
gamma_res <- labels_res <- 
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
phi_res <- delta_res <- array(NA,
                              dim = c(T_final, # numero istanti temporali
                                      n_iter # numero iterazioni
                                      )
                              )
lambda_res <- xi_res <- rep(NA, n_iter)



gamma_temp <- replicate(p,
                        matrix(0,
                               nrow = n, # numero di osservazioni
                               ncol = T_final # numero istanti temporali
                        ),
                        simplify = FALSE)
labels_temp <- replicate(p,
                         matrix(NA,
                                nrow = n, # numero di osservazioni
                                ncol = T_final # numero istanti temporali
                         ),
                         simplify = FALSE)




### ### ### ### ### ###
### Inizializzazioni ##
### ### ### ### ### ###
lambda_res[1] <- rnorm(1, mean = m0, sd = sqrt(s02))
xi_res[1] <- rinvgamma(1, a_xi, b_xi) 

phi_res[, 1] <- rnorm(T_final, lambda_res[1], xi_res[1])
delta_res[, 1] <- rinvgamma(T_final, a_delta, b_delta)

for (j in 1:p){
  theta_res[[j]][ , 1] <- rnorm(T_final, phi_res[, 1], delta_res[, 1])
  tau_res[[j]][ , 1] <- rinvgamma(T_final, a_tau, b_tau)
  for (t in 1:T_final){
    lab <- rho2lab(rCRP(4, 1))
    labels_temp[[j]][, t] <- labels_res[[j]][, t, 1] <- lab
    beta_res[[j]][1:max(lab), t, 1] <- rnorm(max(lab),
                                             mean = theta_res[[j]][t, 1], sd = tau_res[[j]][t, 1])
  }
  
  # Passo a gamma_res[[j]][ , , d] il gamma_temp che ho inizializzato
  gamma_res[[j]][ , , 1] <- gamma_temp[[j]]
}

rm(j); rm(t)



f <- function(){
  for (d in 2:n_iter){ # Ciclo sulle iterazioni
    
    if ((d %% 100) == 0) {cat(d, "\n")}
    
    for (j in 1:p){ # Ciclo sui coefficienti
      for (t in 1:T_final){ # Ciclo sugli istanti
        ### ### ### ### ### ### ####
        ### ### UPDATE GAMMA ### ###
        # Recupero le partizioni al tempo t, t-1 e t+1.
        #   Lo faccio qui perché è uguale per tutte le osservazioni i e
        #   lo userò anche dopo per l'update delle labels.
        # Se t==1 non mi serve rho_tm1 perché non devo fare il confronto
        #   di compatibilità.
        if (t > 1){
          rho_tm1 <- lab2rho(labels_res[[j]][, t-1, d])
        }
        rho_t <- lab2rho(labels_res[[j]][, t, d-1])
        if (t < T_final){
          rho_tp1 <- lab2rho(labels_res[[j]][, t+1, d-1])
        }
        
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
                                                alpha_t = alpha,
                                                rho_t = rho_t,
                                                rho_tm1 = rho_tm1)
          }
        } # Fine ciclo sulle osservazioni per gamma
        gamma_res[[j]][, t, d] <- gamma_temp[[j]][, t]
        # Forse potrei farlo una volta sola, alla fine del ciclo sulle "t".
        
        
        
        ### ### ### ### ### ### ### #
        ### ### UPDATE LABELS ### ###
        # Per la partizione al tempo t+1 posso usare quella calcolata
        #   prima dell'aggiornamento dei gamma (tanto si usano quelle
        #   dell'iterazione d-1).
        # Per la partizione al tempo t devo creare un oggetto temporaneo 
        #   perché devo aggiornarla ogni volta che la label di 
        #   un'osservazione "i" viene aggiornata
        labels_temp <- labels_res[[j]][,t,d-1]
        
        for (i in 1:n){ # Ciclo sulle osservazioni per gamma
          # Per ogni coefficiente trovo la label a cui è assegnata l'unità i
          #   al tempo t all'iterazione più recente. Quindi per i coefficienti
          #   con indice <j possiamo prendere all'iteraz. d perché abbiamo
          #   già fatto l'update.
          beta_actual <- rep(0, p)
          if (j == 1){
            select <- cbind(1:p, 
                            vapply(labels_res, function(x) x[i, t, d - 1], integer(1)))
          } else {
            select_means1 <- cbind(1:(j-1),
                                   vapply(labels_res, function(x) x[i, t, d], FUN.VALUE = integer(1))[1:(j-1)])
            select_means2 <- cbind(j:p, 
                                   vapply(labels_res, function(x) x[i, t, d - 1], FUN.VALUE = integer(1))[j:p])
          }
          
          
          #@@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@#
          #@@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@#
          ### DEVO CAMPIONARE POSSIBILE NUOVO VALORE PER BETA ###
          #@@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@#
          #@@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@#
          
          
          
          #@@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@#
          #@@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@#
          ### DEVO RISCRIVERE LA FUNZIONE up_label_i ###
          #@@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@#
          #@@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@ @@@#
          labels_temp[i] <-  up_label_i(i = i, j = j,
                                        b_it =  b_res[i, , t, d],
                                        b_itm1 = if (t > 1) {b_res[i, , t-1, d]} else {rep(0, 3*p)},
                                        mustarvec_jt = mustar_temp,
                                        muvec_t = muvec_actual,
                                        T_ = param_list[[i]]$T[, , 1],
                                        Q = param_list[[i]]$Q[, , 1],
                                        rho_t = lab2rho(labels_temp),
                                        rho_tp1 = rho_tp1,
                                        gamma_tp1 = if (t == T_final) {'last time'} 
                                        else {gamma_res[[j]][, t+1, d-1]},
                                        newclustermean = possiblenewmean)
          
          # Se la label dell'unità i è maggiore di tutte le altre, vuol dire
          #   che è andato in un nuovo cluster e non in uno di quelli esistenti,
          #   quindi aggiungo la media del nuovo cluster al vettore delle medie
          #   dei cluster...
          if (all(labels_temp[i] > labels_temp[-i])){
            mustar_temp <- append(mustar_temp, possiblenewmean)
          }
          
          
          # Sistemo le labels in modo tale che non ci siano buchi (tipo non 
          # voglio c(1, 2, 2, 4) ma c(1, 2, 2, 4)) e che siano ordinate
          # (tipo non voglio c(2, 1, 2, 3) ma c(1, 2, 1, 3)).
          
          # cat("Prima:", labels_temp, mustar_temp, " ")
          labels_new <- c(1L)
          counter <- 2L
          for (jj in 2:n){
            if (!(labels_temp[jj] %in% labels_temp[1:(jj-1)])){
              labels_new[jj] <- counter
              counter <- counter+1
            } else{
              labels_new[jj] <- labels_new[which(labels_temp == labels_temp[jj])[1]]
            }
          }
          mustar_temp <- mustar_temp[unique(labels_temp)]
          labels_temp <- labels_new
          # cat("Dopo:", labels_temp, mustar_temp, "\n")  
          
        } # Fine ciclo sulle osservazioni per labels
        
        labels_res[[j]][, t, d] <- as.integer(labels_temp)
        
        
        ### ### ### ### ### ### ### #
        ### ### UPDATE MUSTAR ### ###
        
        # Devo recuperare le medie dei b_it - T*b_{i t-1}, che dipende dalle
        #   medie dei cluster e anche dalle label per ogni i.
        # Devo stare attento a che iterazione considerare per le label: se 
        #   considero la label nuova (iteraz. d) potrei avere un cluster in
        #   più dell'iteraz. precedente e quindi non avere la mustar corrispendente.
        #   Perciò devo usare la label nuova per i coeff. per cui ho già aggiornato
        #   le mustar dei cluster (perché ho già il valore anche per un eventuale
        #   nuovo cluster), mentre devo usare le label vecchie per i coeff. per 
        #   cui non ho ancora aggiornato le mustar dei cluster.
        #   Perciò devo fare questa selezione all'interno del ciclo su "j" e non
        #   posso farlo una volta sola fuori dal ciclo.
        mu_sel <- matrix(0, n, 3*p)
        for (jj in 1:p){
          if (jj < j){
            sel_m <- labels_res[[jj]][, t, d]
            mu_sel[, 1 + 3*(jj-1)] <- mustar_res[ , , t, d][jj, sel_m]
          } else if (jj == j){ 
            # Per il coefficiente j corrente uso gli oggetti temporanei creati
            #   durante l'aggiornamento delle labels. 
            #   Per le label potrei anche usare labels_res all'iteraz. d, 
            #   ma è più facile così.
            #   Per le medie invece sono obbligato a fare così perché non ho
            #   aggiornato mustar_res visto che sono in un limbo in cui
            #   mustar_temp non è nè la versione all'iteraz. d-1 nè quella
            #   all'iteraz. d.
            sel_m <- labels_temp
            mu_sel[, 1 + 3*(jj-1)] <- mustar_temp[sel_m]
          } else{
            sel_m <- labels_res[[jj]][, t, d-1]
            mu_sel[, 1 + 3*(jj-1)] <- mustar_res[ , , t, d-1][jj, sel_m]
          }
        }
        
        
        # # Recupero le partizioni al tempo t. Devo aggiornarle rispetto a quelle
        # #   trovate per aggiornare gamma e labels, perché ora ho la versione 
        # #   aggiornata delle label all'iterazione d.
        # rho_t <- lab2rho(labels_res[[j]][, t, d])
        # n_cluster <- max(labels_res[[j]][, t, d])
        
        # La partizione al tempo t aggiornata c'è in labels_temp, quindi posso
        #   farla più semplice
        rho_t <- lab2rho(labels_temp)
        n_cluster <- max(labels_temp)
        
        for (k in 1:n_cluster){ # Inizio ciclo sui cluster per mustar
          mustar_temp[k] <- up_mustar(j = j, k = k, t = t,
                                      b_t = b_res[, , t, d],
                                      b_tm1 = b_res[, , t-1, d],
                                      mu_t = mu_sel,
                                      T_list = lapply(model_list, function(x) x$fit$T[,,1]),
                                      Q_list = lapply(model_list, function(x) x$fit$Q[,,1]),
                                      rho_t = rho_t,
                                      priorm = m0,
                                      priorv = v0)
        } # Fine ciclo sui cluster per mustar
        
        mustar_res[j, 1:n_cluster, t, d] <- mustar_temp
        
      } # Fine ciclo sugli istanti "t"
    } # Fine ciclo sui coefficienti "j"
  } # Fine ciclo sulle iterazioni "d"
}
