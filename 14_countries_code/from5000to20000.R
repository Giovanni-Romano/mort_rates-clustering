rm(list = ls())
load("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/res/dynamic_on_beta/no_hierarchy/14_countries/full_bayes/res_prior_on_M_5000.RData")
n_iter <- 20000
gamma_res_old <- gamma_res
gamma_res <- replicate(p,
                       array(NA,
                             dim = c(n, # numero di osservazioni
                                     T_final, # numero istanti temporali
                                     n_iter # numero iterazioni
                             )
                       ),
                       simplify = FALSE)
labels_res_old <- labels_res
labels_res <- replicate(p,
                        array(NA,
                              dim = c(n, # numero di osservazioni
                                      T_final, # numero istanti temporali
                                      n_iter # numero iterazioni
                              )
                        ),
                        simplify = FALSE)
beta_res_old <- beta_res
beta_res <- replicate(p,
                      array(NA,
                            dim = c(n, # numero di osservazioni
                                    T_final, # numero istanti temporali
                                    n_iter # numero iterazioni
                            )
                      ),
                      simplify = FALSE)
delta_res_old <- delta_res
delta_res <- replicate(p,
                       array(NA,
                             dim = c(n_iter) # numero iterazioni
                       ),
                       simplify = FALSE)
alpha_res_old <- alpha_res
alpha_res <- replicate(p,
                       array(NA,
                             dim = c(n_iter # numero iterazioni
                             )
                       ),
                       simplify = FALSE)
sigma_res_old <- sigma_res
sigma_res <- array(NA,
                   dim = c(n,
                           n_iter # numero iterazioni
                   )
)
M_res_old <- M_res
M_res <- replicate(p, 
                   array(NA, dim = c(n_iter)),
                   simplify = F)
for (j in 1:p){
  
  
  gamma_res[[j]][,,1:5000] <- gamma_res_old[[j]]
  
  labels_res[[j]][,,1:5000] <- labels_res_old[[j]]
  
  
  beta_res[[j]][,,1:5000] <- beta_res_old[[j]]
  
  
  delta_res[[j]][1:5000] <- delta_res_old[[j]]
  alpha_res[[j]][1:5000] <- alpha_res_old[[j]]
  M_res[[j]][1:5000] <- M_res_old[[j]]
}

sigma_res[,1:5000] <- sigma_res_old
set.seed(23)






