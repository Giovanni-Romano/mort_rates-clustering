rm(list=ls())

source("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/BSP_Pavone/source/BSP.R")
source("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/BSP_Pavone/source/setup.R")
load("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/BSP_Pavone/output/mortality.Rdata")

require(tidyverse)
require(KFAS)
require(doParallel)

set.seed(5477)

data_list_man <- list(ita_man = Y_ita_man / N_ita_man,
                      swe_man = Y_swe_man / N_swe_man,
                      uk_man = Y_uk_man / N_uk_man,
                      us_man = Y_us_man / N_us_man)

data_list_woman <- list(ita_woman = Y_ita_woman / N_ita_woman,
                        swe_woman = Y_swe_woman / N_swe_woman,
                        uk_woman = Y_uk_woman / N_uk_woman,
                        us_woman = Y_us_woman / N_us_woman)

# Model definition
model_list_man <- lapply(data_list_man, 
                         . %>% bsp.model(rates = .,
                                         delta = delta,
                                         age_knots = age_knots,
                                         kernel = independent_kernel))

model_list_woman <- lapply(data_list_woman, 
                           . %>% bsp.model(rates = .,
                                           delta = delta,
                                           age_knots = age_knots,
                                           kernel = independent_kernel))


# Model fit
fit_list_man <- lapply(model_list_man,
                       . %>%
                         bsp.fit(., 
                                 rep = 10, 
                                 method = 'Nelder-Mead', 
                                 parallel = TRUE, 
                                 maxcl = 10))

fit_list_woman <- lapply(model_list_woman,
                         . %>%
                           bsp.fit(., 
                                   rep = 10, 
                                   method = 'Nelder-Mead', 
                                   parallel = TRUE, 
                                   maxcl = 10))

  
save(list = c('fit_list_man',
              'fit_list_woman'),
     file = 'fit_indep.RData')

