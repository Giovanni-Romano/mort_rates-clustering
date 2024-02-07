rm(list=ls())
require(tidyverse)
require(magrittr)

# Data loading
load("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/data/death_rates_1x1/14_countries.RData")

# Years of analysis
years <- 1933:2019
ages <- 0:98

# Select Male/Female
select_male <- function(data, years = 1933:2020, ages = 0:98){
  out <- data %>%
    select(Year, Age, Male) %>%
    filter(Year %in% years,
           Age %in% ages) %>%
    pivot_wider(names_from = Age, values_from = Male) %>%
    select(-Year) %>% mutate_if(is.character, as.numeric) %>%
    as.matrix() %>% set_rownames(years)
  out
}

select_female <- function(data, years = 1933:2020, ages = 0:98){
  out <- data %>%
    select(Year, Age, Female) %>%
    filter(Year %in% years,
           Age %in% ages) %>%
    pivot_wider(names_from = Age, values_from = Female) %>%
    select(-Year) %>% mutate_if(is.character, as.numeric) %>%
    as.matrix() %>% set_rownames(years)
  out
}

data_male <- lapply(data, select_male)
data_female <- lapply(data, select_female)

# Data checks
# NA (there aren't)
sapply(data_male, function(x) sum(is.na(x)))
sapply(data_female, function(x) sum(is.na(x)))
# Negative values (there aren't)
sapply(data_male, function(x) sum(x < 0))
sapply(data_female, function(x) sum(x < 0))
# Exact 0 (there are some)
sapply(data_male, function(x) sum(x == 0))
sapply(data_female, function(x) sum(x == 0))


# Remove ISL because there are too many 0 to be modeled with our model
data_noISL <- data[setdiff(names(data), "ISL")]
data_male <- lapply(data_noISL, select_male)
data_female <- lapply(data_noISL, select_female)

# For BEL M, CHE M&F, DNK M&F, NOR M&F, SWE M&F there are some 0 but not 
#   too many. Following Pavone's work, I estimate these values as the mean of 
#   the two closest ages and years average

do_avg <- function(data){
  
  out <- data
  
  idx_0 <- which(data == 0, arr.ind = T)
  
  if (nrow(idx_0) > 0){
    
    M1 <- matrix(c(+1,0, -1,0, 0,+1, 0,-1), nrow = 2)
    M2 <- cbind(M1, 
                matrix(c(+2,0, -2,0, 0,+2, 0,-2, 
                         -1,-1, +1,+1, -1,+1, +1, 1), nrow = 2))
    
    for(it in 1:nrow(idx_0)){
      
      rc <- idx_0[it, ]
      
      tmp <- t(rc + M1)
      sel <- tmp[tmp[ , 1]>0 & tmp[ , 1]<=nrow(data) & tmp[ , 2]>0 & tmp[ , 2]<=ncol(data), ]
      
      data_sel <- data[sel]
      
      if (any(data_sel > 0)){
        out[rc[1], rc[2]] <- mean(data_sel)
      } else {
        # If all first neighbours are 0, then I go to second neighbours
        tmp <- t(rc + M2)
        sel <- tmp[tmp[ , 1]>0 & tmp[ , 1]<=nrow(data) & tmp[ , 2]>0 & tmp[ , 2]<=ncol(data), ]

        data_sel <- data[sel]
        
        out[rc[1], rc[2]] <- mean(data_sel)
      }
      # I hope second neighbours are enough
    }
  }
    
  out
}

data_male_no0 <- lapply(data_male, do_avg)
data_female_no0 <- lapply(data_female, do_avg)
# Exact 0
sapply(data_male_no0, function(x) sum(x == 0))
sapply(data_female_no0, function(x) sum(x == 0))
# There aren't 0 anymore (--> no need to go further than 2nd neighb.)
# NA
sapply(data_male_no0, function(x) sum(is.na(x)))
sapply(data_female_no0, function(x) sum(is.na(x)))



# Objects to save
rates_male <- data_male_no0
rates_female <- data_female_no0

setwd("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/death_rates_1x1")

save(rates_male,
     rates_female,
     years,
     ages,
     file = "14_countries_preproc.RData")

