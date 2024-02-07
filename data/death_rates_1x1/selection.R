rm(list = ls())
library(tidyverse)
library(magrittr)
setwd("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/data/death_rates_1x1")
temp = list.files(pattern="\\.txt$")
myfiles = lapply(temp, function(x) read.delim(x, skip=2, sep="",
                                              na.strings = "."))
names(myfiles) <- sapply(temp, function(x) unlist(strsplit(x, ".", fixed = T,))[1])
myfiles <- myfiles[setdiff(names(myfiles), 
                           c("DEUTE", "DEUTW", "FRACNP", "GBR_NIR", "GBR_SCO", 
                             "GBRCENW", "GBRTENW", "NZL_MA", "NZL_NM"))]
str(myfiles)

# Are there NA?
sapply(myfiles, function(x) sum(is.na(x)))
# A lot, but I want to understand in what years and for what ages

# First of all I keep only countries with data at least in the interval
#   (1945, 2020)

min <- sapply(myfiles, function(x) min(x$Year))
max <- sapply(myfiles, function(x) max(x$Year))

cbind(min, max)

sort(min[max > 2019])
sum(min < 1939 & max > 2020); which(min < 1939 & max > 2020)
sum(min < 1939 & max > 2019); which(min < 1939 & max > 2019)

tokeep <-names(myfiles)[min < 1939 & max > 2019]


data_list <- myfiles[tokeep]
names(data_list) <- sapply(tokeep, function(x) substr(x, 1, 3))


clean_years_ages <- function(data, years = 1933:2020, ages = 0:100){
  out <- data %>%
    select(Year, Age, Male, Female) %>%
    filter(Year %in% years,
           Age %in% ages) #%>%
    # pivot_wider(names_from = Age, values_from = Male) %>%
    # select(-Year) %>% mutate_if(is.character, as.numeric) %>%
    # as.matrix() %>% set_rownames(years)
  out
}

clean <- lapply(data_list, clean_years_ages)

# Are there still NA ?
sapply(clean, function(x) sum(is.na(x)))
# There are 15 in FIN and 20 in ISL

# Look into FIN
clean$FIN[is.na(clean$FIN$Male) | is.na(clean$FIN$Female), ]
# NA are all for 99/100 yrs old Males for some years between 1934 and 1965

# Look into ISL
clean$ISL[is.na(clean$ISL$Male) | is.na(clean$ISL$Female), ]
# NA for 99/100 yrs old Male for some years between 1937 and 1975 and for Female
#   (still 100 yrs old) in 1945 and 1946


clean2 <- lapply(data_list, function(x) clean_years_ages(x, ages = 0:98))

# Are there still NA ?
sapply(clean2, function(x) sum(is.na(x)))
# No

data <- clean2

rm(list = c("clean", "clean2", "data_list", "myfiles", "max", "min", "temp", "tokeep", "clean_years_ages"))
# save.image("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/death_rates_1x1/15_countries.RData")