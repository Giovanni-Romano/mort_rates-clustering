rm(list = ls())
library(tidyverse)

countries_model <- c("AUS", "BEL", "CAN", "DNK", "FIN", "FRA", "ITA", "NLD", 
                     "NOR", "ESP", "SWE", "CHE", "GBR", "USA")

#### GDP ####
gdp <- read.csv("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/covariate_data/OECD/gdp.csv")
str(gdp)

# Fix value for Flag.Codes that are empty "" (i.e. it's true data)
gdp$Flag.Codes[gdp$Flag.Codes == ""] <- "T"

# Filter only countries in the model, pick USD dollars per capita as measure and
# filter for Year <= 2020
gdp %>% filter(LOCATION %in% countries_model,
               TIME <= 2020,
               FREQUENCY == "A",
               MEASURE == "USD_CAP") -> gdp1

# Are there any missing values?
gdp1 %>% 
  group_by(LOCATION) %>% 
  summarise(sum_na = sum(is.na(Value)))

# How many years of data do we have?
table(gdp1$TIME)
# We have all 14 countries for 1970-2020

# How many flags?
# "E" = Estimated data
# "P" = Provisional data
gdp1 %>% 
  filter(TIME >= 1970) %>% 
  group_by(TIME, Flag.Codes) %>%
  summarise(n = n()) %>%
  spread(Flag.Codes, n)%>%
  print(n = Inf)
# We have "E" for 4+ countries every year before 1995, is it ok?

gdp_clean <- gdp1 %>% filter(TIME >= 1970)



#### Health_spending ####
# I checked data availability also w/ "PC_GDP" as MEASURE, but it's the same
health_spend <- read.csv("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/covariate_data/OECD/health_spending.csv")
str(health_spend)

# Fix value for Flag.Codes that are empty "" (i.e. it's true data)
table(health_spend$Flag.Codes)
health_spend$Flag.Codes[health_spend$Flag.Codes == ""] <- "T"
table(health_spend$Flag.Codes)


table(health_spend$MEASURE, health_spend$SUBJECT)
# Filter only countries in the model, pick TOT as subject and USD dollars per 
# capita as measure and filter for Year <= 2020.
# Also PC_GDP could be interesting
health_spend %>% 
  filter(LOCATION %in% countries_model,
                        TIME <= 2020,
                        FREQUENCY == "A",
                        SUBJECT == "TOT",
                        MEASURE %in% c("USD_CAP", "PC_GDP")) %>% 
  group_by(MEASURE) %>% 
  group_split() -> tmp
tmp <- merge(tmp[[1]], tmp[[2]], by = c("LOCATION", "TIME"))
cor(tmp$Value.x, tmp$Value.y) # They are highly correlated, just keep "USD_CAP"
plot(tmp$Value.x, tmp$Value.y)


health_spend %>% filter(LOCATION %in% countries_model,
                        TIME <= 2020,
                        FREQUENCY == "A",
                        SUBJECT == "TOT",
                        MEASURE == "USD_CAP") -> health_spend1

# Are there any missing values?
health_spend1 %>% 
  group_by(LOCATION) %>% 
  summarise(sum_na = sum(is.na(Value)))
#No

# How many years of data do we have?
table(health_spend1$TIME)
# We have all 14 countries for 1990-2020
# We have 12+ countries since 1972 and 11 and 10 respectively in 1971 and 1970

# How many flags?
health_spend1 %>% 
  group_by(TIME, Flag.Codes) %>%
  summarise(n = n()) %>%
  spread(Flag.Codes, n)%>%
  print(n = Inf)
# 1/2 "D" = "Difference methodology" in almost all years
# 1/2/3 "B" = "Break in time series" in some years (even if I'm not sure about its meaning)

table(health_spend1$LOCATION[health_spend1$Flag.Codes == "B"])
table(health_spend1$LOCATION[health_spend1$Flag.Codes == "D"])
# "B" is in many different countries,
# "D" is 16 times in AUS, 25 in BEL, 9 in DNK and 1 in CAN


health_spend_clean <- health_spend1 %>% filter(TIME >= 1990)



#### Income inequality ####
inc_ineq <- read.csv("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/covariate_data/OECD/income_inequality.csv")
str(inc_ineq)

# Fix value for Flag.Codes that are empty "" (i.e. it's true data)
table(inc_ineq$Flag.Codes)
inc_ineq$Flag.Codes[inc_ineq$Flag.Codes == ""] <- "T"

# Filter only countries in the model, pick GINI as SUBJECT and Year <= 2020
table(inc_ineq$INDICATOR)
table(inc_ineq$SUBJECT, inc_ineq$MEASURE)
inc_ineq %>% 
  filter(LOCATION %in% countries_model,
         TIME <= 2020,
         FREQUENCY == "A",
         SUBJECT == "GINI") -> inc_ineq_clean

# Are there any missing values?
inc_ineq_clean %>% 
  group_by(LOCATION) %>% 
  summarise(sum_na = sum(is.na(Value)))

# How many years of data do we have?
table(inc_ineq_clean$TIME)
# We have few data: never have all 14 countries, 10+ countries only since 2012 and 5+ countries only since 2006

# How many flags?
inc_ineq_clean %>% 
  group_by(TIME, Flag.Codes) %>%
  summarise(n = n()) %>%
  spread(Flag.Codes, n)%>%
  print(n = Inf)
# Only true data!




#### Unemployment rate ####
unemp <- read.csv("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/covariate_data/OECD/unemployment_rate.csv")
str(unemp)

# Fix value for Flag.Codes that are empty "" (i.e. it's true data)
table(unemp$Flag.Codes)
unemp$Flag.Codes[unemp$Flag.Codes == ""] <- "T"

# Filter only countries in the model and Year <= 2020 and SUBJECT = "TOT"
table(unemp$INDICATOR)
table(unemp$SUBJECT, unemp$MEASURE)
# "PC_LF" is the % of Labour Force
unemp %>% filter(LOCATION %in% countries_model,
                 FREQUENCY == "A",
                 TIME <= 2020,
                 SUBJECT == "TOT") -> unemp1

# Are there any missing values?
unemp1 %>% 
  group_by(LOCATION) %>% 
  summarise(sum_na = sum(is.na(Value)))

# How many years of data do we have?
table(unemp1$TIME)
# 10+ countries only since 1983, before that there are <=3.
# 13 countries between 1989 and 2009, then all 14
table(unemp1$LOCATION)
# CHE is the last one to kick in, in 2010

# How many flags?
unemp1 %>% 
  group_by(TIME, Flag.Codes) %>%
  summarise(n = n()) %>%
  spread(Flag.Codes, n)%>%
  print(n = Inf)
table(unemp_clean$LOCATION[unemp_clean$Flag.Codes == "P"])
table(unemp_clean$LOCATION[unemp_clean$Flag.Codes == "B"])
# ITA has provisional data between 1983 and 1991
# 14 "B" are divided into different countries and years

unemp_clean <- unemp1 %>% filter(TIME >= 2010)





#### Output ####
covariates <- list(GDP = gdp_clean %>% 
                     select(LOCATION, TIME, Value) %>% 
                     mutate(LOCATION = ifelse(LOCATION == "GBR", "UK", LOCATION)),
                   HEALTH_EXP = health_spend_clean %>% 
                     select(LOCATION, TIME, Value) %>% 
                     mutate(LOCATION = ifelse(LOCATION == "GBR", "UK", LOCATION)),
                   UNEMP = unemp_clean %>% 
                     select(LOCATION, TIME, Value) %>% 
                     mutate(LOCATION = ifelse(LOCATION == "GBR", "UK", LOCATION)))

saveRDS(covariates,
        "C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/covariate_data/covariates_clean.rds")


# #### Put variables together ####
# covariates <- data.frame(TIME = rep(1970:2020, 14),
#                          LOCATION = rep(countries_model, each = 51))
# 
# gdp_clean %>% 
#   select(LOCATION, TIME, Value) %>% 
#   rename(GDP = Value) %>%
#   merge(covariates, ., by = c("TIME", "LOCATION"), all.x = TRUE) %>%
#   merge(., 
#         health_spend_clean %>% select(LOCATION, TIME, Value) %>% rename(HEALTH_EXP = Value), 
#         by = c("TIME", "LOCATION"), all.x = TRUE) %>%
#   merge(.,
#         unemp_clean %>% select(LOCATION, TIME, Value) %>% rename(UNEMP = Value),
#         by = c("TIME", "LOCATION"), all.x = TRUE) -> covariates
# 
# str(covariates)
# 
# 
# # NA check
# covariates %>% 
#   filter(TIME >= 2010) %>%
#   group_by(LOCATION) %>% 
#   summarise(GDP_na = sum(is.na(GDP)),
#             HEA_na = sum(is.na(HEALTH_EXP)),
#             UNE_na = sum(is.na(UNEMP)))
# 
# # Covariates post 2010
# covariates_post2010 <- covariates %>% 
#   filter(TIME >= 2010)
# saveRDS(covariates_post2010, 
#         "C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/covariate_data/covariates_post2010.rds")
