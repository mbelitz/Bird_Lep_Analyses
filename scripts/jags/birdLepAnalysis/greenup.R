library(tidyverse)
library(R2jags)

# read in phenometrics, gdd and annual temp data
load("phenoData.RData")

# combine phenometrics with gdd & ffp
arr_gdd <- left_join(arr, gdd)
arr_gdd <- left_join(arr_gdd, ffp)
arr_gdd <- left_join(arr_gdd, temp)

# select limited number of data columns
greenup <- arr_gdd %>% 
  select(gr_mn, cell, year, spring.gdd, mean_ffp, mean_temp) 
# center gdd for each cell
greenup <- greenup %>% 
  group_by(cell) %>% 
  mutate(spring.gdd.center = scale(spring.gdd, scale = F)) %>% # center gdd by cell
  filter(year >= 2002 & year <= 2017) # restrict to same years
# remove rows without mean greenup estimate
greenup <- filter(greenup, !is.na(gr_mn))
# only include cells with at least three years of sampling
yearsGrouping <- greenup %>% 
  group_by(cell) %>% 
  summarise(nYears = length(unique(year)))
cellsToKeep <- filter(yearsGrouping, nYears >= 3)
greenup <- greenup %>% 
  filter(cell %in% cellsToKeep$cell)
# munge data for future model building
cellSort <- unique(greenup$cell) %>% sort()
cellRanking <- data.frame(cell = cellSort, cellID = 1:length(cellSort))
cellRanking$cellID <- as.character(cellRanking$cellID)
greenup$cell <- factor(greenup$cell)
greenup$cell <- droplevels(greenup$cell)
greenup$cell <- as.integer(greenup$cell)

# trying to build a model looking at response of greenup in a cell to gdd 
# random slopes and intercepts

lm_jags_re_slope <- function(){
  
  # likelihood:
  for(i in 1:N){
    y[i] ~ dnorm(mu[i], tau) # tau is precision (1/variance)
    mu[i] <- alpha + a[cell[i]] + beta + b[cell[i]] * x[i] # equation
  }
  
  #priors:
  alpha ~ dnorm(0, 0.01) # intercept
  sigma_a ~ dunif(0,100) # standard deviation of random effect
  tau_a <- 1 / (sigma_a * sigma_a) #converted to precision
  for(j in 1:Ncells){
    a[j] ~ dnorm(0, tau_a) # random intercept for each cell
  }
  beta ~ dnorm(0, 0.01) # slope
  sigma_b ~ dunif(0,100) # s.d. of random slope
  tau_b <- 1 / (sigma_b * sigma_b)
  for(j in 1:Ncells){
    b[j] ~ dnorm(0, tau_b)
  }
  # set priors for precision in terms of sd
  sigma ~ dunif(0,100) # standard deviation of fixed effect
  tau <- 1 / (sigma * sigma) # convert to precision
}

init_values <- function(){
  list(alpha = rnorm(1),
       sigma_a = runif(1),
       sigma_b = runif(1),
       sigma = runif(1),
       beta = rnorm(1))
}

params <- c("alpha", "beta", "sigma", "sigma_a", "a", "b")

DATA <- list(y = as.vector(greenup$gr_mn),
             x = as.vector(greenup$spring.gdd.center),
             cell = as.vector(greenup$cell),
             N = length(greenup$gr_mn),
             Ncells = length(unique(greenup$cell)))

fit_lm_re_slope <- jags(data = DATA, inits = init_values, parameters.to.save = params, 
               model.file = lm_jags_re_slope,
               n.chains = 3, n.iter = 20000, n.burnin = 5000, n.thin = 10, DIC = F)

# get output summary
gu_est_table <- fit_lm_re_slope$BUGSoutput$summary %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "coef")

gu_est_table <- gu_est_table %>% 
  mutate(cellID = str_extract(coef, "(?<=\\[)[^]]+"),
         ranef = str_extract(coef,'[^\\[]+'))

# assing cell ID to effects
gu_est_table <- left_join(gu_est_table, cellRanking)

lat_df <- dplyr::distinct(arr, cell, .keep_all = T) %>% 
  select(cell, cell_lat)

gu_est_table <- left_join(gu_est_table, lat_df, by = "cell") 

# plot slope for each cell
gu_est_table %>% 
  filter(ranef == "b") %>% 
  ggplot() +
  geom_point(aes(x = cell_lat, y = mean))

# those cells that deviate from the latitude trend all have few number of data points,
# likely providing strong evidence we should increase the minimum year requirement in a cell