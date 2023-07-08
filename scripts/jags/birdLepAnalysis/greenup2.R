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
# add cell_lat to dataframe
lat_df <- dplyr::distinct(arr, cell, .keep_all = T) %>% 
  select(cell, cell_lat)
greenup <- left_join(greenup, lat_df, by = "cell")
# munge data for future model building
cellSort <- unique(greenup$cell) %>% sort()
cellRanking <- data.frame(cell = cellSort, cellID = 1:length(cellSort))
cellRanking$cellID <- as.character(cellRanking$cellID)
greenup$cell <- factor(greenup$cell)
greenup$cell <- droplevels(greenup$cell)
greenup$cell <- as.integer(greenup$cell)


# trying to build a model looking at response of greenup in a cell to gdd 
# random slopes and intercepts

params <- c("alpha", "beta")

DATA <- list(y = as.vector(greenup$gr_mn),
             x = as.vector(greenup$spring.gdd.center),
             group = as.vector(greenup$cell),
             lat = as.vector(greenup$cell_lat),
             N = length(greenup$gr_mn),
             M = length(unique(greenup$cell)))

fit <- jagsUI::jags(model.file = 'scripts/jags/birdLepAnalysis/greenup_latEffect.jags',
                    data = DATA,
                    parameters.to.save = c("alpha", "beta", "theta"),
                    n.chains = 4,
                    n.burnin = 10000, #number of 'burn-in' iterations - so the model can start to find the best areas of parameter space before tracking MCMC values
                    n.iter = 20000, #number of samples to draw from the posterior + burnin
                    parallel = TRUE) #whether to run each chain on a separate core (you should) 

summary(fit)
# get output summary
#Summary
MCMCvis::MCMCsummary(fit, 
                     params = c('alpha', 'beta','theta'),
                     round = 3)
#Look at trace plots
MCMCvis::MCMCtrace(fit, 
                   ind = TRUE, 
                   pdf = FALSE)
#Caterpillar plots
MCMCvis::MCMCplot(fit, 
                  params = c('alpha', 'beta'))
#Extract posterior values for a given parameter
MCMCvis::MCMCchains(fit, params = 'beta')
#Calculate posterior mean for a given parameter
MCMCvis::MCMCpstr(fit, params = 'alpha', fun = mean)[[1]]
MCMCvis::MCMCpstr(fit, params = 'beta', fun = mean)[[1]]
hist(MCMCvis::MCMCpstr(fit, params = 'beta', fun = mean)[[1]])
#Extract posterior values for a given parameter
hist(MCMCvis::MCMCchains(fit, params = 'theta'))
MCMCvis::MCMCpstr(fit, params = 'theta', fun = mean)

# visualize results
cell_df <- distinct(greenup, cell, .keep_all = T) %>% 
  select(cell, cell_lat)
beta_df <- data.frame(beta = MCMCvis::MCMCpstr(fit, params = 'beta', fun = mean)[[1]],
                      cell = 1:69)
alpha_df <- data.frame(alpha = MCMCvis::MCMCpstr(fit, params = 'alpha', fun = mean)[[1]],
                       cell = 1:69)

tdf <- cell_df %>% 
  left_join(beta_df) %>% 
  left_join(alpha_df)

ggplot(tdf, aes(x = cell_lat, y = alpha)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(tdf, aes(x = cell_lat, y = beta)) +
  geom_point() +
  geom_smooth(method = "lm")
