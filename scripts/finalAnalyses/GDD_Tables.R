library(dplyr)
library(ggplot2)
library(MuMIn)
library(dplyr)
library(ggplot2)
library(sf)
library(lme4)
library(lmerTest)
library(sjPlot)
library(merTools)
library(fpp2)

# bird data
## read in phenometrics
arr <- read.csv("Outputs/birds_hopkinsCorrected.csv")
leps <- read.csv("data/phenoEstimates/adult_bfly_phenometrics_noCountCircles_withFull2020Data.csv") %>% 
  dplyr::rename(cell = HEXcell)
fledge <- readRDS("data/phenoEstimates/MAPS-fledge-dates-2022-02-22.rds") %>% 
  mutate(species = stringr::str_replace(sci_name, pattern = " ",
                                        replacement = "_"))
bird_pc <- readRDS("data/bird_PC_vals.rds")

#km2 <- bird_pc %>% 
#  dplyr::select(PC1, PC2)
#km <- kmeans(x = km2, centers = 3)
#
#bird_pc <- cbind(bird_pc, km$cluster) %>% 
#  rename(cluster = 4)
#
#ggplot(bird_pc, mapping = aes(x = PC1, y = PC2, color = as.character(cluster))) +
#  geom_point() +
#  theme_bw()
#
#ggsave(filename = "figures/clustering.png")

# read in gdd data
gdd <- data.table::fread('data/gdd_calcs.txt')
# read in frost free period data
ffp <- read.csv("Outputs/frostFreePeriod_byHex.csv")
# read in annual temp data
temp <- read.csv("Outputs/temp_byHex.csv")

# combine phenometrics with gdd & ffp
arr_gdd <- left_join(arr, gdd)
arr_gdd <- left_join(arr_gdd, ffp)
arr_gdd <- left_join(arr_gdd, temp)
fledge_gdd <- left_join(fledge,gdd)
fledge_gdd <- left_join(fledge_gdd, ffp) %>% 
  filter(!is.na(mean_ffp))
leps_gdd <- left_join(leps, gdd)
leps_gdd <- left_join(leps_gdd, ffp)
leps_gdd <- left_join(leps_gdd, temp)

# what years do we cover over this analysis?
# arrival is 2002-2017
# fledge is 1992 - 2018
# leps is 2000-2020
# greenup is currently from arrival so 
#2002-2017 should be the dates of the current analyses


# add PC based clusters to df
#arr_gdd_pc <- left_join(arr_gdd, bird_pc)
#fledge_gdd_pc <- left_join(fledge_gdd, bird_pc, by = "species")

# look into sensitivty of phenometrics to gdd through space and time
# start with greenup
greenup <- arr_gdd %>% 
  dplyr::distinct(cell, year, .keep_all = T) %>% 
  dplyr::select(cell, year, gr_mn, spring.gdd, summer.gdd, 
                spring.dev, summer.dev, mean_ffp) %>% 
  dplyr::filter(year >= 2002 & year <= 2017)
greenup$cell <- as.character(greenup$cell)
# scale greenup dataframe
greenup_scaled <- greenup %>% 
  mutate(spring.dev = scale(spring.dev),
         mean_ffp = scale(mean_ffp))

# deviation greenup model
gu_m <-  lmer(formula = gr_mn ~ spring.dev + mean_ffp +
                spring.dev:mean_ffp +
                (1|cell), data = greenup_scaled,
              na.action = na.fail, REML = F)

gu_d <- dredge(gu_m)

gu_tm <- get.models(gu_d, subset = 1)[[1]]
gu_tm

gu_tm_unscaled <-  lmer(formula = gr_mn ~ spring.dev + mean_ffp +
                spring.dev:mean_ffp +
                (1|cell), data = greenup,
              na.action = na.fail, REML = F)

## now get this for fledge and leps
fledge_gdd <- fledge_gdd %>% 
  filter(!is.na(juv_meanday),
         !is.na(spring.gdd),
         !is.na(spring.dev)) %>% 
  dplyr::select(sci_name, station, year, cell, juv_meanday, 
                spring.gdd, summer.gdd,
                spring.dev, summer.dev, PC1, PC2, species, mean_ffp) %>% 
  filter(year >= 2002 & year <= 2017)

head(fledge_gdd)

fledge_scaled <- ungroup(fledge_gdd) %>% 
  mutate(spring.dev = scale(spring.dev),
         mean_ffp= scale(mean_ffp)) %>% 
  filter(!is.na(juv_meanday),
         !is.na(spring.gdd),
         !is.na(spring.dev),
         !is.na(year)) 

# deviation fledge model
fledge_m <- lmer(juv_meanday ~ spring.dev + mean_ffp +
                   PC1 +
                   spring.dev:mean_ffp +
                   spring.dev:PC1 +
                   (1|station) + (1|sci_name),
                 data = fledge_scaled, na.action = na.fail, REML = F)

#dredge
fledge_d <- dredge(fledge_m)
tm_fledge <- get.models(fledge_d, subset = 1)[[1]]

## now lep time
# deviation model 
leps_gdd <- leps_gdd %>% 
  filter(!is.na(spring.dev),
         !is.na(year),
         !is.na(mean_ffp),
         !is.na(year),
         !is.na(cell),
         !is.na(code),
         !is.na(q5)) %>% 
  filter(year >= 2002 & year <= 2017)

leps_scaled <- leps_gdd %>% 
  mutate(spring.dev = scale(spring.dev),
         mean_ffp = scale(mean_ffp),
         uniqObsDays = scale(uniqObsDays)) 

# note that we can't look at the interaction between code and year because there 
# is not enough RE data in early years
leps_m <- lmer(q5 ~ spring.dev + mean_ffp +
                 uniqObsDays + code +
                 spring.dev:code +
                 spring.dev:mean_ffp +
                 (1|cell), 
               data = leps_scaled, na.action = na.fail, REML = F)

leps_d <- dredge(leps_m)
leps_tm <- get.models(leps_d, 2)[[1]]

## three top models
summary(gu_tm)
MuMIn::r.squaredGLMM(gu_tm)
summary(tm_fledge)
MuMIn::r.squaredGLMM(tm_fledge)
summary(leps_tm)
MuMIn::r.squaredGLMM(leps_tm)
