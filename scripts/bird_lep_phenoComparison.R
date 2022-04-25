library(dplyr)
library(ggplot2)

## read in phenometrics
birds <- read.csv("Outputs/birds_hopkinsCorrected.csv")
leps <- read.csv("data/phenoEstimates/adult_bfly_phenometrics_noCountCircles_withFull2020Data.csv")

head(birds)
## get benchmark values for group, cell, year
## of years 2016-2020 for leps and 
## of years 2013 - 2017 for birds
## first birds
bird_bm <- birds %>% 
  dplyr::filter(year >= 2013 & year <= 2017) %>% 
  group_by(cell, season) %>% 
  summarise(bmPheno = mean(arr_GAM_mean),
            bmPheno_SD = sd(arr_GAM_mean),
            uniqueYears = length(unique(year)),
            uniqueSpp = length(unique(species)))

## now leps
lep_bm <- leps %>% 
  dplyr::filter(q50 >= 152 & q50 <=243) %>% 
  dplyr::filter(year >= 2016 & year <= 2020) %>% 
  group_by(HEXcell, code) %>% 
  summarise(bmPheno_5 = mean(q5),
            bmPheno_50 = mean(q50),
            bmPheno_5_SD = sd(q5),
            bmPheno_50_SD = sd(q50),
            uniqueYears = length(unique(year)),
            uniqueSpp = length(unique(code)))


## combine these benchmark values with original values
leps2 <- left_join(leps, lep_bm) %>% 
  dplyr::filter(uniqueYears >= 5) %>% 
  dplyr::mutate(dev_5_lep = q5 - bmPheno_5,
                dev_50_lep = q50 - bmPheno_50,
                sd_5_lep = (q5 - bmPheno_5) / bmPheno_5_SD,
                sd_50_lep = (q50 - bmPheno_50) / bmPheno_50_SD) %>% 
  dplyr::select(HEXcell, year, code, dev_5_lep, dev_50_lep, sd_5_lep, sd_50_lep, q5) %>% 
  dplyr::rename(cell = HEXcell)

## for birds, we must first combine phenometrics
birds2 <- birds %>% 
  group_by(cell, year, season) %>% 
  summarise(mean_arr = mean(arr_GAM_mean),
            sd_arr = sd(arr_GAM_mean))

birds3 <- left_join(birds2, bird_bm) %>% 
  dplyr::filter(uniqueYears >= 5) %>% 
  dplyr::mutate(dev_arr_bird = mean_arr - bmPheno,
                sd_arr_bird = (mean_arr - bmPheno) / sd_arr) %>% 
  dplyr::select(cell, year, season, dev_arr_bird, sd_arr_bird, mean_arr)

## combine bird and lepcells
lep_bird <- left_join(leps2, birds3) %>% 
  dplyr::filter(!is.na(season)) 


ggplot(data = lep_bird) +
  geom_point(mapping = aes(x = q5, y = mean_arr)) +
  geom_smooth(mapping = aes(x = q5, y = mean_arr), method = "lm") +
  facet_wrap(code ~ season) +
  labs(x = "Adult butterfly emergence", y = "Mean Bird Arrival") +
  theme_bw()


# look into q50 results
lep_bird <- left_join(leps2, birds3) %>% 
  dplyr::filter(!is.na(season)) %>% 
  dplyr::filter(dev_50_lep < 50,
                dev_50_lep > -50)


ggplot(data = lep_bird) +
  geom_point(mapping = aes(x = dev_50_lep, y = dev_arr_bird)) +
  geom_smooth(mapping = aes(x = dev_50_lep, y = dev_arr_bird), method = "lm") +
  facet_wrap(code ~ season) +
  theme_light()

## make these 9 linear models, plot the residuals spatially and temporally
# egg, early
ee <- dplyr::filter(lep_bird, code == "RE" & season == "early") %>% 
  dplyr::filter(!is.na(dev_5_lep))
ee_m <- lm(dev_arr_bird ~ dev_5_lep, data = ee)
r <- residuals(ee_m)
ee_r <- ee %>% 
  mutate(resid = r)

ggplot(ee_r, mapping = aes(x = year, y = resid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_light()
