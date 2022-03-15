library(dplyr)
library(ggplot2)
library(sf)
library(lme4)
library(lmerTest)
library(sjPlot)

# bird data
## read in phenometrics
arr <- read.csv("Outputs/birds_hopkinsCorrected.csv")
leps <- read.csv("data/phenoEstimates/adult_bfly_phenometrics_noCountCircles_withFull2020Data.csv") %>% 
  dplyr::rename(cell = HEXcell)
fledge <- readRDS("data/phenoEstimates/MAPS-fledge-dates-2022-02-22.rds") %>% 
  mutate(species = stringr::str_replace(sci_name, pattern = " ",
                                        replacement = "_"))
bird_pc <- readRDS("data/bird_PC_vals.rds")

km2 <- bird_pc %>% 
  dplyr::select(PC1, PC2)
km <- kmeans(x = km2, centers = 3)

bird_pc <- cbind(bird_pc, km$cluster) %>% 
  rename(cluster = 4)

ggplot(bird_pc, mapping = aes(x = PC1, y = PC2, color = as.character(cluster))) +
  geom_point() +
  theme_bw()

# read in gdd data
gdd <- data.table::fread('data/gdd_calcs.txt')

# combine phenometrics with gdd
arr_gdd <- left_join(arr, gdd)
fledge_gdd <- left_join(fledge,gdd)
leps_gdd <- left_join(leps, gdd)

# add PC based clusters to df
arr_gdd_pc <- left_join(arr_gdd, bird_pc)
fledge_gdd_pc <- left_join(fledge_gdd, bird_pc, by = "species")

# look into sensitivty of phenometrics to gdd through space and time
# start with greenup
greenup <- arr_gdd_pc %>% 
  dplyr::distinct(cell, year, .keep_all = T) %>% 
  dplyr::select(cell, year, gr_mn, spring.gdd, summer.gdd, spring.dev, summer.dev)
greenup$cell <- as.character(greenup$cell)

greenup_scaled <- greenup %>% 
  mutate(year = scale(year), 
         spring.gdd = scale(spring.gdd),
         spring.dev = scale(spring.dev))


# spatiotemporal greenup model
sp_gu <- lmer(formula = gr_mn ~  year + spring.gdd:year +
                (1|cell), data = greenup)

summary(sp_gu)

# deviation greenup model
t_gu <-  lmer(formula = gr_mn ~ spring.dev + year +
                  (1|cell), data = greenup)

summary(t_gu)

# let's do birds next
arr_gdd_pc$cluster <- as.character(arr_gdd_pc$cluster)
arr_gdd_pc <- arr_gdd_pc %>% 
  filter(!is.na(arr_GAM_mean)) %>% 
  dplyr::select(species, cell, year, arr_GAM_mean, 
                arr_IAR_mean, arr_IAR_sd, spring.gdd, summer.gdd,
                spring.dev, summer.dev, cluster)

arr_gdd_pc_group <- arr_gdd_pc %>% 
  group_by(cell, year, cluster) %>% 
  summarise(mean_arr = mean(arr_IAR_mean))

arr_gdd_pc_group <- left_join(arr_gdd_pc_group, gdd)

arr_gdd_scaled <- arr_gdd_pc_group %>% 
  mutate(year2 = t, 
         spring.gdd2 = t2,
         spring.dev2 = t3)

# spatiotemporal arrival model
sp_ar <- lmer(mean_arr ~ spring.gdd + year + spring.gdd:year +
                cluster + 
                spring.gdd:cluster +
                year:cluster +
                spring.gdd:year:cluster +
                (1|cell), data = arr_gdd_pc_group)

summary(sp_ar)

plot_model(sp_ar, type = "pred", terms = c("year", "spring.gdd", "cluster"))

# deviation model 
t_ar <- lmer(mean_arr ~ spring.dev + year + 
                cluster + 
                spring.dev:cluster +
                year:cluster +
                (1|cell), data = arr_gdd_pc_group)

summary(t_ar)

plot_model(t_ar, type = "pred", terms = c("year", "cluster"))
plot_model(t_ar, type = "pred", terms = c("spring.dev", "cluster"))


## now let's do bird fledging
fledge_gdd_pc$cluster <- as.character(fledge_gdd_pc$cluster)
fledge_gdd_pc <- fledge_gdd_pc %>% 
  filter(!is.na(juv_meanday),
         !is.na(spring.gdd),
         !is.na(spring.dev)) %>% 
  dplyr::select(species, cell, year, juv_meanday, 
                spring.gdd, summer.gdd,
                spring.dev, summer.dev, cluster)

fledge_gdd_pc_group <- fledge_gdd_pc %>% 
  group_by(cell, year, cluster) %>% 
  summarise(mean_fledge = mean(juv_meanday))

fledge_gdd_pc_group <- left_join(fledge_gdd_pc_group, gdd)

# spatiotemporal arrival model
sp_fledge <- lmer(mean_fledge ~ spring.gdd + year + spring.gdd:year +
                cluster + 
                spring.gdd:cluster +
                year:cluster +
                spring.gdd:year:cluster +
                (1|cell), data = fledge_gdd_pc_group)

summary(sp_fledge)

plot_model(sp_fledge, type = "pred", terms = c("year", "spring.gdd", "cluster"))

# deviation model 
t_fledge <- lmer(mean_fledge ~ spring.dev + year + 
               cluster + 
               spring.dev:cluster +
               year:cluster +
               (1|cell), data = fledge_gdd_pc_group)

summary(t_fledge)

plot_model(t_fledge, type = "pred", terms = c("year", "cluster"))
plot_model(t_fledge, type = "pred", terms = c("spring.dev", "cluster"))

## lep sensitivity time
fledge_gdd_pc <- fledge_gdd_pc %>% 
  filter(!is.na(juv_meanday),
         !is.na(spring.gdd),
         !is.na(spring.dev)) %>% 
  dplyr::select(species, cell, year, juv_meanday, 
                spring.gdd, summer.gdd,
                spring.dev, summer.dev, cluster)

fledge_gdd_pc_group <- fledge_gdd_pc %>% 
  group_by(cell, year, cluster) %>% 
  summarise(mean_fledge = mean(juv_meanday))

fledge_gdd_pc_group <- left_join(fledge_gdd_pc_group, gdd)

# spatiotemporal arrival model
sp_fledge <- lmer(mean_fledge ~ spring.gdd + year + spring.gdd:year +
                    cluster + 
                    spring.gdd:cluster +
                    year:cluster +
                    spring.gdd:year:cluster +
                    (1|cell), data = fledge_gdd_pc_group)

summary(sp_fledge)

plot_model(sp_fledge, type = "pred", terms = c("year", "spring.gdd", "cluster"))

# deviation model 
t_fledge <- lmer(mean_fledge ~ spring.dev + year + 
                   cluster + 
                   spring.dev:cluster +
                   year:cluster +
                   (1|cell), data = fledge_gdd_pc_group)

summary(t_fledge)

plot_model(t_fledge, type = "pred", terms = c("year", "cluster"))
plot_model(t_fledge, type = "pred", terms = c("spring.dev", "cluster"))
