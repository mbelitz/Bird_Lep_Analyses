library(dplyr)
library(ggplot2)
library(sf)
library(lme4)
library(lmerTest)
library(sjPlot)
library(merTools)

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

ggsave(filename = "figures/clustering.png")

# read in gdd data
gdd <- data.table::fread('data/gdd_calcs.txt')

# read in frost free perid data
ffp <- read.csv("Outputs/frostFreePeriod_byHex.csv")

# combine phenometrics with gdd
arr_gdd <- left_join(arr, gdd)
fledge_gdd <- left_join(fledge,gdd)
leps_gdd <- left_join(leps, gdd)

# combine phenometrics & gdd with ffp
arr_gdd <- left_join(arr_gdd, ffp)
fledge_gdd <- left_join(fledge_gdd, ffp)
leps_gdd <- left_join(leps_gdd, ffp)

# add PC based clusters to df
arr_gdd_pc <- left_join(arr_gdd, bird_pc)
fledge_gdd_pc <- left_join(fledge_gdd, bird_pc, by = "species")

# look into sensitivty of phenometrics to gdd through space and time
# start with greenup
greenup <- arr_gdd_pc %>% 
  dplyr::distinct(cell, year, .keep_all = T) %>% 
  dplyr::select(cell, year, gr_mn, spring.gdd, summer.gdd, 
                spring.dev, summer.dev, FFD)
greenup <- left_join(greenup, ffp)
greenup$cell <- as.character(greenup$cell)


# deviation greenup model
t_gu <-  lmer(formula = gr_mn ~ spring.dev + year + FFD +
                spring.dev:FFD +
                  (1|cell), data = greenup)

t_gu_s <- step(t_gu)
t_gu_s

summary(t_gu)
car::vif(t_gu)

b <- plot_model(t_gu, type = "pred", terms = "spring.dev")
b
plot_model(t_gu, type = "eff", terms = c("FFD","spring.dev"))
plot_model(t_gu, type = "eff", terms = c("year"))

greenup_spdev <- t_gu@beta[2]
greenup_year <- t_gu@beta[3]

## now predict green up across first a low spring.gdd
length(unique(greenup$cell))

GFFD <- greenup %>% 
  group_by(cell) %>% 
  summarise(meanFFD = mean(FFD))

gu_pdf_low <- data.frame(spring.dev = rep(mean(greenup$spring.dev) - sd(greenup$spring.dev), 73),
                  year = rep(mean(greenup$year),73),
                  cell = GFFD$cell,
                  FFD = GFFD$meanFFD)
gu_pdf_mid <- data.frame(spring.dev = rep(mean(greenup$spring.dev), 73),
                         year = rep(mean(greenup$year),73),
                         cell =  GFFD$cell,
                         FFD = GFFD$meanFFD)
gu_pdf_high <- data.frame(spring.dev = rep(mean(greenup$spring.dev) + sd(greenup$spring.dev), 73),
                         year = rep(mean(greenup$year),73),
                         cell =  GFFD$cell,
                         FFD = GFFD$meanFFD)

gu_pred_low <- predictInterval(t_gu, newdata = gu_pdf_low, which = "full")
gu_pred_mid <- predictInterval(t_gu, newdata = gu_pdf_mid, which = "full")
gu_pred_high <- predictInterval(t_gu, newdata = gu_pdf_high, which = "full")

gu_pred_low_df <- gu_pdf_low  %>% 
  cbind(gu_pred_low)
gu_pred_mid_df <- gu_pdf_mid %>% 
  cbind(gu_pred_mid)
gu_pred_high_df <- gu_pdf_high  %>% 
  cbind(gu_pred_high)

library(sf)
hex_sf <- raster::shapefile("data/hex_grid_crop.shp") %>% 
  st_as_sf()

hex_sf <- st_transform(hex_sf, "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")
hex_sf$cell <- as.character(hex_sf$cell)

gu_pred_high_sf <- left_join(hex_sf, gu_pred_high_df) %>% 
  filter(!is.na(fit))
gu_pred_mid_sf <- left_join(hex_sf, gu_pred_mid_df) %>% 
  filter(!is.na(fit))
gu_pred_low_sf <- left_join(hex_sf, gu_pred_low_df) %>% 
  filter(!is.na(fit))

ggplot() +
  geom_sf(gu_pred_mid_sf, mapping = aes(fill = fit))
ggplot() +
  geom_sf(gu_pred_low_sf, mapping = aes(fill = fit))

# let's do birds next
arr_gdd_pc$cluster <- as.character(arr_gdd_pc$cluster)
arr_gdd_pc <- arr_gdd_pc %>% 
  filter(!is.na(arr_GAM_mean)) %>% 
  dplyr::select(species, cell, year, arr_GAM_mean, 
                arr_IAR_mean, arr_IAR_sd, spring.gdd, summer.gdd,
                spring.dev, summer.dev, cluster, mean_ffp, FFD)

arr_gdd_pc_group <- arr_gdd_pc %>% 
  group_by(cell, year, cluster) %>% 
  summarise(mean_arr = mean(arr_IAR_mean))

arr_gdd_pc_group <- left_join(arr_gdd_pc_group, gdd)
arr_gdd_pc_group <- left_join(arr_gdd_pc_group, ffp)
arr_gdd_pc_group$cluster <- as.character(arr_gdd_pc_group$cluster)

# deviation model 
t_ar1 <- lmer(mean_arr ~ spring.dev + year + FFD +
                spring.dev:cluster +
                year:cluster +
               spring.dev:FFD +
                (1|cell), data = arr_gdd_pc_group)

t_ar2 <- lmer(mean_arr ~ spring.dev + year + FFD +
                spring.dev:FFD:cluster +
                year:cluster +
                spring.dev:FFD +
                (1|cell), data = arr_gdd_pc_group)


summary(t_ar1)
summary(t_ar2)
car::vif(t_ar1)
car::vif(t_ar2)

MuMIn::Weights(AIC(t_ar1,t_ar2))
summary(t_ar1)

arr_spring.dev <-    t_ar1@beta[2]
arr_spring.dev_c2 <- t_ar1@beta[5]
arr_spring.dev_c3 <- t_ar1@beta[6]

plot_model(t_ar1, type = "eff", terms = c("year", "cluster"))
plot_model(t_ar1, type = "eff", terms = c("spring.dev", "FFD"))

## now predict arrival across first a low spring.gdd
source('scripts/functions_toPredictModels.R')

ar_high_c1_sf <- generate_arrival_predicted_sf(lowmidhigh = 1, grouping = 1)
ar_high_c2_sf <- generate_arrival_predicted_sf(lowmidhigh = 1, grouping = 2)
ar_high_c3_sf <- generate_arrival_predicted_sf(lowmidhigh = 1, grouping = 3)

ar_mid_c1_sf <- generate_arrival_predicted_sf(lowmidhigh = 0, grouping = 1)
ar_mid_c2_sf <- generate_arrival_predicted_sf(lowmidhigh = 0, grouping = 2)
ar_mid_c3_sf <- generate_arrival_predicted_sf(lowmidhigh = 0, grouping = 3)

ar_low_c1_sf <- generate_arrival_predicted_sf(lowmidhigh = -1, grouping = 1)
ar_low_c2_sf <- generate_arrival_predicted_sf(lowmidhigh = -1, grouping = 2)
ar_low_c3_sf <- generate_arrival_predicted_sf(lowmidhigh = -1, grouping = 3)

## now let's do bird fledging
fledge_gdd_pc$cluster <- as.character(fledge_gdd_pc$cluster)
fledge_gdd_pc <- fledge_gdd_pc %>% 
  filter(!is.na(juv_meanday),
         !is.na(spring.gdd),
         !is.na(spring.dev)) %>% 
  dplyr::select(species, cell, year, juv_meanday, 
                spring.gdd, summer.gdd,
                spring.dev, summer.dev, cluster, FFD)

fledge_gdd_pc_group <- fledge_gdd_pc %>% 
  group_by(cell, year, cluster) %>% 
  summarise(mean_fledge = mean(juv_meanday))

fledge_gdd_pc_group <- left_join(fledge_gdd_pc_group, gdd)
fledge_gdd_pc_group <- left_join(fledge_gdd_pc_group, ffp)

# deviation arrival model
t_fledge1 <- lmer(mean_fledge ~ spring.dev + year + FFD +
                    spring.dev:cluster +
                    year:cluster +
                    spring.dev:FFD +
                (1|cell), data = fledge_gdd_pc_group)

t_fledge2 <- lmer(mean_fledge ~ spring.dev + year + FFD +
                     spring.dev:cluster +
                     year:cluster +
                     spring.dev:FFD:cluster +
                     (1|cell), data = fledge_gdd_pc_group)

MuMIn::Weights(AIC(t_fledge1, t_fledge2))

summary(t_fledge1)
car::vif(t_fledge1)

plot_model(t_fledge1, type = "eff", terms = c("spring.dev", "FFD"))
plot_model(t_fledge1, type = "eff", terms = c("year", "cluster"))
plot_model(t_fledge1, type = "eff", terms = c("spring.dev", "cluster"))

summary(t_fledge1)

fledge_spring.dev <- t_fledge1@beta[2]
fledge_spring.dev_c2 <- t_fledge1@beta[5]
fledge_spring.dev_c3 <- t_fledge1@beta[6]

## get sf objects of fleding
fledge_high_c1_sf <- generate_fledge_predicted_sf(lowmidhigh = 1, grouping = 1)
fledge_high_c2_sf <- generate_fledge_predicted_sf(lowmidhigh = 1, grouping = 2)
fledge_high_c3_sf <- generate_fledge_predicted_sf(lowmidhigh = 1, grouping = 3)

fledge_mid_c1_sf <- generate_fledge_predicted_sf(lowmidhigh = 0, grouping = 1)
fledge_mid_c2_sf <- generate_fledge_predicted_sf(lowmidhigh = 0, grouping = 2)
fledge_mid_c3_sf <- generate_fledge_predicted_sf(lowmidhigh = 0, grouping = 3)

fledge_low_c1_sf <- generate_fledge_predicted_sf(lowmidhigh = -1, grouping = 1)
fledge_low_c2_sf <- generate_fledge_predicted_sf(lowmidhigh = -1, grouping = 2)
fledge_low_c3_sf <- generate_fledge_predicted_sf(lowmidhigh = -1, grouping = 3)

## lep sensitivity time
leps_gdd <- leps_gdd %>% 
  filter(!is.na(q5),
         !is.na(code),
         !is.na(spring.gdd),
         !is.na(spring.dev)) %>% 
  dplyr::select(code, cell, year, q5, uniqObsDays,
                spring.gdd, summer.gdd,
                spring.dev, summer.dev, FFD)

# deviation model 
t_leps1 <- lmer(q5 ~ spring.dev + year + FFD +
                   uniqObsDays +
                   spring.dev:code +
                   year:code +
                  spring.dev:FFD +
                   (1|cell), data = leps_gdd)

t_leps2 <- lmer(q5 ~ spring.dev + year + FFD +
                  uniqObsDays +
                  spring.dev:code +
                  year:code +
                  spring.dev:FFD:code +
                  (1|cell), data = leps_gdd)

MuMIn::Weights(AIC(t_leps1, t_leps2))

summary(t_leps1)
car::vif(t_leps1)

summary(t_leps1)

plot_model(t_leps1, type = "slope", terms = c("year", "code"))
plot_model(t_leps1, type = "eff", terms = c("spring.dev", "code")) 
plot_model(t_leps1, type = "eff", terms = c("spring.dev", "FFD")) 


## generate sf objects for lep predictions
leps_high_RL_sf <- generate_lep_predicted_sf(lowmidhigh = 1, grouping = "RL")
leps_high_RP_sf <- generate_lep_predicted_sf(lowmidhigh = 1, grouping = "RP")
leps_high_RE_sf <- generate_lep_predicted_sf(lowmidhigh = 1, grouping = "RE")

leps_mid_RL_sf <- generate_lep_predicted_sf(lowmidhigh = 0, grouping = "RL")
leps_mid_RP_sf <- generate_lep_predicted_sf(lowmidhigh = 0, grouping = "RP")
leps_mid_RE_sf <- generate_lep_predicted_sf(lowmidhigh = 0, grouping = "RE")

leps_low_RL_sf <- generate_lep_predicted_sf(lowmidhigh = -1, grouping = "RL")
leps_low_RP_sf <- generate_lep_predicted_sf(lowmidhigh = -1, grouping = "RP")
leps_low_RE_sf <- generate_lep_predicted_sf(lowmidhigh = -1, grouping = "RE")

## plot together
library(ggpubr)

## first lets plot the sensitivities
arr_pred <- predict_arrival_allGDD(1) %>% 
  mutate(Pheno.phase = "Bird arrival")
arr_int <- filter(arr_pred, spring.dev == 0)$fit

ggplot() +
  geom_ribbon(arr_pred, mapping = aes(x = spring.dev, 
                                      ymin = lwr - arr_int, 
                                      ymax = upr - arr_int),
              fill = "blue", alpha = 0.2) +
  geom_smooth(arr_pred, mapping = aes(x = spring.dev,
                                      y = fit - arr_int),
              method = "lm") 

# fledge
fledge_pred_1 <- predict_fledge_allGDD(1) %>% 
  mutate(Pheno.phase = "Fledge")
fledge_int_1 <- filter(fledge_pred_1, spring.dev == 0)$fit
fledge_pred_1  <- fledge_pred_1  %>% 
  mutate(fit2 = fit - fledge_int_1,
         upr2 = upr - fledge_int_1,
         lwr2 = lwr - fledge_int_1)

fledge_pred_2 <- predict_fledge_allGDD(2) %>% 
  mutate(Pheno.phase = "Fledge")
fledge_int_2 <- filter(fledge_pred_2, spring.dev == 0)$fit
fledge_pred_2  <- fledge_pred_2  %>% 
  mutate(fit2 = fit - fledge_int_2,
         upr2 = upr - fledge_int_2,
         lwr2 = lwr - fledge_int_2)

fledge_pred_3 <- predict_fledge_allGDD(3) %>% 
  mutate(Pheno.phase = "Fledge")
fledge_int_3 <- filter(fledge_pred_3, spring.dev == 0)$fit
fledge_pred_3  <- fledge_pred_3  %>% 
  mutate(fit2 = fit - fledge_int_3,
         upr2 = upr - fledge_int_3,
         lwr2 = lwr - fledge_int_3)

# grab green up
gu_pred <- predict_greenup_allGDD(cellz = "725") %>% 
  dplyr::mutate(cluster = "",
                Pheno.phase = "Greenup") %>% 
  dplyr::select(spring.dev, cluster, cell, year, FFD, fit, upr, lwr, Pheno.phase)
gu_int <- filter(gu_pred, spring.dev == 0)$fit
gu_pred <- gu_pred %>% 
  mutate(fit2 = fit - gu_int,
         upr2 = upr - gu_int,
         lwr2 = lwr - gu_int)

# leps
lep_pred_RE <- predict_lep_allGDD("RE") %>% 
  mutate(Pheno.phase = "Emergence") %>% 
  dplyr::rename(cluster = code) %>% 
  dplyr::select(-uniqObsDays)
lep_int_RE <- filter(lep_pred_RE, spring.dev == 0)$fit
lep_pred_RE <- lep_pred_RE %>% 
  mutate(fit2 = fit - lep_int_RE,
         upr2 = upr - lep_int_RE,
         lwr2 = lwr - lep_int_RE)

lep_pred_RP <- predict_lep_allGDD("RP") %>% 
  mutate(Pheno.phase = "Emergence") %>% 
  dplyr::rename(cluster = code) %>% 
  dplyr::select(-uniqObsDays)
lep_int_RP <- filter(lep_pred_RP, spring.dev == 0)$fit
lep_pred_RP <- lep_pred_RP %>% 
  mutate(fit2 = fit - lep_int_RP,
         upr2 = upr - lep_int_RP,
         lwr2 = lwr - lep_int_RP)

lep_pred_RL <- predict_lep_allGDD("RL") %>% 
  mutate(Pheno.phase = "Emergence") %>% 
  dplyr::rename(cluster = code) %>% 
  dplyr::select(-uniqObsDays)
lep_int_RL <- filter(lep_pred_RL, spring.dev == 0)$fit
lep_pred_RL <- lep_pred_RL %>% 
  mutate(fit2 = fit - lep_int_RL,
         upr2 = upr - lep_int_RL,
         lwr2 = lwr - lep_int_RL)

## rbind this all together
tdf <- rbind(gu_pred,
             fledge_pred_1, fledge_pred_2, fledge_pred_3,
             lep_pred_RE, lep_pred_RL, lep_pred_RP) %>% 
  mutate(id = paste(Pheno.phase, cluster))


# plot sensitivity
ggplot() +
  geom_ribbon(tdf, mapping = aes(x = spring.dev, 
                                     ymin = lwr2, 
                                     ymax = upr2,
                                     fill = id),
              alpha = 0.2) +
  geom_smooth(tdf, mapping = aes(x = spring.dev,
                                     y = fit2,
                                     color = id),
               method = "lm") +
  scale_color_manual(values = c(
    rev(c("#3E92CC", "#2A628F", "#13293D")),
    "#F95738", "#EE964B", "#F4D35E",
    "forest green"
  )) +
  scale_fill_manual(values = c(
    rev(c("#3E92CC", "#2A628F", "#13293D")),
    "#F95738", "#EE964B", "#F4D35E",
    "forest green"
  )) +
  theme_bw()











x = c(-100:100)

ggplot() +
  geom_line(mapping = aes(x = x, y = greenup_sens*x, color = "Greenup"), 
            size = 1.1, linetype = "dashed") + 
  geom_line(mapping = aes(x = x, y = arr_spring.dev*x, color = "Arrival (C1)"), 
            size = 1.1) + 
  geom_line(mapping = aes(x = x, y = arr_spring.dev*x + arr_spring.dev_c2*x, color = "Arrival (C2)"),
            size = 1.1) + 
  geom_line(mapping = aes(x = x, y = arr_spring.dev*x + arr_spring.dev_c3*x, color = "Arrival (C3)"),
             size = 1.1) + 
  geom_line(mapping = aes(x = x, y = lep_spring.dev*x, color = "Emergence (Egg)"), 
            size = 1.1) + 
  geom_line(mapping = aes(x = x, y = lep_spring.dev*x + lep_spring.dev_RL*x, color = "Emergence (Larvae)"),
            size = 1.1) + 
  geom_line(mapping = aes(x = x, y = lep_spring.dev*x + lep_spring.dev_RP*x, color = "Emergence (Pupae)"),
            size = 1.1) + 
  labs(x = "GDD deviation", y = "Phenological sensitivity", color = "") + 
  scale_color_manual(values = 
                       c(
                         "#F95738", "#EE964B", "#F4D35E",
                         rev(c("#3E92CC", "#2A628F", "#13293D")),
    "forest green")) +
  theme_bw()

ggsave(filename = "figures/arrival_sensitivity.png", plot = arr_sens, width = 7.5, height = 4)


fledge_sens <- ggplot() +
  geom_line(mapping = aes(x = x, y = greenup_sens*x, color = "Greenup"), 
            size = 1.1, linetype = "dashed") + 
  geom_line(mapping = aes(x = x, y = fledge_spring.dev*x, color = "Fledge (C1)"), 
            size = 1.1) + 
  geom_line(mapping = aes(x = x, y = fledge_spring.dev*x + fledge_spring.dev_c2*x, color = "Fledge (C2)"),
            size = 1.1) + 
  geom_line(mapping = aes(x = x, y = fledge_spring.dev*x + fledge_spring.dev_c3*x, color = "Fledge (C3)"),
            size = 1.1) + 
  geom_line(mapping = aes(x = x, y = lep_spring.dev*x, color = "Emergence (Egg)"), 
            size = 1.1) + 
  geom_line(mapping = aes(x = x, y = lep_spring.dev*x + lep_spring.dev_RL*x, color = "Emergence (Larvae)"),
            size = 1.1) + 
  geom_line(mapping = aes(x = x, y = lep_spring.dev*x + lep_spring.dev_RP*x, color = "Emergence (Pupae)"),
            size = 1.1) + 
  labs(x = "GDD deviation", y = "Phenological sensitivity", color = "") + 
  scale_color_manual(values = c(
    rev(c("#3E92CC", "#2A628F", "#13293D")),
    "#F95738", "#EE964B", "#F4D35E",
    "forest green"
                       )) +
  theme_bw()

ggsave(filename = "figures/fledge_sensitivity.png", plot = fledge_sens, width = 7.5, height = 4)

## get resids
lep_r <- resid(t_leps1)
leps_resids <- leps_gdd %>%
mutate(resids = lep_r)

fledge_r <- resid(t_fledge1)
fledge_resids <- as.data.frame(fledge_gdd_pc_group) %>%
mutate(resids = fledge_r)

green_r <- resid(t_gu)
green_resids <- as.data.frame(greenup) %>%
mutate(resids = green_r)

#Temporal residuals
leps_r_p <- ggplot(leps_resids, mapping = aes(x = year, y = resids)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()

arrival_r_p <- ggplot(arr_resids, mapping = aes(x = year, y = resids)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()

fledge_r_p <- ggplot(fledge_resids, mapping = aes(x = year, y = resids)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()

green_r_p <- ggplot(green_resids, mapping = aes(x = year, y = resids)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()

# spatial residuals
library(sf)
h <- raster::shapefile("data/hex_grid_crop.shp")

h_sf <- st_as_sf(h) %>% 
  st_transform(crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
h_sf$cell <- as.integer(h_sf$cell)

## leps
lep_h_resids <- left_join(h_sf, leps_resids) %>% 
  dplyr::filter(!is.na(resids))

lep_sp_resid <- ggplot(data = lep_h_resids) +
  geom_sf(mapping = aes(fill = resids)) +
  scale_fill_gradient2() +
  theme_classic() 

# birds
arr_h_resids <- left_join(h_sf, arr_resids) %>% 
  dplyr::filter(!is.na(resids))

arr_sp_res <- ggplot(data = arr_h_resids) +
  geom_sf(mapping = aes(fill = resids)) +
  scale_fill_gradient2() +
  theme_classic() 

fledge_h_resids <- left_join(h_sf, fledge_resids) %>% 
  dplyr::filter(!is.na(resids))

fledge_sp_res <- ggplot(data = fledge_h_resids) +
  geom_sf(mapping = aes(fill = resids)) +
  scale_fill_gradient2() +
  theme_classic() 

# greenup
green_resids$cell <- as.integer(green_resids$cell)
green_h_resids <- left_join(h_sf, green_resids) %>% 
  dplyr::filter(!is.na(resids))

green_sp_res <- ggplot(data = green_h_resids) +
  geom_sf(mapping = aes(fill = resids)) +
  scale_fill_gradient2() +
  theme_classic() 

cp <- cowplot::plot_grid(green_sp_res, lep_sp_resid,
                         arr_sp_res, fledge_sp_res,
                         labels = c("Greenup", "Butterfly emergence",
                                    "Bird arrival", "Bird fledge"))

cp
ggsave(plot = cp, filename = "figures/spatial_residuals.png", width = 12, height = 8)

## now we are plotting spatial locations of estimates in regards to greenup
arr_fig_c1_high <- left_join(ar_high_c1_sf, gu_pred_high_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
arr_fig_c1_high_sf <- left_join(hex_sf, arr_fig_c1_high) %>% 
  filter(!is.na(Difference))

p1 <- ggplot() +
  geom_sf(arr_fig_c1_high_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("C1 Arrival, Warm year") +
  theme_bw()

arr_fig_c2_high <- left_join(ar_high_c2_sf, gu_pred_high_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
arr_fig_c2_high_sf <- left_join(hex_sf, arr_fig_c2_high) %>% 
  filter(!is.na(Difference))

p2 <- ggplot() +
  geom_sf(arr_fig_c2_high_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("C2 Arrival, Warm year") +
  theme_bw()

arr_fig_c3_high <- left_join(ar_high_c3_sf, gu_pred_high_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
arr_fig_c3_high_sf <- left_join(hex_sf, arr_fig_c3_high) %>% 
  filter(!is.na(Difference))

p3 <- ggplot() +
  geom_sf(arr_fig_c3_high_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("C3 Arrival, Warm year") +
  theme_bw()

# now move on to fledge
fledge_fig_c1_high <- left_join(fledge_high_c1_sf, gu_pred_high_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
fledge_fig_c1_high_sf <- left_join(hex_sf, fledge_fig_c1_high) %>% 
  filter(!is.na(Difference))

p4 <- ggplot() +
  geom_sf(fledge_fig_c1_high_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("C1 Fledge, Warm year") +
  theme_bw()

fledge_fig_c2_high <- left_join(fledge_high_c2_sf, gu_pred_high_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
fledge_fig_c2_high_sf <- left_join(hex_sf, fledge_fig_c2_high) %>% 
  filter(!is.na(Difference))

p5 <- ggplot() +
  geom_sf(fledge_fig_c2_high_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("C2 Fledge, Warm year") +
  theme_bw()

fledge_fig_c3_high <- left_join(fledge_high_c3_sf, gu_pred_high_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
fledge_fig_c3_high_sf <- left_join(hex_sf, fledge_fig_c3_high) %>% 
  filter(!is.na(Difference))

p6 <- ggplot() +
  geom_sf(fledge_fig_c3_high_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("C3 Fledge, Warm year") +
  theme_bw()

## moving on to leps
leps_fig_RE_high <- left_join(leps_high_RE_sf, gu_pred_high_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
leps_fig_RE_high_sf <- left_join(hex_sf, leps_fig_RE_high) %>% 
  filter(!is.na(Difference))%>% 
  filter(Difference < 100)

p7 <- ggplot() +
  geom_sf(leps_fig_RE_high_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("Egg, Warm year") +
  theme_bw()

leps_fig_RP_high <- left_join(leps_high_RP_sf, gu_pred_high_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
leps_fig_RP_high_sf <- left_join(hex_sf, leps_fig_RP_high) %>% 
  filter(!is.na(Difference)) %>% 
  filter(Difference < 100)

p8 <- ggplot() +
  geom_sf(leps_fig_RP_high_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("Pupae, Warm year") +
  theme_bw()

leps_fig_RL_high <- left_join(leps_high_RL_sf, gu_pred_high_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
leps_fig_RL_high_sf <- left_join(hex_sf, leps_fig_RL_high) %>% 
  filter(!is.na(Difference))%>% 
  filter(Difference < 100)

p9 <- ggplot() +
  geom_sf(leps_fig_RL_high_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("Larvae, Warm year") +
  theme_bw()

library(cowplot)
cp <- plot_grid(p1, p2, p3,
                p4, p5, p6,
                p7, p8, p9,
                nrow = 3, ncol = 3)

ggsave(filename = "figures/warmYear_difference.png", plot = cp,
       width = 12, height = 12)




## do spatial difference but for mid
arr_fig_c1_mid <- left_join(ar_mid_c1_sf, gu_pred_mid_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
arr_fig_c1_mid_sf <- left_join(hex_sf, arr_fig_c1_mid) %>% 
  filter(!is.na(Difference))

p1_mid <- ggplot() +
  geom_sf(arr_fig_c1_mid_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("C1 Arrival, Average year") +
  theme_bw()

arr_fig_c2_mid <- left_join(ar_mid_c2_sf, gu_pred_mid_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
arr_fig_c2_mid_sf <- left_join(hex_sf, arr_fig_c2_mid) %>% 
  filter(!is.na(Difference))

p2_mid <- ggplot() +
  geom_sf(arr_fig_c2_mid_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("C2 Arrival, Average year") +
  theme_bw()

arr_fig_c3_mid <- left_join(ar_mid_c3_sf, gu_pred_mid_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
arr_fig_c3_mid_sf <- left_join(hex_sf, arr_fig_c3_mid) %>% 
  filter(!is.na(Difference))

p3_mid <- ggplot() +
  geom_sf(arr_fig_c3_mid_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("C3 Arrival, Average year") +
  theme_bw()

# now move on to fledge
fledge_fig_c1_mid <- left_join(fledge_mid_c1_sf, gu_pred_mid_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
fledge_fig_c1_mid_sf <- left_join(hex_sf, fledge_fig_c1_mid) %>% 
  filter(!is.na(Difference))

p4_mid <- ggplot() +
  geom_sf(fledge_fig_c1_mid_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("C1 Fledge, Average year") +
  theme_bw()

fledge_fig_c2_mid <- left_join(fledge_mid_c2_sf, gu_pred_mid_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
fledge_fig_c2_mid_sf <- left_join(hex_sf, fledge_fig_c2_mid) %>% 
  filter(!is.na(Difference))

p5_mid <- ggplot() +
  geom_sf(fledge_fig_c2_mid_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("C2 Fledge, Average year") +
  theme_bw()

fledge_fig_c3_mid <- left_join(fledge_mid_c3_sf, gu_pred_mid_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
fledge_fig_c3_mid_sf <- left_join(hex_sf, fledge_fig_c3_mid) %>% 
  filter(!is.na(Difference))

p6_mid <- ggplot() +
  geom_sf(fledge_fig_c3_mid_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("C3 Fledge, Average year") +
  theme_bw()

## moving on to leps
leps_fig_RE_mid <- left_join(leps_mid_RE_sf, gu_pred_mid_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
leps_fig_RE_mid_sf <- left_join(hex_sf, leps_fig_RE_mid) %>% 
  filter(!is.na(Difference))%>% 
  filter(Difference < 100)

p7_mid <- ggplot() +
  geom_sf(leps_fig_RE_mid_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("Egg, Average year") +
  theme_bw()

leps_fig_RP_mid <- left_join(leps_mid_RP_sf, gu_pred_mid_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
leps_fig_RP_mid_sf <- left_join(hex_sf, leps_fig_RP_mid) %>% 
  filter(!is.na(Difference)) %>% 
  filter(Difference < 100)

p8_mid <- ggplot() +
  geom_sf(leps_fig_RP_mid_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("Pupae, Average year") +
  theme_bw()

leps_fig_RL_mid <- left_join(leps_mid_RL_sf, gu_pred_mid_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
leps_fig_RL_mid_sf <- left_join(hex_sf, leps_fig_RL_mid) %>% 
  filter(!is.na(Difference))%>% 
  filter(Difference < 100)

p9_mid <- ggplot() +
  geom_sf(leps_fig_RL_mid_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("Larvae, Average year") +
  theme_bw()




## do spatial difference but for low
arr_fig_c1_low <- left_join(ar_low_c1_sf, gu_pred_low_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
arr_fig_c1_low_sf <- left_join(hex_sf, arr_fig_c1_low) %>% 
  filter(!is.na(Difference))

p1_low <- ggplot() +
  geom_sf(arr_fig_c1_low_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("C1 Arrival, Cold year") +
  theme_bw()

arr_fig_c2_low <- left_join(ar_low_c2_sf, gu_pred_low_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
arr_fig_c2_low_sf <- left_join(hex_sf, arr_fig_c2_low) %>% 
  filter(!is.na(Difference))

p2_low <- ggplot() +
  geom_sf(arr_fig_c2_low_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("C2 Arrival, Cold year") +
  theme_bw()

arr_fig_c3_low <- left_join(ar_low_c3_sf, gu_pred_low_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
arr_fig_c3_low_sf <- left_join(hex_sf, arr_fig_c3_low) %>% 
  filter(!is.na(Difference))

p3_low <- ggplot() +
  geom_sf(arr_fig_c3_low_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("C3 Arrival, Cold year") +
  theme_bw()

# now move on to fledge
fledge_fig_c1_low <- left_join(fledge_low_c1_sf, gu_pred_low_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
fledge_fig_c1_low_sf <- left_join(hex_sf, fledge_fig_c1_low) %>% 
  filter(!is.na(Difference))

p4_low <- ggplot() +
  geom_sf(fledge_fig_c1_low_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("C1 Fledge, Cold year") +
  theme_bw()

fledge_fig_c2_low <- left_join(fledge_low_c2_sf, gu_pred_low_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
fledge_fig_c2_low_sf <- left_join(hex_sf, fledge_fig_c2_low) %>% 
  filter(!is.na(Difference))

p5_low <- ggplot() +
  geom_sf(fledge_fig_c2_low_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("C2 Fledge, Cold year") +
  theme_bw()

fledge_fig_c3_low <- left_join(fledge_low_c3_sf, gu_pred_low_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
fledge_fig_c3_low_sf <- left_join(hex_sf, fledge_fig_c3_low) %>% 
  filter(!is.na(Difference))

p6_low <- ggplot() +
  geom_sf(fledge_fig_c3_low_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("C3 Fledge, Cold year") +
  theme_bw()

## moving on to leps
leps_fig_RE_low <- left_join(leps_low_RE_sf, gu_pred_low_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
leps_fig_RE_low_sf <- left_join(hex_sf, leps_fig_RE_low) %>% 
  filter(!is.na(Difference))%>% 
  filter(Difference < 100)

p7_low <- ggplot() +
  geom_sf(leps_fig_RE_low_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("Egg, Cold year") +
  theme_bw()

leps_fig_RP_low <- left_join(leps_low_RP_sf, gu_pred_low_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
leps_fig_RP_low_sf <- left_join(hex_sf, leps_fig_RP_low) %>% 
  filter(!is.na(Difference)) %>% 
  filter(Difference < 100)

p8_low <- ggplot() +
  geom_sf(leps_fig_RP_low_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("Pupae, Cold year") +
  theme_bw()

leps_fig_RL_low <- left_join(leps_low_RL_sf, gu_pred_low_df, by = "cell") %>% 
  mutate(Difference = fit.x - fit.y) 
leps_fig_RL_low_sf <- left_join(hex_sf, leps_fig_RL_low) %>% 
  filter(!is.na(Difference))%>% 
  filter(Difference < 100)

p9_low <- ggplot() +
  geom_sf(leps_fig_RL_low_sf, mapping = aes(fill = Difference)) +
  scale_fill_gradient2() +
  ggtitle("Larvae, Cold year") +
  theme_bw()


library(cowplot)
cp <- plot_grid(p1_low, p1_mid, p1,
                p4_low, p4_mid, p4,
                nrow = 2, ncol = 3)

ggsave(filename = "figures/bird_C1_difference.png", plot = cp,
       width = 12, height = 8)


cp <- plot_grid(p2_low, p2_mid, p2,
                p5_low, p5_mid, p5,
                nrow = 2, ncol = 3)

ggsave(filename = "figures/bird_C2_difference.png", plot = cp,
       width = 12, height = 8)


cp <- plot_grid(p3_low, p3_mid, p3,
                p6_low, p6_mid, p6,
                nrow = 2, ncol = 3)

ggsave(filename = "figures/bird_C3_difference.png", plot = cp,
       width = 12, height = 8)

cp <- plot_grid(p7_low, p7_mid, p7,
                p8_low, p8_mid, p8,
                p9_low, p9_mid, p9,
                nrow = 3, ncol = 3)

ggsave(filename = "figures/leps_difference.png", plot = cp,
       width = 12, height = 8)


# plotting the delta values

rest <- st_join(leps_fig_RE_high_sf, leps_fig_RE_low_sf, by = "cell") %>% 
  mutate(delta = Difference.y - Difference.x)

re_delta <- ggplot() +
  geom_sf(rest, mapping = aes(fill = delta)) +
  scale_fill_gradient2(limits = c(-50,37)) +
  ggtitle("Eggs") +
  theme_bw()

rlst <- st_join(leps_fig_RL_high_sf, leps_fig_RL_low_sf, by = "cell") %>% 
  mutate(delta = Difference.y - Difference.x)

rl_delta <- ggplot() +
  geom_sf(rlst, mapping = aes(fill = delta)) +
  scale_fill_gradient2(limits = c(-50,37)) +
  ggtitle("Larvae") +
  theme_bw()

rpst <- st_join(leps_fig_RP_high_sf, leps_fig_RP_low_sf, by = "cell") %>% 
  mutate(delta = Difference.y - Difference.x)

rp_delta <- ggplot() +
  geom_sf(rpst, mapping = aes(fill = delta)) +
  scale_fill_gradient2(limits = c(-50,37)) +
  ggtitle("Pupae") +
  theme_bw()

## bird time
f1st <- st_join(fledge_fig_c1_high_sf, fledge_fig_c1_low_sf, by = "cell") %>% 
  mutate(delta = Difference.y - Difference.x)

f1_delta <- ggplot() +
  geom_sf(f1st, mapping = aes(fill = delta)) +
  scale_fill_gradient2(limits = c(-50,37)) +
  ggtitle("Fledge (C1)") +
  theme_bw()


f2st <- st_join(fledge_fig_c2_high_sf, fledge_fig_c2_low_sf, by = "cell") %>% 
  mutate(delta = Difference.y - Difference.x)

f2_delta <- ggplot() +
  geom_sf(f2st, mapping = aes(fill = delta)) +
  scale_fill_gradient2(limits = c(-50,37)) +
  ggtitle("Fledge (C2)") +
  theme_bw()

f3st <- st_join(fledge_fig_c3_high_sf, fledge_fig_c3_low_sf, by = "cell") %>% 
  mutate(delta = Difference.y - Difference.x)

f3_delta <- ggplot() +
  geom_sf(f3st, mapping = aes(fill = delta)) +
  scale_fill_gradient2(limits = c(-50,37)) +
  ggtitle("Fledge (C3)") +
  theme_bw()

## combine birds and leps

ga <- ggpubr::ggarrange(f1_delta, f2_delta, f3_delta,
                         re_delta, rl_delta, rp_delta,
                         nrow = 2 , ncol = 3,
                        common.legend = TRUE, legend = "right")

ga

ggsave(ga, filename = "figures/delta_phenophase.png")

## grenup time
us <- rnaturalearth::ne_countries(country = "United States of America", returnclass = "sf") %>% 
  st_transform(crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

gust <- st_join(gu_pred_high_sf, gu_pred_low_sf, by = "cell") %>% 
  mutate(delta = Difference.y - Difference.x)

