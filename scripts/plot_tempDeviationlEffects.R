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
# caclulate baseline from 2016-2016
temp_baselineYears <- temp %>% 
  filter(year %in% 2016:2019)
tempBaseline <- temp_baselineYears %>% 
  group_by(cell) %>% 
  summarise(baseline_temp = mean(mean_temp))
temp <- left_join(temp, tempBaseline)
temp <- temp %>% 
  mutate(temp.dev = mean_temp - baseline_temp)
# combine phenometrics with gdd & ffp
arr_gdd <- left_join(arr, gdd)
arr_gdd <- left_join(arr_gdd, ffp)
arr_gdd <- left_join(arr_gdd, temp)
fledge_gdd <- left_join(fledge,gdd)
fledge_gdd <- left_join(fledge_gdd, ffp) %>% 
  filter(!is.na(mean_ffp))
fledge_gdd <- left_join(fledge_gdd, temp)
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
                summer.dev, mean_ffp, temp.dev) %>% 
  dplyr::filter(year >= 2002 & year <= 2017)
greenup$cell <- as.character(greenup$cell)
# scale greenup dataframe
greenup_scaled <- greenup %>% 
  mutate(temp.dev = scale(temp.dev),
         mean_ffp = scale(mean_ffp),
         temp.dev = scale(temp.dev))

# deviation greenup model
gu_m <-  lmer(formula = gr_mn ~ temp.dev + mean_ffp +
                temp.dev:mean_ffp +
                (1|cell), data = greenup_scaled,
              na.action = na.fail, REML = F)

gu_d <- dredge(gu_m)

gu_tm <- get.models(gu_d, subset = 1)[[1]]
gu_tm

summary(gu_tm)
car::vif(gu_tm)

MuMIn::r.squaredGLMM(gu_tm)


plot_model(gu_tm, type = "pred", terms = "temp.dev")
plot_model(gu_tm, type = "pred", terms = c("mean_ffp"))
plot_model(gu_tm, type = "pred", terms = c("temp.dev", "mean_ffp"))

gu_gdd_plot <- plot_model(gu_tm, type = "pred", terms = c("temp.dev", "mean_ffp"))

gu_gdd_plot +
  theme_classic() +
  labs(x = "Temperature anomoly", y = "Greenup") +
  ggtitle("")

# get ACF for green up plot
resids <- residuals(gu_tm)
resid_greenup <- greenup_scaled %>% 
  mutate(resid = resids)

gu_resid_plot <- ggplot(resid_greenup, mapping = aes(x = year, y = resid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Year", y = "Residual")+
  ggtitle("Greenup") +
  theme_bw()

resid_greenup_sort <- resid_greenup[order(resid_greenup$year),]

acfPlot <- acf(resid_greenup_sort$resid) 

gu_acf_plot <- ggAcf(resid_greenup_sort$resid) +
  theme_bw() +
  ggtitle("")

## spatial autocorrelation and qqplot of green-up model
hist(resids)
qqnorm(resids)
qqline(resids)

# read in 
r <- raster("data/bioclim/Normal_1991_2020_FFP.tif")
hex_sf <- raster::shapefile("data/hex_grid_crop.shp") %>% 
  st_as_sf()
hex_sf_ea <- st_transform(hex_sf, crs = crs(r))
hex_coords <- st_coordinates(st_centroid(hex_sf_ea))
hex_coords_df <- data.frame(cell = hex_sf_ea$cell, 
                            X = hex_coords[,1],
                            Y = hex_coords[,2])

resid_greenup <- left_join(resid_greenup, hex_coords_df)

library(ncf)
ncf_fit <- correlog(x = resid_greenup$X, y = resid_greenup$Y, z = resid_greenup$resid,
                    increment = 250000, latlon = F, resamp = 100)
plot(ncf_fit)

ncf_gu_df <- data.frame(distance = ncf_fit$mean.of.class, 
                      correlation = ncf_fit$correlation,
                      p.value = if_else(ncf_fit$p < 0.025, 
                                        true = "Sig",
                                        false = "Not Sig")) %>% 
  filter(distance < 5000000 & distance > 1)


gu_correlogram_plot <- ggplot() +
  geom_line(ncf_gu_df,
            mapping = aes(x = distance, y = correlation)) +
  geom_point(ncf_gu_df,
             mapping = aes(x = distance, y = correlation, fill = p.value),
             size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("white", "black")) +
  labs(x = "Distance", y = "Moran's I", fill = "P-Value") +
  theme_classic()

gu_correlogram_plot # no autocorrelation found at closest spatial lag, so spatial model not made

## now get this for fledge and leps
fledge_gdd <- fledge_gdd %>% 
  filter(!is.na(juv_meanday),
         !is.na(spring.gdd),
         !is.na(temp.dev)) %>% 
  dplyr::select(sci_name, station, year, cell, juv_meanday, 
                spring.gdd, summer.gdd,
                temp.dev, summer.dev, PC1, PC2, species, mean_ffp) %>% 
  filter(year >= 2002 & year <= 2017)

head(fledge_gdd)

fledge_scaled <- fledge_gdd %>% 
  mutate(temp.dev = scale(temp.dev),
         mean_ffp = scale(mean_ffp)) %>% 
  filter(!is.na(juv_meanday),
         !is.na(spring.gdd),
         !is.na(temp.dev),
         !is.na(year)) 

# deviation fledge model
fledge_m <- lmer(juv_meanday ~ temp.dev + mean_ffp +
                   PC1 +
                   temp.dev:mean_ffp +
                   temp.dev:PC1 +
                   (1|station) + (1|sci_name),
                 data = fledge_scaled, na.action = na.fail, REML = F)

#dredge
fledge_d <- dredge(fledge_m)
tm_fledge <- get.models(fledge_d, subset = 1)[[1]]

summary(tm_fledge)
car::vif(tm_fledge)
MuMIn::r.squaredGLMM(tm_fledge)

plot_model(tm_fledge, type = "pred", terms = "temp.dev")
plot_model(tm_fledge, type = "pred", terms = c("PC1"))
plot_model(tm_fledge, type = "pred", terms = c("mean_ffp"))
plot_model(tm_fledge, type = "pred", terms = c("temp.dev", "mean_ffp"))

fledge_gdd_plot <- plot_model(tm_fledge, type = "pred", terms = c("temp.dev", "mean_ffp"))

fledge_gdd_plot +
  theme_classic() +
  labs(x = "GDD anomaly", y = "Greenup") +
  ggtitle("")

fp_int <- plot_model(tm_fledge, type = "pred", terms = c("temp.dev", "mean_ffp"))

fp_int +
  theme_classic() +
  labs(x = "GDD anomaly", y = "Fledge",
       color = "Frost free period", fill = "Frost free period") +
  ggtitle("")

ggsave("figures/ffp_temp.dev_greenup.png", width = 6, height = 3.5)


# there is an important interaction between mean_ffp and temp.dev

# get ACF for fledge plot
f_resid <- residuals(tm_fledge) %>% as.numeric()
resid_fledge <- ungroup(fledge_scaled) %>% 
  mutate(resid = f_resid)

fledge_resid_plot <- ggplot(resid_fledge, mapping = aes(x = year, y = resid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Year", y = "Residual")+
  ggtitle("Fledge") +
  theme_bw()

resid_fledge_sort <- resid_fledge[order(resid_fledge$year),]

acfPlot <- acf(resid_fledge_sort$resid) 

fledge_acf_plot <- ggAcf(resid_fledge_sort$resid) +
  theme_bw() +
  ggtitle("")

## spatial autocorrelation and qqplot of green-up model
hist(f_resid)
qqnorm(f_resid)
qqline(f_resid)

resid_fledge$cell <- as.character(resid_fledge$cell)
resid_fledge <- left_join(resid_fledge, hex_coords_df)

library(ncf)
ncf_fit <- correlog(x = resid_fledge$X, y = resid_fledge$Y, z = resid_fledge$resid,
                    increment = 250000, latlon = F, resamp = 100)
plot(ncf_fit)

ncf_fledge_df <- data.frame(distance = ncf_fit$mean.of.class, 
                        correlation = ncf_fit$correlation,
                        p.value = if_else(ncf_fit$p < 0.025, 
                                          true = "Sig",
                                          false = "Not Sig")) %>% 
  filter(distance < 5000000 & distance > 1)


fledge_correlogram_plot <- ggplot() +
  geom_line(ncf_fledge_df,
            mapping = aes(x = distance, y = correlation)) +
  geom_point(ncf_fledge_df,
             mapping = aes(x = distance, y = correlation, fill = p.value),
             size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("white", "black")) +
  labs(x = "Distance", y = "Moran's I", fill = "P-Value") +
  theme_classic()

fledge_correlogram_plot # no autocorrelation found at closest spatial lag, so spatial model not made


## now lep time
# deviation model 
leps_gdd <- leps_gdd %>% 
  filter(!is.na(temp.dev),
         !is.na(year),
         !is.na(mean_ffp),
         !is.na(year),
         !is.na(cell),
         !is.na(code),
         !is.na(q5)) %>% 
  filter(year >= 2002 & year <= 2017)

leps_scaled <- leps_gdd %>% 
  mutate(temp.dev = scale(temp.dev),
         mean_ffp = scale(mean_ffp),
         uniqObsDays = scale(uniqObsDays)) 

# note that we can't look at the interaction between code and year because there 
# is not enough RE data in early years
with(leps_scaled, table(code, year))
with(leps_scaled, table(code, temp.dev))

leps_m <- lmer(q5 ~ temp.dev + mean_ffp +
                 uniqObsDays + code +
                 temp.dev:code +
                 temp.dev:mean_ffp +
                 (1|cell), 
               data = leps_scaled, na.action = na.fail, REML = F)

leps_d <- dredge(leps_m)
leps_tm <- get.models(leps_d, 3)[[1]]

summary(leps_tm)

car::vif(leps_tm)
summary(leps_tm)
MuMIn::r.squaredGLMM(leps_tm)

plot_model(leps_tm, type = "pred", terms= "temp.dev")
plot_model(leps_tm, type = "pred", terms= c("mean_ffp"))

leps_gdd_plot <- plot_model(leps_tm, type = "pred", terms = c("temp.dev", "mean_ffp")) 

# get ACF for leps plot
resids_l <- residuals(leps_tm)
resid_leps <- leps_scaled %>% 
  mutate(resid = resids_l)

leps_resid_plot <- ggplot(resid_leps, mapping = aes(x = year, y = resid)) +
  geom_point(aes(color = code)) +
  geom_smooth(aes(color = code), method = "lm") +
  labs(x = "Year", y = "Residual")+
  theme_bw()

#no grouping
leps_resid_plot <- ggplot(resid_leps, mapping = aes(x = year, y = resid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Year", y = "Residual")+
  ggtitle('Leps') +
  theme_bw()

resid_leps_sort <- resid_leps[order(resid_leps$year),]

acfPlot <- acf(resid_leps_sort$resid) 

leps_acf_plot <- ggAcf(resid_leps_sort$resid) +
  theme_bw() +
  ggtitle("")

## spatial autocorrelation and qqplot of green-up model
hist(resids_l)
qqnorm(resids_l)
qqline(resids_l)

resids_l$cell <- as.character(resids_l$cell)
resid_leps$cell <- as.character(resid_leps$cell)
resid_leps <- left_join(resid_leps, hex_coords_df)

library(ncf)
ncf_fit <- correlog(x = resid_leps$X, y = resid_leps$Y, z = resid_leps$resid,
                    increment = 250000, latlon = F, resamp = 100)
plot(ncf_fit)

ncf_leps_df <- data.frame(distance = ncf_fit$mean.of.class, 
                            correlation = ncf_fit$correlation,
                            p.value = if_else(ncf_fit$p < 0.025, 
                                              true = "Sig",
                                              false = "Not Sig")) %>% 
  filter(distance < 5000000 & distance > 1)


leps_correlogram_plot <- ggplot() +
  geom_line(ncf_leps_df,
            mapping = aes(x = distance, y = correlation)) +
  geom_point(ncf_leps_df,
             mapping = aes(x = distance, y = correlation, fill = p.value),
             size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("white", "black")) +
  labs(x = "Distance", y = "Moran's I", fill = "P-Value") +
  theme_classic()

leps_correlogram_plot # no autocorrelation found at closest spatial lag, so spatial model not made


## plot correlograms
cp <- cowplot::plot_grid(gu_correlogram_plot, leps_correlogram_plot,fledge_correlogram_plot, 
                         labels = c("A", "B", "C", nrow = 2, ncol = 2))

ggsave(filename = "figures/correlelogram.png", plot = cp, width = 8, height = 5)


fledge_gdd_plot 
leps_gdd_plot
gu_gdd_plot

# combine temporal effects plot into single data frame
tdf <- list(
  data.frame(leps_gdd_plot$data, line = "leps"),
  data.frame(gu_gdd_plot$data, line = "gu"),
  data.frame(fledge_gdd_plot$data, line = "fledge")
) %>% 
  bind_rows()

tdf <- tdf %>% 
  mutate(group2 = rep(-1:1,24))

lep_int_north <- dplyr::filter(tdf, x == 0 & group2 == -1, line == "leps")$predicted
lep_int_mid <- dplyr::filter(tdf, x == 0 & group2 == 0, line == "leps")$predicted
lep_int_south <- dplyr::filter(tdf, x == 0 & group2 == 1, line == "leps")$predicted

fledge_int_north <- dplyr::filter(tdf, x == 0 & group2 == -1, line == "fledge")$predicted
fledge_int_mid <- dplyr::filter(tdf, x == 0 & group2 == 0, line == "fledge")$predicted
fledge_int_south <- dplyr::filter(tdf, x == 0 & group2 == 1, line == "fledge")$predicted

gu_int_north <- dplyr::filter(tdf, x == 0 & group2 == -1, line == "gu")$predicted
gu_int_mid<- dplyr::filter(tdf, x == 0 & group2 == 0, line == "gu")$predicted
gu_int_south <- dplyr::filter(tdf, x == 0 & group2 == 1, line == "gu")$predicted



tdf <- tdf %>% 
  mutate(predicted2 = case_when(
    line == "gu" & group2 == -1 ~ predicted - gu_int_north,
    line == "gu" & group2 == 0 ~ predicted - gu_int_mid,
    line == "gu" & group2 == 1 ~ predicted - gu_int_south,
    
    line == "fledge" & group2 == -1 ~ predicted - fledge_int_north,
    line == "fledge" & group2 == 0 ~ predicted - fledge_int_mid,
    line == "fledge" & group2 == 1 ~ predicted - fledge_int_south,
    
    line == "leps" & group2 == -1 ~ predicted - lep_int_north,
    line == "leps" & group2 == 0 ~ predicted - lep_int_mid,
    line == "leps" & group2 == 1 ~ predicted - lep_int_south,
    
  )) %>% 
  mutate(conf.low2 = case_when(
    line == "gu" & group2 == -1 ~ conf.low - gu_int_north,
    line == "gu" & group2 == 0 ~ conf.low - gu_int_mid,
    line == "gu" & group2 == 1 ~ conf.low - gu_int_south,
    
    line == "fledge" & group2 == -1 ~ conf.low - fledge_int_north,
    line == "fledge" & group2 == 0 ~ conf.low - fledge_int_mid,
    line == "fledge" & group2 == 1 ~ conf.low - fledge_int_south,
    
    line == "leps" & group2 == -1 ~ conf.low - lep_int_north,
    line == "leps" & group2 == 0 ~ conf.low - lep_int_mid,
    line == "leps" & group2 == 1 ~ conf.low - lep_int_south
  )) %>% 
  mutate(conf.high2 = case_when(
    line == "gu" & group2 == -1 ~ conf.high - gu_int_north,
    line == "gu" & group2 == 0 ~ conf.high - gu_int_mid,
    line == "gu" & group2 == 1 ~ conf.high - gu_int_south,
    
    line == "fledge" & group2 == -1 ~ conf.high - fledge_int_north,
    line == "fledge" & group2 == 0 ~ conf.high - fledge_int_mid,
    line == "fledge" & group2 == 1 ~ conf.high - fledge_int_south,
    
    line == "leps" & group2 == -1 ~ conf.high - lep_int_north,
    line == "leps" & group2 == 0 ~ conf.high - lep_int_mid,
    line == "leps" & group2 == 1 ~ conf.high - lep_int_south
  )) 


# what are the sd for gdd.dev for greenup and lep dfs?
gdd_sd_gu <- sd(greenup$temp.dev)
gdd_mean_gu <- mean(greenup$temp.dev)
gdd_sd_leps <- sd(leps_gdd$temp.dev)
gdd_mean_leps <- mean(leps_gdd$temp.dev)
gdd_sd_fledge <- sd(fledge_gdd$temp.dev)
gdd_mean_fledge <- mean(fledge_gdd$temp.dev)

# what are the sd for mean_ffp for greenup and lep dfs?
ffp_sd_gu <- sd(greenup$mean_ffp)
ffp_mean_gu <- mean(greenup$mean_ffp)
ffp_sd_leps <- sd(leps_gdd$mean_ffp)
ffp_mean_leps <- mean(leps_gdd$mean_ffp)
ffp_sd_fledge <- sd(fledge_gdd$mean_ffp)
ffp_mean_fledge <- mean(fledge_gdd$mean_ffp)


head(tdf)

tdf2 <- tdf %>% 
  mutate(x2 = case_when(line == "gu" ~ (x * gdd_sd_gu) + gdd_mean_gu,
                        line == "leps" ~ (x * gdd_sd_leps) + gdd_mean_leps,
                        line == "fledge" ~ (x * gdd_sd_fledge) + gdd_mean_fledge)) %>% 
  mutate(
    group2 = case_when(
      group2 == 1 ~ "Southerly",
      group2 == -1 ~ "Northerly",
      group2 == 0 ~ "Central"),
    group2 = factor(group2, levels = c("Southerly", "Central", "Northerly"))
  )



ggplot(tdf2, mapping = aes(x = x2)) +
  geom_ribbon(mapping = aes(x = x2, ymin = conf.low2, ymax = conf.high2, 
                            fill = line), color = NA, alpha = 0.08) +
  geom_line(mapping = aes(x = x2, y=predicted2, color = line, linetype = line),  size = 1.2) +
  scale_linetype_manual(values = c(2,1,1), 
                        breaks = c("gu", "leps", "fledge"),
                        labels = c("Greenup (Plants)", "Emergence (Leps)", "Fledge (Birds)"),
                        guide = "none"
  ) +
  scale_color_manual(values = c("#316119","#1F8690", "#A73EBF"),
                     breaks = c("gu", "leps", "fledge"),
                     labels = c("Greenup (Plants)", "Emergence (Leps)", "Fledge (Birds)")) +
  scale_fill_manual(values = c("#316119","#1F8690", "#A73EBF"),
                    breaks = c("gu", "leps", "fledge"),
                    labels = c("Greenup (Plants)", "Emergence (Leps)", "Fledge (Birds)")) +
  labs(color = "", fill = "",
       x = "Temperature anomaly", y = "Phenology anomaly") +
  #scale_x_continuous(limits = c(-150,150)) +
  theme_bw() +
  guides(linetype = element_blank()) +
  facet_wrap(~group2)


ggsave(filename = "figures/tempDeviationEffects.png", width = 8, height = 3.5)

### get temporal residuals
cp <- cowplot::plot_grid(gu_resid_plot, leps_resid_plot, fledge_resid_plot, 
                         gu_acf_plot, leps_acf_plot, fledge_acf_plot)

cp

ggsave(filename = "figures/tempResiduals.png", plot = cp,
       width = 12, height = 8)

#R2 of the models
r.squaredGLMM(leps_tm)
# R2m       R2c
# [1,] 0.5591229 0.7328355
r.squaredGLMM(tm_fledge)
# R2m       R2c
# [1,] 0.1906912 0.5539464
r.squaredGLMM(gu_tm)
# R2m       R2c
# [1,] 0.7044038 0.9547725
