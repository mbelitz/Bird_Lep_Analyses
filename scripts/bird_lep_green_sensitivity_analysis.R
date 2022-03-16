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

ggsave(filename = "figures/clustering.png")

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

sp_gu_s <- step(sp_gu)

summary(sp_gu)


# deviation greenup model
t_gu <-  lmer(formula = gr_mn ~ spring.dev + year +
                  (1|cell), data = greenup)

t_gu_s <- step(t_gu)

summary(t_gu)

b <- plot_model(t_gu, type = "pred", terms = "spring.dev")
b

nd <- data.frame(spring.dev = -152.5:137.5, year = mean(greenup$year), cell = 1)

predict(t_gu,nd)

greenup_sens <- t_gu@beta[2]

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
sp_ar <- lmer(mean_arr ~ spring.gdd + year +  cluster + 
                spring.gdd:cluster +
                year:cluster +
                (1|cell), data = arr_gdd_pc_group)

sp_ar_s <- step(sp_ar)

car::vif(sp_ar)

sp_ar <- lmer(mean_arr ~ spring.gdd + year +  
                spring.gdd:cluster +
                year:cluster +
                (1|cell), data = arr_gdd_pc_group)

car::vif(sp_ar)

summary(sp_ar)

plot_model(sp_ar, type = "pred", terms = "spring.gdd")

# deviation model 
t_ar <- lmer(mean_arr ~ spring.dev + year + 
                spring.dev:cluster +
                year:cluster +
                (1|cell), data = arr_gdd_pc_group)

t_ar_s <- step(t_ar)
car::vif(t_ar)

summary(t_ar)


arr_spring.dev <-    t_ar@beta[2]
arr_spring.dev_c2 <- t_ar@beta[4]
arr_spring.dev_c3 <- t_ar@beta[5]

plot_model(t_ar, type = "pred", terms = c("year", "cluster"), ci.lvl = NA)
plot_model(t_ar, type = "pred", terms = c("spring.dev", "cluster"), ci.lvl = NA)

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
               spring.dev:cluster +
               year:cluster +
               (1|cell), data = fledge_gdd_pc_group)

t_fledge_s <- step(t_fledge)

summary(t_fledge)

fledge_spring.dev <- t_fledge@beta[2]
fledge_spring.dev_c2 <- t_fledge@beta[4]
fledge_spring.dev_c3 <- t_fledge@beta[5]

plot_model(t_fledge, type = "pred", terms = c("year", "cluster"), ci.lvl = NA)
plot_model(t_fledge, type = "pred", terms = c("spring.dev", "cluster"), ci.lvl = NA)


## lep sensitivity time
leps_gdd <- leps_gdd %>% 
  filter(!is.na(q5),
         !is.na(code),
         !is.na(spring.gdd),
         !is.na(spring.dev)) %>% 
  dplyr::select(code, cell, year, q5, uniqObsDays,
                spring.gdd, summer.gdd,
                spring.dev, summer.dev)


# spatiotemporal arrival model
sp_leps <- lmer(q5 ~ spring.gdd + year + spring.gdd:year +
                    code + 
                  uniqObsDays +
                    spring.gdd:code +
                    year:code +
                    spring.gdd:year:code +
                    (1|cell), data = leps_gdd)

summary(sp_leps)

plot_model(sp_leps, type = "pred", terms = c("year", "spring.gdd", "code"))

# deviation model 
t_leps <- lmer(q5 ~ spring.dev + year + 
                  uniqObsDays +
                   spring.dev:code +
                   year:code +
                   (1|cell), data = leps_gdd)

summary(t_leps)

lep_spring.dev <- t_leps@beta[2]
lep_spring.dev_RL <- t_leps@beta[4]
lep_spring.dev_RP <- t_leps@beta[5]


plot_model(t_leps, type = "pred", terms = c("year", "code"), ci.lvl = NA)
p <- plot_model(t_leps, type = "pred", terms = c("spring.dev", "code"), ci.lvl = NA) 


## plot together
library(ggpubr)

arr_sens <- ggplot() +
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

# get redids
lep_resid <- resid(t_leps)
head(t_leps)
t_leps <- resids <- t_leps %>%
mutate(resids = t_leps)
leps_resids <- leps_gdd %>%
mutate(resids = t_leps)
leps_resids <- leps_gdd %>%
mutate(resids = lep_resid)
## get resids
lep_resid <- resid(t_leps)
leps_resids <- leps_gdd %>%
mutate(resids = lep_r)
## get resids
lep_r <- resid(t_leps)
leps_resids <- leps_gdd %>%
mutate(resids = lep_r)
arr_r <- resids(t_arr)
arr_r <- resid(t_arr)
t_ar
arr_r <- resid(t_ar)
arr_resids <- arr_gdd_pc_group %>%
mutate(resids = arr_r)
arr_gdd_pc_group
arr_resids <- as.data.frame(arr_gdd_pc_group) %>%
mutate(resids = arr_r)
fledge_r <- resid(t_fledge)
fledge_resids <- as.data.frame(fledge_gdd_pc_group) %>%
mutate(resids = fledge_r)
green_r <- resid(t_gu)
green_resids <- as.data.frame(greenup) %>%
mutate(resids = green_r)
head(green_resids)

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
