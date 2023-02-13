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

ggplot(fledge_scaled, mapping = aes(x = spring.dev, y = juv_meanday, color = sci_name)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  theme(legend.position = "none")

#dredge
fledge_d <- dredge(fledge_m)
tm_fledge <- get.models(fledge_d, subset = 1)[[1]]

summary(tm_fledge)
car::vif(tm_fledge)
MuMIn::r.squaredGLMM(tm_fledge)

plot_model(tm_fledge, type = "pred", terms = "spring.dev")
plot_model(tm_fledge, type = "pred", terms = c("PC1"))
plot_model(tm_fledge, type = "pred", terms = c("mean_ffp"))
plot_model(tm_fledge, type = "pred", terms = c("spring.dev", "mean_ffp"))
plot_model(tm_fledge, type = "pred", terms = c("spring.dev", "PC1"))

mylabs <- c("Short distance / early breeding", "Mid distance / average breeding", "Long distance / late breeding")
plot_model(tm_fledge, type = "pred", terms = c("spring.dev", "PC1 [-1,0,1]")) +
  scale_color_manual(values = c("#C3DBC5","#2A628F","#13293D"), breaks = c(-1,0,1), labels = mylabs) +
  scale_fill_manual(values = c("#C3DBC5","#2A628F","#13293D"), breaks = c(-1,0,1), labels = mylabs) +
  guides(
    color = guide_legend(reverse = T), fill = guide_legend(reverse = T) ) +
  labs(x = "GDD anomaly", y = "Fledge date", 
       color = "Trait score", fill = "Trait score") +
  ggtitle("") +
  theme_classic()

ggsave(filename = "figures/fledge_sensitivity_byTraitScore.png",
       width = 6, height = 3.5)
