library(dplyr)
library(ggplot2)
library(cowplot)

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

# green up data
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

# scale fledge data
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

# combine fledge_gu
fledge_scaled$cell <- as.character(fledge_scaled$cell)
fledge_gdd_gu <- left_join(fledge_scaled, greenup_scaled, by = c("year", "cell"))

a <- ggplot(fledge_gdd_gu, aes(x = gr_mn, y = juv_meanday, color = sci_name)) +
 # geom_point() +
  geom_smooth(method = "lm", se = F, size = 0.7) +
  scale_color_manual(values = rep("black", 43)) +
  theme_classic() +
  labs(x = "Green-up date", y = "Fledge date") +
  theme(legend.position = "none")

# bring leps in too
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

# combine fledge_gdd_gu
leps_scaled$cell <- as.character(leps_scaled$cell)
fledge_leps_gdd_gu <- left_join(fledge_gdd_gu, leps_scaled, by = c("year", "cell"))

fledge_leps_gdd_gu <- fledge_leps_gdd_gu %>% 
  mutate(OWS = case_when(code == "RE" ~ "Egg",
                         code == "RL" ~ "Larvae",
                         code == "RP" ~ "Pupae"))

b <- ggplot(fledge_leps_gdd_gu, aes(x = q5, y = juv_meanday, color = sci_name)) +
  # geom_point() +
  geom_smooth(method = "lm", se = F, size = 0.7) +
  scale_color_manual(values = rep("black", 43)) +
  theme_classic() +
  labs(x = "Emergence (Leps)", y = "Fledge date") +
  theme(legend.position = "none")

####### look at leps to to greenup sensitivity
c <- ggplot(fledge_leps_gdd_gu, aes(x = gr_mn, y = q5, color = OWS)) +
  # geom_point() +
  geom_smooth(method = "lm", se = F) +
  scale_color_manual(values = c("#575366", "#EDEEC0", "#433E0E")) +
  theme_classic() +
  labs(x = "Green-up date", y = "Emergence (Leps)", color = "Overwinter stage") +
  theme(legend.position = "bottom")

cp <- plot_grid(a,b,c, labels = c("A", "B", "C"), ncol = 1, nrow = 3)
cp  

ggsave(cp, filename = "figures/Figure3.png", width = 4, height = 8)
