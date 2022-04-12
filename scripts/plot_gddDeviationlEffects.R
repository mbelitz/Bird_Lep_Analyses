library(dplyr)
library(ggplot2)
library(sf)
library(lme4)
library(lmerTest)
library(sjPlot)
library(merTools)
library(MuMIn)

# bird data
## read in phenometrics\
arr <- readRDS("data/phenoEstimates/pheno-data-2020-08-25.rds")
leps <- read.csv("data/phenoEstimates/adult_bfly_phenometrics_noCountCircles_withFull2020Data.csv") %>% 
  dplyr::rename(cell = HEXcell)
fledge <- readRDS("data/phenoEstimates/MAPS-fledge-dates-2022-02-22.rds") %>% 
  mutate(species = stringr::str_replace(sci_name, pattern = " ",
                                        replacement = "_"))

# read in bird pc values
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
  dplyr::select(cell, year, gr_mn, spring.gdd, summer.gdd, 
                spring.dev, summer.dev, FFD)
greenup$cell <- as.character(greenup$cell)
greenup <- na.omit(greenup)

# deviation greenup model
t_gu <-  lmer(formula = gr_mn ~ spring.dev + year + FFD +
                spring.dev:FFD +
                (1|cell), data = greenup, na.action = na.fail)

t_gu_d <- dredge(t_gu)

t_gu <- get.models(t_gu_d, subset = 1)[[1]]
t_gu

summary(t_gu)
car::vif(t_gu)

## top model does not have interaction betwen Spring GDD and FFD)
t_gu <- lmer(formula = gr_mn ~ spring.dev + year + FFD +
               (1|cell), data = greenup)

MuMIn::r.squaredGLMM(t_gu)

gu_gdd_plot <- plot_model(t_gu, type = "eff", terms = c("spring.dev"))

gu_gdd_plot +
  theme_classic() +
  labs(x = "GDD Deviation", y = "Greenup") +
  ggtitle("")

## now get this for fledge and leps
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


# deviation arrival model
t_fledge1 <- lmer(mean_fledge ~ spring.dev + year + FFD +
                    spring.dev:cluster +
                    (1|cell), data = fledge_gdd_pc_group)


summary(t_fledge1)
car::vif(t_fledge1)

## look into stepwise model for fledge
tm_fledge <-  t_fledge1

summary(tm_fledge)
MuMIn::r.squaredGLMM(tm_fledge)

fledge_gdd_plot <- plot_model(tm_fledge, type = "eff", terms = c("spring.dev", "cluster"))

fledge_gdd_plot +
  theme_classic() +
  ggtitle("")



## now lep time
# deviation model 
t_leps1 <- lmer(q5 ~ spring.dev + year + FFD +
                  uniqObsDays +
                  spring.dev:code +
                  (1|cell), data = leps_gdd)


summary(t_leps1)
car::vif(t_leps1)

leps_tm <- t_leps1

summary(leps_tm)
MuMIn::r.squaredGLMM(leps_tm)

leps_gdd_plot <- plot_model(t_leps1, type = "eff", terms = c("spring.dev", "code")) 

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

head(tdf)

re_int <- dplyr::filter(tdf, x == 0, group_col == "RE")$predicted
rl_int <- dplyr::filter(tdf, x == 0, group_col == "RL")$predicted
rp_int <- dplyr::filter(tdf, x == 0, group_col == "RP")$predicted

fledge1_int <- dplyr::filter(tdf, x == 0, group_col == "1", line == "fledge")$predicted
fledge2_int <- dplyr::filter(tdf, x == 0, group_col == "2", line == "fledge")$predicted
fledge3_int <- dplyr::filter(tdf, x == 0, group_col == "3", line == "fledge")$predicted

gu_int <- dplyr::filter(tdf, x == 0, line == "gu")$predicted

tdf <- tdf %>% 
  mutate(predicted2 = case_when(
    line == "gu" ~ predicted - gu_int,
    line == "fledge" & group_col == "1" ~ predicted - fledge1_int,
    line == "fledge" & group_col == "2" ~ predicted - fledge2_int,
    line == "fledge" & group_col == "3" ~ predicted - fledge3_int,
    group_col == "RE" ~ predicted - re_int,
    group_col == "RL" ~ predicted - rl_int,
    group_col == "RP" ~ predicted - rp_int
  )) %>% 
  mutate(conf.low2 = case_when(
    line == "gu" ~ conf.low - gu_int,
    line == "fledge" & group_col == "1" ~ conf.low - fledge1_int,
    line == "fledge" & group_col == "2" ~ conf.low - fledge2_int,
    line == "fledge" & group_col == "3" ~ conf.low - fledge3_int,
    group_col == "RE" ~ conf.low - re_int,
    group_col == "RL" ~ conf.low - rl_int,
    group_col == "RP" ~ conf.low - rp_int
  )) %>% 
  mutate(conf.high2 = case_when(
    line == "gu" ~ conf.high - gu_int,
    line == "fledge" & group_col == "1" ~ conf.high - fledge1_int,
    line == "fledge" & group_col == "2" ~ conf.high - fledge2_int,
    line == "fledge" & group_col == "3" ~ conf.high - fledge3_int,
    group_col == "RE" ~ conf.high - re_int,
    group_col == "RL" ~ conf.high - rl_int,
    group_col == "RP" ~ conf.high - rp_int
  )) 



ggplot() +
  geom_ribbon(filter(tdf, line != "gu"),
              mapping = aes(x = x, ymin = conf.low2, ymax = conf.high2, 
                            fill = interaction(line,group)), color = NA, alpha = 0.2) +
  geom_line(filter(tdf, line != "gu"),
            mapping = aes(x = x, y=predicted2, color = interaction(line,group)) ) +
  geom_ribbon(filter(tdf, line == "gu"),
              mapping = aes(x = x, ymin = conf.low2, ymax = conf.high2, 
                            fill = interaction(line,group)), color = NA, alpha = 0.2) +
  geom_line(filter(tdf, line == "gu"),
            mapping = aes(x = x, y=predicted2, color = interaction(line,group)),
            linetype = "dashed", size = 1.2) +
  scale_color_manual(values = c(
    rev(c("#3E92CC", "#2A628F", "#13293D")),"forest green",
    "#F95738", "#EE964B", "#F4D35E"
    
  )) +
  scale_fill_manual(values = c(
    rev(c("#3E92CC", "#2A628F", "#13293D")),"forest green",
    "#F95738", "#EE964B", "#F4D35E"
    
  )) +
  labs(color = "", fill = "",
       x = "GDD Deviation", y = "Pheno-phase") +
  theme_bw()

ggsave(filename = "figures/gddDeviationEffects.png")
