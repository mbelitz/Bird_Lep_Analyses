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
# read in frost free period data
ffp <- read.csv("Outputs/frostFreePeriod_byHex.csv")

# combine phenometrics with gdd & ffp
arr_gdd <- left_join(arr, gdd)
arr_gdd <- left_join(arr_gdd, ffp)
fledge_gdd <- left_join(fledge,gdd)
fledge_gdd <- left_join(fledge_gdd, ffp) %>% 
  filter(!is.na(mean_ffp))
leps_gdd <- left_join(leps, gdd)
leps_gdd <- left_join(leps_gdd, ffp)

# add PC based clusters to df
arr_gdd_pc <- left_join(arr_gdd, bird_pc)
#fledge_gdd_pc <- left_join(fledge_gdd, bird_pc, by = "species")

# look into sensitivty of phenometrics to gdd through space and time
# start with greenup
greenup <- arr_gdd_pc %>% 
  dplyr::distinct(cell, year, .keep_all = T) %>% 
  dplyr::select(cell, year, gr_mn, spring.gdd, summer.gdd, 
                spring.dev, summer.dev, mean_ffp)
greenup$cell <- as.character(greenup$cell)
# scale greenup dataframe
greenup_scaled <- greenup %>% 
  mutate(spring.dev = scale(spring.dev),
         year = scale(year),
         mean_ffp = scale(mean_ffp))

# deviation greenup model
gu_m <-  lmer(formula = gr_mn ~ spring.dev + year + mean_ffp +
                spring.dev:mean_ffp +
                (1|cell), data = greenup_scaled, na.action = na.fail)

gu_d <- dredge(gu_m)

gu_tm <- get.models(gu_d, subset = 1)[[1]]
gu_tm

summary(gu_tm)
car::vif(gu_tm)

MuMIn::r.squaredGLMM(gu_tm)

plot_model(gu_tm, type = "pred", terms = "spring.dev")
plot_model(gu_tm, type = "pred", terms = c("year"))
plot_model(gu_tm, type = "pred", terms = c("mean_ffp"))
plot_model(gu_tm, type = "pred", terms = c("spring.dev", "mean_ffp"))

gu_gdd_plot <- plot_model(gu_tm, type = "pred", terms = c("spring.dev"))

gu_gdd_plot +
  theme_classic() +
  labs(x = "GDD anomoly", y = "Greenup") +
  ggtitle("")

## now get this for fledge and leps
fledge_gdd <- fledge_gdd %>% 
  filter(!is.na(juv_meanday),
         !is.na(spring.gdd),
         !is.na(spring.dev)) %>% 
  dplyr::select(sci_name, station, year, cell, juv_meanday, 
                spring.gdd, summer.gdd,
                spring.dev, summer.dev, PC1, PC2, species, mean_ffp)

head(fledge_gdd)

fledge_scaled <- fledge_gdd %>% 
  mutate(spring.dev = scale(spring.dev),
         year = scale(year),
         mean_ffp = scale(mean_ffp)) %>% 
  filter(!is.na(juv_meanday),
         !is.na(spring.gdd),
         !is.na(spring.dev),
         !is.na(year)) 

# deviation arrival model
fledge_m <- lmer(juv_meanday ~ spring.dev + year + mean_ffp +
                   PC1 +
                   spring.dev:mean_ffp +
                   spring.dev:PC1 +
                   year:PC1 +
                   (1|station) + (1|sci_name),
                 data = fledge_scaled, na.action = na.fail)

#dredge
fledge_d <- dredge(fledge_m)
tm_fledge <- get.models(fledge_d, subset = 1)[[1]]

summary(tm_fledge)
MuMIn::r.squaredGLMM(tm_fledge)

plot_model(tm_fledge, type = "pred", terms = "spring.dev")
plot_model(tm_fledge, type = "pred", terms = c("PC1"))
plot_model(tm_fledge, type = "pred", terms = c("mean_ffp"))
plot_model(tm_fledge, type = "pred", terms = c("spring.dev", "mean_ffp"))

fledge_gdd_plot <- plot_model(tm_fledge, type = "pred", terms = c("spring.dev"))

fledge_gdd_plot +
  theme_classic() +
  labs(x = "GDD anomaly", y = "Greenup") +
  ggtitle("")

fp_int <- plot_model(tm_fledge, type = "pred", terms = c("spring.dev", "mean_ffp"))

fp_int +
  theme_classic() +
  labs(x = "GDD anomaly", y = "Fledge",
       color = "Frost free period", fill = "Frost free period") +
  ggtitle("")

ggsave("figures/ffp_spring.dev_greenup.png", width = 6, height = 3.5)


# there is an important interaction between mean_ffp and spring.dev
## now lep time
# deviation model 
leps_gdd <- leps_gdd %>% 
  filter(!is.na(spring.dev),
         !is.na(year),
         !is.na(mean_ffp),
         !is.na(year),
         !is.na(cell),
         !is.na(code),
         !is.na(q5))

leps_scaled <- leps_gdd %>% 
  mutate(spring.dev = scale(spring.dev),
         year = scale(year),
         mean_ffp = scale(mean_ffp),
         uniqObsDays = scale(uniqObsDays)) 

# note that we can't look at the interaction between code and year because there 
# is not enough RE data in early years
with(leps_scaled, table(code, year))
with(leps_scaled, table(code, spring.dev))

leps_m <- lmer(q5 ~ spring.dev + year + mean_ffp +
                 uniqObsDays + code +
                 spring.dev:code +
                 spring.dev:mean_ffp +
                 (1|cell), 
               data = leps_scaled, na.action = na.fail)

leps_d <- dredge(leps_m)
leps_tm <- get.models(leps_d, 3)[[1]]

summary(leps_tm)
car::vif(leps_tm)
summary(leps_tm)
MuMIn::r.squaredGLMM(leps_tm)

plot_model(leps_tm, type = "pred", terms= "spring.dev")
plot_model(leps_tm, type = "pred", terms= c("mean_ffp"))
plot_model(leps_tm, type = "pred", terms= "year")

leps_gdd_plot <- plot_model(leps_tm, type = "pred", terms = c("spring.dev")) 

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

lep_int <- dplyr::filter(tdf, x == 0, line == "leps")$predicted
fledge_int <- dplyr::filter(tdf, x == 0, line == "fledge")$predicted
gu_int <- dplyr::filter(tdf, x == 0, line == "gu")$predicted

tdf <- tdf %>% 
  mutate(predicted2 = case_when(
    line == "gu" ~ predicted - gu_int,
    line == "fledge" ~ predicted - fledge_int,
    line == "leps" ~ predicted - lep_int
  )) %>% 
  mutate(conf.low2 = case_when(
    line == "gu" ~ conf.low - gu_int,
    line == "fledge" ~ conf.low - fledge_int,
    line == "leps" ~ conf.low - lep_int
  )) %>% 
  mutate(conf.high2 = case_when(
    line == "gu" ~ conf.high - gu_int,
    line == "fledge" ~ conf.high - fledge_int,
    line == "leps" ~ conf.high - lep_int
  )) 


# what are the sd for gdd.dev for greenup and lep dfs?
gdd_sd_gu <- sd(greenup$spring.dev)
gdd_mean_gu <- mean(greenup$spring.dev)
gdd_sd_leps <- sd(leps_gdd$spring.dev)
gdd_mean_leps <- mean(leps_gdd$spring.dev)
gdd_sd_fledge <- sd(fledge_gdd$spring.dev)
gdd_mean_fledge <- mean(fledge_gdd$spring.dev)

head(tdf)

tdf2 <- tdf %>% 
  mutate(x2 = case_when(line == "gu" ~ (x * gdd_sd_gu) + gdd_mean_gu,
                        line == "leps" ~ (x * gdd_sd_leps) + gdd_mean_leps,
                        line == "fledge" ~ (x * gdd_sd_fledge) + gdd_mean_fledge))



ggplot() +
  geom_ribbon(filter(tdf2, line != "gu"),
              mapping = aes(x = x2, ymin = conf.low2, ymax = conf.high2, 
                            fill = interaction(line,group)), color = NA, alpha = 0.2) +
  geom_line(filter(tdf2, line != "gu"),
            mapping = aes(x = x2, y=predicted2, color = interaction(line,group)),
            size = 1.0) +
  geom_ribbon(filter(tdf2, line == "gu"),
              mapping = aes(x = x2, ymin = conf.low2, ymax = conf.high2, 
                            fill = interaction(line,group)), color = NA, alpha = 0.2) +
  geom_line(filter(tdf2, line == "gu"),
            mapping = aes(x = x2, y=predicted2, color = interaction(line,group)),
            linetype = "dashed", size = 1.2) +
  scale_color_manual(values = c("#2A628F","forest green", "#EE964B"),
                     labels = c("Fledge (Birds)", "Greenup (Plants)", "Emergence (Leps)")) +
  scale_fill_manual(values = c("#2A628F","forest green", "#EE964B"),
                    labels = c("Fledge (Birds)", "Greenup (Plants)", "Emergence (Leps)"))+
  labs(color = "", fill = "",
       x = "GDD anomaly", y = "Pheno-phase") +
  theme_bw()

ggsave(filename = "figures/gddDeviationEffects.png", width = 6, height = 3.5)
