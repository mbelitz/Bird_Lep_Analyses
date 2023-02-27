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

summary(gu_tm)
car::vif(gu_tm)

MuMIn::r.squaredGLMM(gu_tm)


plot_model(gu_tm, type = "pred", terms = "spring.dev")
plot_model(gu_tm, type = "pred", terms = c("mean_ffp"))
plot_model(gu_tm, type = "pred", terms = c("spring.dev", "mean_ffp"))

gu_gdd_plot <- plot_model(gu_tm, type = "pred", terms = c("spring.dev", "mean_ffp"))

gu_gdd_plot +
  theme_classic() +
  labs(x = "GDD anomoly", y = "Greenup") +
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

plot_model(tm_fledge, type = "pred", terms = c("spring.dev", "PC1")) +
  scale_color_manual(values = c("#C3DBC5","#2A628F","#13293D")) +
  scale_fill_manual(values = c("#C3DBC5","#2A628F","#13293D")) +
  labs(x = "GDD anomaly", y = "Fledge date", color = "Trait score", fill = "Trait score") +
  ggtitle("") +
  theme_classic()

ggsave(filename = "figures/fledge_sensitivity_byTraitScore.png",
       width = 4, height = 2.5)


fledge_gdd_plot <- plot_model(tm_fledge, type = "pred", terms = c("spring.dev", "mean_ffp"))

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
with(leps_scaled, table(code, year))
with(leps_scaled, table(code, spring.dev))

leps_m <- lmer(q5 ~ spring.dev + mean_ffp +
                 uniqObsDays + code +
                 spring.dev:code +
                 spring.dev:mean_ffp +
                 (1|cell), 
               data = leps_scaled, na.action = na.fail, REML = F)

leps_d <- dredge(leps_m)
leps_tm <- get.models(leps_d, 2)[[1]]

summary(leps_tm)
car::vif(leps_tm)
summary(leps_tm)
MuMIn::r.squaredGLMM(leps_tm)

plot_model(leps_tm, type = "pred", terms= "spring.dev")
plot_model(leps_tm, type = "pred", terms= c("mean_ffp"))

leps_gdd_plot <- plot_model(leps_tm, type = "pred", terms = c("spring.dev", "mean_ffp")) 

# get ACF for leps plot
resids_l <- residuals(leps_tm)
resid_leps <- leps_scaled %>% 
  mutate(resid = resids_l)

leps_resid_plot <- ggplot(resid_leps, mapping = aes(x = year, y = resid)) +
  geom_point(aes(color = code)) +
  geom_smooth(aes(color = code), method = "lm") +
  labs(x = "Year", y = "Residual")+
  theme_classic()

#no grouping
leps_resid_plot <- ggplot(resid_leps, mapping = aes(x = year, y = resid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Year", y = "Residual")+
  ggtitle('Lepidoptera') +
  theme_bw()

resid_leps_sort <- resid_leps[order(resid_leps$year),]

acfPlot <- acf(resid_leps_sort$resid) 

leps_acf_plot <- ggAcf(resid_leps_sort$resid) +
  theme_bw() +
  ggtitle("")

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
  mutate(group2 = rep(-1:1,22))

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
gdd_sd_gu <- sd(greenup$spring.dev)
gdd_mean_gu <- mean(greenup$spring.dev)
gdd_sd_leps <- sd(leps_gdd$spring.dev)
gdd_mean_leps <- mean(leps_gdd$spring.dev)
gdd_sd_fledge <- sd(fledge_gdd$spring.dev)
gdd_mean_fledge <- mean(fledge_gdd$spring.dev)

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
    group2 = factor(group2, levels = c("Northerly", "Central", "Southerly"))
  )


gdd_plot <- ggplot(tdf2, mapping = aes(x = x2)) +
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
       x = "GDD anomaly", y = "Phenology anomaly") +
  #scale_x_continuous(limits = c(-150,150)) +
  theme_classic() +
  guides(linetype = element_blank(),
         color=guide_legend(override.aes=list(fill=NA))) +
  theme(legend.position = "bottom") +
  facet_wrap(~group2, nrow = 3, ncol = 1)

gdd_plot

ggsave(filename = "figures/gddDeviationEffects.png", width = 8, height = 3.5)

### get temporal residuals
cp <- cowplot::plot_grid(gu_resid_plot, leps_resid_plot, fledge_resid_plot, 
                         gu_acf_plot, leps_acf_plot, fledge_acf_plot, 
                         labels = c("A", "B", "C",
                                    "D", "E", "F"))

cp

ggsave(filename = "figures/gddResiduals.png", plot = cp,
       width = 12, height = 8)


# source intro conceptual figure
source("scripts/generate_introConceptualFigure.R")

pp2 <- pp +
  ggtitle("Hypotheses") +
  theme_classic() +  
  theme(axis.ticks = element_blank()) 

gdd_plot2 <- gdd_plot +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "right")

pp2
gdd_plot2

ggsave(plot = pp2, filename = "figures/hypotheses.png", 
       width = 10,height = 3,dpi = 450)

ggsave(plot = gdd_plot2, filename = "figures/results2.png", 
       width = 5,height = 7,dpi = 450)

pp3 <- pp +
  theme_classic() +  
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) 

pp3

ggsave(plot = pp3, filename = "figures/hypotheses2.png", 
       width = 8,height = 4,dpi = 450)

# save to powerpoint dimensions
pp3 <- pp +
  ggtitle("Hypotheses") +
  theme_classic() +  
  theme(axis.ticks = element_blank(), 
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) 

pp3

ggsave(plot = pp3, filename = "figures/hypotheses3.png", 
       width = 12,height = 3,dpi = 450)


# combine all the plots together of gdd dev and gdd_plot

library(ggpubr)

cp <- ggarrange(pp2, gdd_plot2, nrow = 2, ncol = 1,
                common.legend = T, 
                label.x = -0.01, 
                label.y = 1.01, 
                widths = c(3,0.3),
                heights = c(2,3))
cp

ggsave(filename = "figures/Figure1.png", width = 8, height = 7)
