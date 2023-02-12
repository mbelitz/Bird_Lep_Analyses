source("scripts/plot_gddDeviationlEffects.R")

# predict green up models into low and high GDDs
cell <-  data.frame(cell = unique(greenup_scaled$cell))
# read in frost free period data
cell <- left_join(cell, greenup_scaled) %>% 
  dplyr::select(cell, mean_ffp) %>% 
  distinct(cell, mean_ffp)

## greenup
gu_pdf <- data.frame(spring.dev = rep(0, 73),
                     mean_ffp = cell$mean_ffp,
                     cell = cell$cell)

gu_pred <- predictInterval(gu_tm, newdata = gu_pdf, which = "full", n.sims = 1000)

gu_pred_df <- gu_pdf  %>% 
  cbind(gu_pred)

## leps
cell <-  data.frame(cell = unique(leps_scaled$cell))
cell <- left_join(cell, leps_scaled) %>% 
  dplyr::select(cell, mean_ffp) %>% 
  distinct(cell, mean_ffp)

leps_pdf <- data.frame(spring.dev = rep(0, 51),
                       mean_ffp = cell$mean_ffp,
                       cell = cell$cell,
                       code = rep("RL", 51),
                       uniqObsDays = rep(0, 51))

leps_pred <- predictInterval(leps_tm, newdata = leps_pdf, which = "full", n.sims = 1000)

leps_pred_df <- leps_pdf  %>% 
  cbind(leps_pred)

#Fledge
station <- distinct(ungroup(fledge_scaled), station, cell)
station$cell <- as.character(station$cell)
fledge_scaled$cell <- as.character(fledge_scaled$cell)
station <- left_join(station, fledge_scaled) %>% 
  dplyr::select(cell,station, mean_ffp) %>% 
  distinct(cell, station, mean_ffp)
sci_name <- "Vireo olivaceus"

fledge_pdf <- data.frame(spring.dev = rep(0, 256),
                         mean_ffp = station$mean_ffp,
                         station = station$station,
                         PC1 = rep(0.7328464, 256),
                         sci_name = rep(sci_name, 256))


fledge_pred <- predictInterval(tm_fledge, newdata = fledge_pdf, which = "full", n.sims = 1000)

fledge_pred_df <- fledge_pdf  %>% 
  cbind(fledge_pred)

## now make model predictions spatial
library(sf)
hex_sf <- raster::shapefile("data/hex_grid_crop.shp") %>% 
  st_as_sf()
hex_sf <- st_transform(hex_sf, crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

usca <- rnaturalearth::ne_countries(country = c("United States of America", "Canada"),
                                    returnclass = "sf", scale = 10) %>% 
  st_transform(crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

hbc <- raster::shapefile("data/hudsonBayCover.shp") %>% 
  st_as_sf()
hbc_sf <- st_transform(hbc, crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

gu_sf <- left_join(hex_sf, gu_pred_df) %>% 
  filter(!is.na(fit))

gu_plot <- ggplot() +
  geom_sf(usca, mapping = aes(), fill = "grey99") +
  geom_sf(hbc_sf, mapping = aes(), fill = "grey99", color = NA) +
  geom_sf(gu_sf, mapping = aes(fill = fit), alpha = 0.75) +
  scale_fill_viridis_c(option = "mako", limits = c(38, 220)) +
  coord_sf(xlim = c(-22260, 3003338),
           ylim = c(-1664985,1754536)) +
  ggtitle("Green-up") +
  labs(fill = "Phenophase") +
  theme_void()+
  theme(legend.position = "top", legend.title.align = 0.5) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

## lep time
leps_pred_df$cell <- as.character(leps_pred_df$cell)

leps_sf <- left_join(hex_sf, leps_pred_df) %>% 
  filter(!is.na(fit))

leps_plot <- ggplot() +
  geom_sf(usca, mapping = aes(), fill = "grey99") +
  geom_sf(hbc_sf, mapping = aes(), fill = "grey99", color = NA) +
  geom_sf(leps_sf, mapping = aes(fill = fit), alpha = 0.75) +
  scale_fill_viridis_c(option = "mako", limits = c(38, 220)) +
  coord_sf(xlim = c(-22260, 3003338),
           ylim = c(-1664985,1754536)) +
  ggtitle("Emergence") +
  labs(fill = "Phenophase") +
  theme_void() +
  theme(legend.position = "top", legend.title.align = 0.5) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

## fledge time
fledge_pred_df_U <- left_join(fledge_pred_df, distinct(fledge_scaled, station, .keep_all = T)
                              , by = "station")

fledge_pred_df_U$cell <- as.character(fledge_pred_df_U$cell)

fledge_sf <- left_join(hex_sf, fledge_pred_df_U) %>% 
  filter(!is.na(fit))

fledge_sf_gb <- fledge_sf %>% 
  group_by(cell) %>% 
  summarise(fit = mean(fit))

fledge_plot <- ggplot() +
  geom_sf(usca, mapping = aes(), fill = "grey99") +
  geom_sf(hbc_sf, mapping = aes(), fill = "grey99", color = NA) +
  geom_sf(fledge_sf_gb, mapping = aes(fill = fit), alpha = 0.75) +
  scale_fill_viridis_c(option = "mako", limits = c(38, 220)) +
  coord_sf(xlim = c(-22260, 3003338),
           ylim = c(-1664985,1754536)) +
  ggtitle("Fledge") +
  labs(fill = "Phenophase") +
  theme_void() +
  theme(legend.position = "top", legend.title.align = 0.5) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

fledge_plot

# common legend with ggarrange
library(ggpubr)

phenophaseLegend <- ggpubr::get_legend(gu_plot + guides(  fill = guide_colourbar(title.position = "left")) + theme(legend.position = "right", legend.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5), legend.key.width = unit(0.5, "cm")))
ga1 <- ggarrange(gu_plot, leps_plot, fledge_plot,
                 ncol = 1, nrow = 3,
                 common.legend = TRUE,
                 legend = "bottom")

ga1 <- ggarrange(
  gu_plot     + theme(legend.position = "none", title = element_text(size = 10)), 
  leps_plot   + theme(legend.position = "none", title = element_text(size = 10)),
  fledge_plot + theme(legend.position = "none", title = element_text(size = 10)),
                 ncol = 1, nrow = 3,
                 common.legend = F)




g2 <- ggarrange(ga1, gdd_plot2)

g2

source("scripts/generate_introConceptualFigure.R")

pp <- ggplot() +
  geom_smooth(sdf, mapping = aes(x = ga, y = pa, color = tl),method = "lm", se = F) +
  scale_color_manual(values = c("#316119","#1F8690", "#A73EBF"),
                     breaks = c("Greenup (Plants)", "Emergence (Leps)", "Fledge (Birds)")) +
  labs(color = "", 
       x = "GDD anomaly", y = "Phenology anomaly") +
  scale_x_continuous(breaks = c(-0.7, 0.7), labels = c("Cool", "Warm")) +
  scale_y_continuous(breaks = c(-0.9, 0.9), labels = c("Advanced", "Delayed")) +
  #scale_x_continuous(limits = c(-150,150)) +
  theme_bw() +
  theme(axis.ticks = element_blank()) +
  facet_wrap(~scenario)

pp

ga1 # vertical maps
gdd_plot2 # faceted results
pp # hypothesis

library(gridExtra)
bigPlot <- arrangeGrob(grobs = list(
  {pp + theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    plot.margin = margin(0.8,0,0,0.3, unit = "cm")
  )},
  {ga1 + theme(plot.margin = margin(0.8,0,0,0, unit = "cm"))}, 
  {gdd_plot2 + theme(
    legend.position = "none",
    panel.grid = element_blank(), 
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    plot.margin = margin(0.8,0,0,0.3, unit = "cm")
    )}
  ),
  
            layout_matrix = matrix(ncol = 2, nrow = 2, byrow = T, data = c(1,1,2,3)),
            widths = c(1,1.5),
            heights = c(1,3)
)

library(cowplot)
bigPlot2 <- cowplot::ggdraw() +
  cowplot::draw_grob(bigPlot) +
  draw_grob(phenophaseLegend, width = 0.1, height = 0.1, x = 0.29, y = 0.045) +
  draw_text("A - Hypotheses", hjust = 0, vjust = 1, x = 0.01, y =0.99) +
  draw_text("B - Spatial predictions", hjust = 0, vjust = 0, x = 0.01, y = 0.71) +
  draw_text("C - GDD sensitivity results", hjust = 0, vjust = 0, x = 0.4, y = 0.71)
  
ggsave(bigPlot2, filename = "figures/bigPlot.png", width = 6, height = 8, bg = "white")
