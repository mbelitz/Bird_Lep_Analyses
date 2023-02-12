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
sci_name <- "Catharus fuscescens"

fledge_pdf <- data.frame(spring.dev = rep(0, 256),
                         mean_ffp = station$mean_ffp,
                         station = station$station,
                         PC1 = rep(1.923689, 256),
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

gu_plot2 <- ggplot() +
  theme_void()

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

leps_pred_df2 <- left_join(leps_pred_df, gu_pred_df, by = "cell") %>% 
  dplyr::mutate(diff = fit.x - fit.y)

leps_sf2 <- left_join(hex_sf, leps_pred_df2) %>% 
  filter(!is.na(diff))

leps_plot2 <- ggplot() +
  geom_sf(usca, mapping = aes(), fill = "grey99") +
  geom_sf(hbc_sf, mapping = aes(), fill = "grey99", color = NA) +
  geom_sf(leps_sf2, mapping = aes(fill = diff), alpha = 0.75) +
  scale_fill_gradient2() +
  coord_sf(xlim = c(-22260, 3003338),
           ylim = c(-1664985,1754536)) +
  labs(fill = "Emergence lag to green-up") +
  theme_void() +
  theme(legend.position = "bottom", legend.title.align = 0.5) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))


## fledge time
fledge_pred_df <- left_join(fledge_pred_df, fledge_scaled, by = "station")

fledge_pred_df$cell <- as.character(fledge_pred_df$cell)

fledge_sf <- left_join(hex_sf, fledge_pred_df) %>% 
  filter(!is.na(fit))

fledge_plot <- ggplot() +
  geom_sf(usca, mapping = aes(), fill = "grey99") +
  geom_sf(hbc_sf, mapping = aes(), fill = "grey99", color = NA) +
  geom_sf(fledge_sf, mapping = aes(fill = fit), alpha = 0.75) +
  scale_fill_viridis_c(option = "mako", limits = c(38, 220)) +
  coord_sf(xlim = c(-22260, 3003338),
           ylim = c(-1664985,1754536)) +
  ggtitle("Fledge") +
  labs(fill = "Phenophase") +
  theme_void() +
  theme(legend.position = "top", legend.title.align = 0.5) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

fledge_plot

fledge_pred_df2 <- left_join(fledge_pred_df, leps_pred_df, by = "cell") %>% 
  dplyr::mutate(diff = fit.x - fit.y)

fledge_sf2 <- left_join(hex_sf, fledge_pred_df2) %>% 
  filter(!is.na(diff))

fledge_plot2 <- ggplot() +
  geom_sf(usca, mapping = aes(), fill = "grey99") +
  geom_sf(hbc_sf, mapping = aes(), fill = "grey99", color = NA) +
  geom_sf(fledge_sf2, mapping = aes(fill = diff), alpha = 0.75) +
  scale_fill_gradient2() +
  coord_sf(xlim = c(-22260, 3003338),
           ylim = c(-1664985,1754536)) +
  labs(fill = "Fledge lag to emergence") +
  theme_void() +
  theme(legend.position = "bottom", legend.title.align = 0.5) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

fledge_plot2

# common legend with ggarrange
library(ggpubr)

phenophaseLegend <- ggpubr::get_legend(gu_plot)
ga1 <- ggarrange(gu_plot, leps_plot, fledge_plot,
                 ncol = 3,
                 common.legend = TRUE, legend = "none" )

ga1


library(gridExtra)
ga2 <- grid.arrange(phenophaseLegend,ga1, leps_plot2, fledge_plot2,
                    layout_matrix = matrix(c(1,1,1,2,2,2,NA,3,4), nrow = 3, ncol = 3, byrow = T),
                    heights = c(0.2,1,1)
)

plot(ga2)

ggsave(ga2, filename = "figures/spatialPrediction_withLag_Catharus_fuscescens_RL.png",
       width = 12, height = 8)
