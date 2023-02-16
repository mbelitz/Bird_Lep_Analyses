library(ncf)

source("scripts/plot_gddDeviationlEffects.R")

#spatial autocorrelation
# read in 
library(raster)
r <- raster::raster("data/bioclim/Normal_1991_2020_FFP.tif")
hex_sf <- raster::shapefile("data/hex_grid_crop.shp") %>% 
  st_as_sf()
hex_sf_ea <- st_transform(hex_sf, crs = crs(r))
hex_coords <- st_coordinates(st_centroid(hex_sf_ea))
hex_coords_df <- data.frame(cell = hex_sf_ea$cell, 
                            X = hex_coords[,1],
                            Y = hex_coords[,2])

# greenup
resid_greenup <- left_join(resid_greenup, hex_coords_df)

# fit correlog
ncf_fit <- correlog(x = resid_greenup$X, y = resid_greenup$Y, z = resid_greenup$resid,
                    increment = 250000, latlon = F, resamp = 100)
plot(ncf_fit)

ncf_gu_df <- data.frame(distance = ncf_fit$mean.of.class, 
                        correlation = ncf_fit$correlation,
                        p.value = if_else(ncf_fit$p < 0.025, 
                                          true = "Sig",
                                          false = "Not Sig")) %>% 
  filter(distance < 4000000 & distance > 1)
ncf_gu_df

gu_correlogram_plot <- ggplot() +
  geom_line(ncf_gu_df,
            mapping = aes(x = distance, y = correlation)) +
  geom_point(ncf_gu_df,
             mapping = aes(x = distance, y = correlation, fill = p.value),
             size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("white", "black")) +
  labs(x = "Distance", y = "Moran's I", fill = "P-Value") +
  theme_classic() +
  theme(legend.position = "none")

gu_correlogram_plot # no autocorrelation found at closest spatial lag, so spatial model not made


# fledge time
resid_fledge$cell <- as.character(resid_fledge$cell)
resid_fledge <- left_join(resid_fledge, hex_coords_df)

ncf_fit_fledge <- correlog(x = resid_fledge$X, y = resid_fledge$Y, z = resid_fledge$resid,
                    increment = 250000, latlon = F, resamp = 100)
plot(ncf_fit_fledge)

ncf_fledge_df <- data.frame(distance = ncf_fit_fledge$mean.of.class, 
                            correlation = ncf_fit_fledge$correlation,
                            p.value = if_else(ncf_fit_fledge$p < 0.025, 
                                              true = "Sig",
                                              false = "Not Sig")) %>% 
  filter(distance < 4000000 & distance > 1)


fledge_correlogram_plot <- ggplot() +
  geom_line(ncf_fledge_df,
            mapping = aes(x = distance, y = correlation)) +
  geom_point(ncf_fledge_df,
             mapping = aes(x = distance, y = correlation, fill = p.value),
             size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("white", "black")) +
  labs(x = "Distance", y = "Moran's I", fill = "P-Value") +
  theme_classic() +
  theme(legend.position = "none")

fledge_correlogram_plot # no autocorrelation found at closest spatial lag

# leps time
resid_leps$cell <- as.character(resid_leps$cell)
resid_leps <- left_join(resid_leps, hex_coords_df)

ncf_fit_leps <- correlog(x = resid_leps$X, y = resid_leps$Y, z = resid_leps$resid,
                    increment = 250000, latlon = F, resamp = 100)
plot(ncf_fit_leps)

ncf_leps_df <- data.frame(distance = ncf_fit_leps$mean.of.class, 
                          correlation = ncf_fit_leps$correlation,
                          p.value = if_else(ncf_fit_leps$p < 0.025, 
                                            true = "Sig",
                                            false = "Not Sig")) %>% 
  filter(distance < 4000000 & distance > 1)



leps_correlogram_plot <- ggplot() +
  geom_line(ncf_leps_df,
            mapping = aes(x = distance, y = correlation)) +
  geom_point(ncf_leps_df,
             mapping = aes(x = distance, y = correlation, fill = p.value),
             size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("white", "black")) +
  labs(x = "Distance", y = "Moran's I", fill = "P-Value") +
  theme_classic() +
  theme(legend.position = "bottom")

leps_correlogram_plot # no autocorrelation found at closest spatial lag, 


#autocorrelation plots
cp <-  cowplot::plot_grid(gu_resid_plot, leps_resid_plot, fledge_resid_plot,
                          gu_acf_plot, leps_acf_plot, fledge_acf_plot, 
                          gu_correlogram_plot, leps_correlogram_plot, fledge_correlogram_plot,
                          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                          ncol = 3, nrow = 3)
cp


ga <-  ggarrange(gu_resid_plot, leps_resid_plot, fledge_resid_plot,
                          gu_acf_plot, leps_acf_plot, fledge_acf_plot, 
                          gu_correlogram_plot, leps_correlogram_plot, fledge_correlogram_plot,
                          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                          ncol = 3, nrow = 3,
                 common.legend = T, legend = "bottom")


ggsave(filename = "figures/FigureS2.png", plot = ga, width = 7, height = 6)
