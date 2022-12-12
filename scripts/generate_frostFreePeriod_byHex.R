library(raster)
library(sf)
library(dplyr)

#  read in frost free period
r <- raster("data/bioclim/Normal_1991_2020_FFP.tif")
r

## aggregate r
r <- raster::aggregate(r, 25)
r
plot(r)

# read in 
hex_sf <- raster::shapefile("data/hex_grid_crop.shp") %>% 
  st_as_sf()

hex_sf_ea <- st_transform(hex_sf, crs = crs(r))
plot(hex_sf_ea)

# do a spatial join after making Frost free period into a df
ffp_df <- raster::as.data.frame(r, xy = T)
memory.limit(size=56000)
ffp_df_sf <- st_as_sf(ffp_df, coords = c("x","y"), crs = crs(hex_sf_ea))

ffp_hex <- st_join(ffp_df_sf, hex_sf_ea) %>% 
  filter(!is.na(cell)) %>% 
  filter(!is.na(Normal_1991_2020_FFP))

ffp_df <- as.data.frame(st_drop_geometry(ffp_hex)) %>% 
  group_by(cell) %>% 
  summarise(mean_ffp = mean(Normal_1991_2020_FFP, na.rm = T))
head(ffp_df)

write.csv(ffp_df, "Outputs/frostFreePeriod_byHex.csv", row.names = F)


## get residuals?

