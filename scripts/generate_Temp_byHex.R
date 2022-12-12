library(terra)
library(sf)
library(dplyr)

# list of daymet files
years <- seq(2000,2020)

#function to make daymet averages by hex
# x = year
daymetByCell <- function(x){
  
  fp = file.path("g:/BigData/Daymet/monthly/")
  fp2 <- paste0(fp,"/", "daymet_v4_tmax_monavg_na_", x , ".tif")
  r <- rast(fp2)
  #select march, april, may, june
  r2 <- r[[3:6]]
  r2 <- terra::aggregate(r2,25)
  
  r3 <- mean(r2)
  plot(r3)
  
  
  # read in hex
  hex <- raster::shapefile("data/hex_grid_crop.shp") %>% 
    st_as_sf() 
  
  hex <- hex %>% 
    st_transform(crs = crs(r3))
  
  # do a spatial join after making r3 into a df
  r3_df <- raster::as.data.frame(r3, xy = T)
  memory.limit(size = 80000)
  r3_df_sf <- st_as_sf(r3_df, coords = c("x","y"), crs = crs(hex))
  
  r3_hex <- st_join(r3_df_sf, hex) %>% 
    dplyr::rename(temp = mean) %>% 
    dplyr::filter(!is.na(cell)) %>% 
    dplyr::filter(!is.na(temp))
  
  r3_df <- as.data.frame(st_drop_geometry(r3_hex)) %>% 
    group_by(cell) %>% 
    summarise(mean_temp = mean(temp, na.rm = T)) %>% 
    na.omit() %>% 
    mutate(Year = x)
  
  write.csv(x = r3_df,
            file = paste0("Outputs/daymetAverages/annualTemp", x, ".csv"))
}

lapply(X = years, FUN = daymetByCell)

# read all csvs into single csv
f <- list.files("Outputs/daymetAverages/", full.names = T)
tbl <- lapply(f, readr::read_csv) %>% 
  bind_rows()

tbl <- tbl %>% 
  dplyr::select(cell, mean_temp, Year) %>% 
  dplyr::rename(year = Year)

write.csv(x = tbl, file = "Outputs/temp_byHex.csv", row.names = F)
