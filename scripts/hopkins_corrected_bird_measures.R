library(dplyr)
library(ggplot2)

## read in phenometrics
birds <- readRDS("data/phenoEstimates/pheno-data-2020-08-25.rds")

head(birds)
## remove bird dates where arr_GAM_mean is valid
birds <- birds %>% 
  dplyr::filter(!is.na(arr_GAM_mean))
head(birds)

## find cell with lowest latitude
lowest_cell <- birds %>% 
  filter(cell_lat == min(birds$cell_lat))

lowest_cell_lon <- lowest_cell$cell_lng[1]
lowest_cell_lat <- lowest_cell$cell_lat[1]

birds2 <- birds %>% 
  mutate(lon_cor = cell_lng - lowest_cell_lon,
         lat_cor = cell_lat - lowest_cell_lat)

# for every degree North, phenology delays 4 days
# for every degree east, phenology delays (4 days
# for every foot higher in elevation, phenology delays (1/100) days # but no elevation stuff today

mdf3 <- birds2 %>% 
  mutate(hopkins_cor =  (lat_cor * -4) + (lon_cor * -4))

ggplot() + 
  geom_point(mdf3, mapping = aes(x = cell_lng, y = cell_lat, color = hopkins_cor)) +
  scale_color_gradient2() +
  theme_classic()

# add hopkins correction to original onset estimates
mdf4 <- mdf3 %>% 
  dplyr::mutate(hopkins_arrival = arr_GAM_mean + hopkins_cor)

spp_onset <- mdf4 %>% 
  group_by(species) %>% 
  summarise(mean_onset = mean(arr_GAM_mean), sd_onset = sd(arr_GAM_mean))


## Can't use equinox to group species into traits. let's try quantiles
quantile(spp_onset$mean_onset, c(0.33,0.66))

ggplot() + 
  geom_histogram(spp_onset, mapping = aes(x = mean_onset), fill = "turquoise") +
  geom_vline(xintercept = c(114, 123)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() 

mdf5 <- left_join(mdf4, spp_onset)

## add in trait of early, mid, late
birds3 <- mdf5 %>% 
  mutate(season = case_when(mean_onset < 114 ~ "early",
                            mean_onset >= 114 & mean_onset <= 123 ~ "mid",
                            mean_onset > 123 ~ "late"))

write.csv(birds3, "Outputs/birds_hopkinsCorrected.csv", row.names = F)
