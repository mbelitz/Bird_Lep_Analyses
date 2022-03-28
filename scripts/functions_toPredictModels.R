generate_arrival_predicted_sf <- function(lowmidhigh, grouping){
  
  pdf <- data.frame(spring.dev = rep(mean(arr_gdd_pc_group$spring.dev) - 
                      (sd(arr_gdd_pc_group$spring.dev) * lowmidhigh),73),
                    cluster = as.character(rep(grouping,73)),
                    cell = unique(arr_gdd_pc_group$cell),
                    year = rep(mean(arr_gdd_pc_group$year),73))
  
  pred <- predictInterval(t_ar, pdf, which = "full")
  
  df <- cbind(pdf, pred)
  df$cell <- as.character(df$cell)
  
  sf_lj <- left_join(hex_sf, df) %>% 
    filter(!is.na(fit))
                
}
