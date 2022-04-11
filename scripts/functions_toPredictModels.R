library(merTools)
library(dplyr)


predict_greenup_allGDD <- function(cellz){
  
  pdf <- data.frame(spring.dev = -100:100,
                    cell = rep(cellz, 201),
                    year = rep(mean(greenup$year),201),
                    FFD = rep(mean(greenup$FFD,201)))
  
  pred <- predictInterval(t_gu, pdf, which = "full",
                          n.sims = 10000, level = 0.8)
  
  df <- cbind(pdf, pred)
  df$cell <- as.character(df$cell)
  return(df)         
}

generate_arrival_predicted_sf <- function(lowmidhigh, grouping){
  
  MFFD <- arr_gdd_pc_group %>% 
    group_by(cell) %>% 
    summarise(meanFFD = mean(FFD))
  
  pdf <- data.frame(spring.dev = rep(mean(greenup$spring.dev) - 
                      (sd(greenup$spring.dev) * lowmidhigh),73),
                    cluster = as.character(rep(grouping,73)),
                    cell = MFFD$cell,
                    year = rep(mean(arr_gdd_pc_group$year),73),
                    FFD = MFFD$meanFFD)
  
  pred <- predictInterval(t_ar1, pdf, which = "full")
  
  df <- cbind(pdf, pred)
  df$cell <- as.character(df$cell)
  return(df)         
}

predict_arrival_allGDD <- function(grouping){
  
  pdf <- data.frame(spring.dev = -100:100,
                    cluster = as.character(rep(grouping,201)),
                    cell = rep("725", 201),
                    year = rep(mean(arr_gdd_pc_group$year),201),
                    FFD = rep(mean(arr_gdd_pc_group$FFD,201)))
  
  pred <- predictInterval(t_ar1, pdf, which = "full",
                          n.sims = 10000, level = 0.8)
  
  df <- cbind(pdf, pred)
  df$cell <- as.character(df$cell)
  return(df)         
}

generate_fledge_predicted_sf <- function(lowmidhigh, grouping){
  
  MFFD <- fledge_gdd_pc_group %>% 
    group_by(cell) %>% 
    summarise(meanFFD = mean(FFD))
  
  pdf <- data.frame(spring.dev = rep(mean(greenup$spring.dev) - 
                                       (sd(greenup$spring.dev) * lowmidhigh),68),
                    cluster = as.character(rep(grouping,68)),
                    cell = MFFD$cell,
                    year = rep(mean(fledge_gdd_pc_group$year),68),
                    FFD = MFFD$meanFFD)
  
  pred <- predictInterval(tm_fledge, pdf, which = "full")
  
  df <- cbind(pdf, pred)
  df$cell <- as.character(df$cell)
  return(df)
}

predict_fledge_allGDD <- function(grouping){
  
  pdf <- data.frame(spring.dev = -100:100,
                    cluster = as.character(rep(grouping,201)),
                    cell = rep("725", 201),
                    year = rep(mean(fledge_gdd_pc_group$year),201),
                    FFD = rep(mean(fledge_gdd_pc_group$FFD,201)))
  
  pred <- predictInterval(t_fledge1, pdf, which = "full",
                          n.sims = 10000, level = 0.8)
  
  df <- cbind(pdf, pred)
  df$cell <- as.character(df$cell)
  return(df)         
}

generate_lep_predicted_sf <- function(lowmidhigh, grouping){
  
  MFFD <- leps_gdd %>% 
    group_by(cell) %>% 
    summarise(meanFFD = mean(FFD, na.rm = T))
  
  pdf <- data.frame(spring.dev = rep(mean(greenup$spring.dev) - 
                                       (sd(greenup$spring.dev)*lowmidhigh), 54),
                    year = rep(mean(leps_gdd$year),54),
                    uniqObsDays = rep(mean(leps_gdd$uniqObsDays),54),
                    code = rep(grouping, 54),
                    cell = MFFD$cell,
                    FFD = MFFD$meanFFD)
  
  pred <- predictInterval(leps_tm, pdf, which = "full")
  
  df <- cbind(pdf, pred)
  df$cell <- as.character(df$cell)
  return(df)
}

predict_lep_allGDD <- function(grouping){
  
  pdf <- data.frame(spring.dev = -100:100,
                    code = as.character(rep(grouping,201)),
                    uniqObsDays = rep(mean(leps_gdd$uniqObsDays),201),
                    cell = rep("725", 201),
                    year = rep(mean(leps_gdd$year),201),
                    FFD = rep(mean(leps_gdd$FFD,201)))
  
  pred <- predictInterval(t_leps1, pdf, which = "full",
                          n.sims = 10000, level = 0.8)
  
  df <- cbind(pdf, pred)
  df$cell <- as.character(df$cell)
  return(df)         
}


