library(dplyr)
library(ggplot2)

# create dataframe for plotting predictive conceptual fig
df <- data.frame(
  ga = c(-1,0,1,
         -1,0,1,
         -1,0,1), # gdd anomaly
  pa = c(1,0,-1,
         0.98, 0 , -0.98,
         1.02, 0 , -1.02), #phenology anomaly
  tl = c("Greenup (Plants)", "Greenup (Plants)","Greenup (Plants)",
         "Emergence (Leps)","Emergence (Leps)","Emergence (Leps)", 
         "Fledge (Birds)","Fledge (Birds)","Fledge (Birds)"),
  scenario = rep("Scenario A", 9)
)


# now plot with green up and leps togther, birds diff
df2 <- data.frame(
  ga = c(-1,0,1,
         -1,0,1,
         -1,0,1), # gdd anomaly
  pa = c(1,0,-1,
         0.98, 0 , -0.98,
         0.5, 0 , -0.5), #phenology anomaly
  tl = c("Greenup (Plants)", "Greenup (Plants)","Greenup (Plants)",
         "Emergence (Leps)","Emergence (Leps)","Emergence (Leps)", 
         "Fledge (Birds)","Fledge (Birds)","Fledge (Birds)"),
  scenario = rep("Scenario B", 9)
)


# now plot with all three having differing sensitivities 
df3 <- data.frame(
  ga = c(-1,0,1,
         -1,0,1,
         -1,0,1), # gdd anomaly
  pa = c(1,0,-1,
         0.75, 0 , -0.75,
         0.5, 0 , -0.5), #phenology anomaly
  tl = c("Greenup (Plants)", "Greenup (Plants)","Greenup (Plants)",
         "Emergence (Leps)","Emergence (Leps)","Emergence (Leps)", 
         "Fledge (Birds)","Fledge (Birds)","Fledge (Birds)"),
  scenario = rep("Scenario C", 9)
)


## combine dataframes together into super dataframe
sdf <- rbind(df, df2, df3)

pp <- ggplot() +
  geom_smooth(sdf, mapping = aes(x = ga, y = pa, color = tl),method = "lm") +
  scale_color_manual(values = c("#316119","#1F8690", "#A73EBF"),
                     breaks = c("Greenup (Plants)", "Emergence (Leps)", "Fledge (Birds)")) +
  labs(color = "", 
       x = "GDD anomaly", y = "Phenology anomaly") +
  scale_x_continuous(labels = c("Cold", "", "Average", "", "Warm")) +
  scale_y_continuous(labels = c("Advanced", "", "Average", "", "Delayed")) +
  #scale_x_continuous(limits = c(-150,150)) +
  theme_bw() +
  theme(axis.ticks = element_blank()) +
  facet_wrap(~scenario)

pp


