library(ggplot2)
library(tidyterra)
library(cowplot)
library(viridis)
library(ggspatial)
library(stringr)
library(sf)

# Total Richness Map
# Load and plot map
richness_map <- rast("Richness_map.tif")

total_BIN_richness <- ggplot() +
  geom_spatraster(data=richness_map, aes(fill=sum))+
  theme_cowplot()+
  scale_fill_viridis_c(name="BIN \nRichness", na.value = NA)+
  theme(axis.line= element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank())+
  annotation_scale()+
  annotation_north_arrow(location= "tl", height=unit(1, "cm"), width=unit(1, "cm"))

ggsave(total_BIN_richness, file= "total_BIN_richness.pdf")



# Order Richness Maps
California <- tigris::states() %>%
  filter(STUSPS == "CA")

Order_map <- list()

Order_files <- list.files("Order_maps/", full.names = TRUE)

for(o in Order_files){
  Order_name <- str_remove(str_remove(o, "Order_maps//"), ".tif")
  
  Order_raster <- rast(o)
  
  Order_map[[Order_name]] <- ggplot() +
    geom_sf(data= California, fill= "transparent")+
    geom_spatraster(data=Order_raster, aes(fill=sum))+
    theme_cowplot()+
    scale_fill_viridis_c(name="BIN \nRichness", na.value = NA)+
    theme(axis.line= element_blank(), 
          axis.ticks = element_blank(),
          axis.text = element_blank())+
    ggtitle(Order_name)
  
}

Order_map$Acerentomata

for(i in seq(1,29,4))
  ggsave(paste0("Temp",i,".pdf"),
         plot_grid(Order_map[[i]],
                   Order_map[[i+1]],
                   Order_map[[i+2]],
                   Order_map[[i+3]],
                   nrow = 2),
         width = 8, height = 11, units = "in")


ggsave(paste0("Temp",29,".pdf"),
       plot_grid(Order_map[[29]],
                 nrow = 2, ncol = 2),
       width = 8, height = 11, units = "in")


