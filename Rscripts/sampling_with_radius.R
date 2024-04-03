# This script was used to sample predictor datasets at 0.1 deg at the recharge observation points

library(raster)

# parent directory is the working directory
data <- read.csv("Data/High_res_data_01/input_01.csv", header=T)

coordinates(data) <- c("Long", "Lat")
proj4string(data) <- CRS("+proj=longlat +datum=WGS84")

data.coord <- as.data.frame(data@coords)

chirps <- raster("Data/High_res_data_01/LTA_CHIRPS_Afr_01.tif")

#cell <- raster::extract(chirps, data.coord, method='bilinear')
cell <- raster::extract(chirps, data.coord, buffer=20000, fun='mean')

data$LTA_CHIRPS_mmpa <- cell

pet <- raster("Data/High_res_data_01/LTA_PET_Afr_01.tif")

cell <- raster::extract(pet, data.coord, buffer=20000, fun='mean')

data$LTA_PET_mm <- cell

ai <- raster("Data/High_res_data_01/LTA_AI_Afr_01_scaled.tif")

cell <- raster::extract(ai, data.coord, buffer=20000, fun='mean')
cell <- cell/10000
data$Aridity <- cell

ndvi <- raster("Data/High_res_data_01/LTA_NDVI_Afr_01.tif")

cell <- raster::extract(ndvi, data.coord, buffer=20000, fun='mean')

data$NDVI <- cell

sm <- raster("Data/High_res_data_01/LTA_soil_moisture_Afr_01.tif")

cell <- raster::extract(sm, data.coord, buffer=20000, fun='mean')

data$SM10_m3m3 <- cell

data

write.csv(data,"sampled_input_01.csv", row.names = FALSE)

