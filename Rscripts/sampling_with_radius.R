# This script was used to sample predictor datasets at 0.1 deg at the recharge observation points

library(raster)

data <- read.csv("High_res_data_01/poc_01.csv", header=T)

coordinates(data) <- c("Long", "Lat")
proj4string(data) <- CRS("+proj=longlat +datum=WGS84")

data.coord <- as.data.frame(data@coords)

chirps <- raster("High_res_data_01/LTA_CHIRPS_clipped_01.tif")

#cell <- raster::extract(chirps, data.coord, method='bilinear')
cell <- raster::extract(chirps, data.coord, buffer=20000, fun='mean')

data$LTA_CHIRPS_mmpa <- cell

pet <- raster("High_res_data_01/LTA_PET_Afr_01.tif")

cell <- raster::extract(pet, data.coord, buffer=20000, fun='mean')

data$LTA_PET_mm <- cell

ai <- raster("High_res_data_01/LTA_AI_Afr_01_scaled.tif")

cell <- raster::extract(ai, data.coord, buffer=20000, fun='mean')
cell <- cell/10000
data$Aridity <- cell

ndvi <- raster("High_res_data_01/NDVI_Afr_01.tif")

cell <- raster::extract(ndvi, data.coord, buffer=20000, fun='mean')

data$NDVI <- cell

sm <- raster("High_res_data_01/Soil_moisture_Afr_01.tif")

cell <- raster::extract(sm, data.coord, buffer=20000, fun='mean')

data$SM10_m3m3 <- cell

data

write.csv(data,"poc_01_m.csv", row.names = FALSE)

