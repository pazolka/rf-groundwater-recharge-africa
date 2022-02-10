######################################################################################################
# R codes to reproduce groundwater recharge interpolation map of Africa
######################################################################################################

library(nlme)
library(lme4)
library(sp)
library(gstat)
library(rgdal)
library(raster)

library(sf)
library(tidyverse)
library(giscoR)
library(maptools)
library(ggmap)
library(rasterize)

data <- read.csv("poc.csv", header=T)

# project coordinates to meters
coordinates(data) <- c("Long", "Lat")
proj4string(data) <- CRS("+proj=longlat +datum=WGS84")

data.coord <- as.data.frame(data@coords)

# Transform the proection to WGS 1984 World Mercator projection for African Continent
#data.prj <- spTransform(data, CRS("+init=epsg:3395 +units=km")) # World Mercator
# keep for now even though it's unnecessary
data.prj <- spTransform(data, CRS("+proj=longlat +datum=WGS84")) 
data.coord <- as.data.frame(data.prj@coords)

# log transformation of recharge values
data.prj$Recharge_mmpa <- log(data.prj$Recharge_mmpa) 

## fixed effects model
fmodel <- function(pcp) {-5 + 1.388*log(pcp)}

# residuals
res <- as.data.frame(data.prj$Recharge_mmpa - fmodel(data.prj$LTA_P_mmpa))
colnames(res) <- c("res")

# put residuals into a spatial data frame
dt <- SpatialPointsDataFrame(coordinates(data.prj), res)
#proj4string(dt) <- CRS("+init=epsg:3395 +units=km")
proj4string(dt) <-CRS("+proj=longlat +datum=WGS84")
dt <- dt[!duplicated(dt@coords),]              # remove any points with duplicated coordinates!

# empirical variogram
TheVariogram <- variogram(res~1, locations=coordinates(dt), data=dt)

# variogram model
vario.model <- vgm(psill=0.489, model="Mat", nugget=0.759, range=288, kappa=0.5)

# fit variogram
# FittedModel <- fit.variogram.reml(formula = res~1, locations=coordinates(dt), 
#                                   data=dt, model=vario.model)  

# !!!!!! not fitting the model anymore, taking optimal values from the paper
#FittedModel <- fit.variogram(TheVariogram, model=vario.model)
#vario.model <- fit.variogram(TheVariogram, model=vario.model)    
plot(TheVariogram, model=vario.model)

# get African grid
study_area <- readOGR("~/Desktop/rf-groundwater-recharge-africa/CRU_precip_CGIAR_AI_data/Africa_continent_shape.shp")
#study_area <- readOGR("African_continent.shp")

# read 0.5 precipitation raster
# BGS raster
precip <- raster("~/Desktop/rf-groundwater-recharge-africa/Recharge_files/Africa_bgs_LTA_AnnPrecip.tif")
#precip <- raster("Africa_bgs_LTA_AnnPrecip.tif")

proj4string(precip) <- CRS("+proj=longlat +datum=WGS84")

precip_Afr <- crop(precip, study_area)
plot(precip_Afr)
plot(study_area, add = TRUE)

# convert raster to data frame
precip_df <- as(precip_Afr, 'SpatialPointsDataFrame')
plot(precip_df)

#precip_df <- spTransform(precip_df, CRS("+init=epsg:3395 +units=km"))
plot(precip_df)

# calculate fixed effect for each pixel
fixed_afr <- fmodel(precip_df@data[["Africa_bgs_LTA_AnnPrecip"]])
is.na(fixed_afr) <- sapply(fixed_afr, is.infinite)                  # to replace -Inf with NA
fixed_afr[is.na(fixed_afr)] <- mean(fixed_afr, na.rm=T)             # replace NA with mean precip

# get kriging prediction for a residual

rand_afr <- gstat::krige(formula=res ~ 1, dt, newdata=precip_df, model=vario.model)

spplot(rand_afr["var1.pred"], main = "ordinary kriging predictions")
spplot(rand_afr["var1.var"],  main = "ordinary kriging variance")

# combine and back-transform
rand_afr@data$var1.pred <- exp(fixed_afr + rand_afr@data$var1.pred) 

spplot(rand_afr["var1.pred"],  main = "final recharge prediction")

recharge_df <- spTransform(rand_afr["var1.pred"], CRS("+proj=longlat +datum=WGS84"))
spplot(recharge_df)

final.gwr <- as.data.frame(recharge_df)

my.raster <- raster(xmn=-25.25, xmx=63.25, ymn=-34.75, ymx=37.25, resolution=c(0.5001,0.5001), 
                    crs="+proj=longlat +datum=WGS84")   # 0.5-deg does not work well!

final.raster <- rasterize(x=final.gwr[,2:3], y=my.raster, field=final.gwr[,1], fun=mean) 

rgb.pal <- colorRampPalette(c("light blue","blue","yellow","orange","red2"), space="rgb")   

plot(final.raster, col=rgb.pal(200))

writeRaster(final.raster, filename="BGS_Africa_gw_recharge_krig.tif", format="GTiff", overwrite=T)
