
######################################################################################################
# R codes to reproduce groundwater recharge interpolation map of Africa
######################################################################################################

library(sp)
library(gstat)
library(rgdal)
library(raster)
library(spaMM)
library(nlme)
library(sf)
library(tidyverse)
library(giscoR)
library(maptools)
library(ggmap)
library(rasterize)

data <- read.csv("High_res_data_01/poc_01.csv", header=T)

coordinates(data) <- c("Long", "Lat")
proj4string(data) <- CRS("+proj=longlat +datum=WGS84")

data.coord <- as.data.frame(data@coords)

## fixed effects model
rain <- log(data$LTA_CHIRPS_mmpa)
recharge <- log(data$Recharge_mmpa)

new.data <- data.frame(recharge=recharge, rain=rain, x=data$Lat, y=data$Long)

# head(new.data)
mod.lmm <- fitme(recharge ~ rain + Matern(1 | x + y), data=new.data, family="gaussian")

# Model summary
summary(mod.lmm)

pred <- as.numeric(predict(mod.lmm, new.data, re.form=NA))  # re.form = NA used to remove spatial effects

# residuals 
res <- as.data.frame(log(data$Recharge_mmpa) - pred)

colnames(res) <- c("res")

# put residuals into a spatial data frame
dt <- SpatialPointsDataFrame(coordinates(data), res)
#proj4string(dt) <- CRS("+init=epsg:3395 +units=km")
proj4string(dt) <-CRS("+proj=longlat +datum=WGS84")
dt <- dt[!duplicated(dt@coords),]              # remove any points with duplicated coordinates!

# empirical variogram
TheVariogram <- variogram(res~1, locations=coordinates(dt), data=dt)

# variogram model
vario.model <- vgm(psill=1.184630, model="Mat", nugget=0.526024, range=603.6305, kappa=0.5)
# fit variogram
# FittedModel <- fit.variogram.reml(formula = res~1, locations=coordinates(dt), 
#                                   data=dt, model=vario.mo\del)  
vario.model <- fit.variogram(TheVariogram, model=vario.model)    
plot(TheVariogram, model=vario.model)

new.data['res'] <- res

##### reestimate fixed effect parameters with GLS
# Matern model with kappa=0.5 is an Exp model
mod.lmm.gls<-gls(recharge~rain, data=new.data, correlation=corExp(form=~res))

pred.gls <- as.numeric(predict(mod.lmm.gls, new.data, re.form=NA))  # re.form = NA used to remove spatial effects
# repeat -> recalculate residuals and create variogram
res.gls <- as.data.frame(log(data$Recharge_mmpa) - pred.gls)
colnames(res.gls) <- c("res")
dt.gls <- SpatialPointsDataFrame(coordinates(data), res.gls)
proj4string(dt.gls) <-CRS("+proj=longlat +datum=WGS84")
dt.gls <- dt.gls[!duplicated(dt.gls@coords),]
TheVariogramGLS <- variogram(res~1, locations=coordinates(dt.gls), data=dt.gls)
vario.model.gls <- vgm(psill=1.184630, model="Mat", nugget=0.526024, range=603.6305, kappa=0.5)
vario.model.gls <- fit.variogram(TheVariogramGLS, model=vario.model.gls)    
plot(TheVariogramGLS, model=vario.model.gls)

# get African grid
study_area <- readOGR("CRU_precip_CGIAR_AI_data/Africa_continent_shape.shp")
#study_area <- readOGR("African_continent.shp")

# read 0.1 precipitation raster
# CHIRPS raster
precip <- raster("~/Desktop/rf-groundwater-recharge-africa/High_res_data_01/LTA_CHIRPS_clipped_01.tif")
#precip <- raster("Africa_bgs_LTA_AnnPrecip.tif")

proj4string(precip) <- CRS("+proj=longlat +datum=WGS84")

precip_Afr <- crop(precip, study_area)
plot(precip_Afr)
plot(study_area, add = TRUE)

# convert raster to data frame
precip_df <- as(precip_Afr, 'SpatialPointsDataFrame')
plot(precip_df)

# calculate fixed effect for each pixel
#fixed_afr <- fmodel(precip_df@data[["LTA_CHIRPS_clipped_01"]])
new.data <- data.frame(rain=log(precip_df@data[["LTA_CHIRPS_clipped_01"]]), x=precip_df$x, y=precip_df$y)
fixed_afr <- as.numeric(predict(mod.lmm.gls, new.data, re.form=NA))

is.na(fixed_afr) <- sapply(fixed_afr, is.infinite)                  # to replace -Inf with NA
#fixed_afr[is.na(fixed_afr)] <- mean(fixed_afr, na.rm=T)             # replace NA with mean precip
fixed_afr[is.na(fixed_afr)] <- 0             # replace NA with 0 precip: Sahara

# get kriging predictions for residuals

rand_afr <- gstat::krige(formula=res ~ 1, dt.gls, newdata=precip_df, model=vario.model.gls)

spplot(rand_afr["var1.pred"], main = "ordinary kriging predictions")
spplot(rand_afr["var1.var"],  main = "ordinary kriging variance")

# combine and back-transform
rand_afr@data$var1.pred <- exp(fixed_afr + rand_afr@data$var1.pred) 

spplot(rand_afr["var1.pred"],  main = "final recharge prediction")

recharge_df <- spTransform(rand_afr["var1.pred"], CRS("+proj=longlat +datum=WGS84"))

spplot(recharge_df)

final.gwr <- as.data.frame(recharge_df)

my.raster <- raster(xmn=-25.25, xmx=63.25, ymn=-34.75, ymx=37.25, resolution=c(0.1001,0.1001), 
                    crs="+proj=longlat +datum=WGS84")

final.raster <- rasterize(x=final.gwr[,2:3], y=my.raster, field=final.gwr[,1], fun=mean) 

rgb.pal <- colorRampPalette(c("light blue","blue","yellow","orange","red2"), space="rgb")   

plot(final.raster, col=rgb.pal(200))

writeRaster(final.raster, filename="Rscripts/LMM_recharge_01.tif", format="GTiff", overwrite=T)
