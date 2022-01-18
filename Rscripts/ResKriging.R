library(nlme)
library(lme4)
library(sp)
library(gstat)
library(rgdal)
library(raster)

data <- read.csv("~/Desktop/rf-groundwater-recharge-africa/poc.csv")
# project coordinates to meters
coordinates(data) <- c("Long", "Lat")
proj4string(data) <- CRS("+proj=longlat +datum=WGS84")
# custom Chamberlin trimetric projected CRS
# https://gis.stackexchange.com/questions/364268/choosing-appropriate-projected-coordinate-system-for-africa
data <- spTransform(data, CRS("+proj=chamb +lat_1=22 +lon_1=0 +lat_2=22 +lon_2=45 +lat_3=-22 +lon_3=22.5 +datum=WGS84 +type=crs +units=km"))

# log transformation
data$Recharge_mmpa <- log(1+data$Recharge_mmpa)

## fixed effects model
fmodel <- function(pcp) {-5 + 1.388*log(pcp)}

# residuals
res = as.data.frame(data$Recharge_mmpa - fmodel(data$LTA_P_mmpa))
colnames(res) <- c('res')

# put residuals into a spatial data frame
dt <- SpatialPointsDataFrame(coordinates(data), res)
proj4string(dt) <- CRS("+proj=chamb +lat_1=22 +lon_1=0 +lat_2=22 +lon_2=45 +lat_3=-22 +lon_3=22.5 +datum=WGS84 +type=crs +units=km")

# empirical variogram
TheVariogram=variogram(res~1, locations=coordinates(dt), data=dt)
plot(TheVariogram)

# variogram model
TheVariogramModel <- vgm(psill=0.789, model="Mat", nugget=0.459, range=1200, kappa=0.5)
plot(TheVariogram, model=TheVariogramModel) 

# fit variogram
FittedModel <- fit.variogram.reml(formula = res~1, locations=coordinates(dt), 
                                  data=dt, model=TheVariogramModel)  
#FittedModel <- fit.variogram(TheVariogram, model=TheVariogramModel)    
plot(TheVariogram, model=FittedModel)

# get African grid
study_area <- readOGR("~/Desktop/rf-groundwater-recharge-africa/CRU_precip_CGIAR_AI_data/Africa_continent_shape.shp")

# read 0.5 precipitation raster
# BGS raster
precip <- raster("~/Desktop/rf-groundwater-recharge-africa/Recharge_files/Africa_bgs_LTA_AnnPrecip.tif")
precip_Afr <- crop(precip, study_area)
plot(precip_Afr)
plot(study_area, add = TRUE)

# convert raster to data frame
precip_df <- as(precip_Afr, 'SpatialPointsDataFrame')
plot(precip_df)

# reproject this data frame (reprojecting raster didnt work)
precip_df <- spTransform(precip_df, CRS("+proj=chamb +lat_1=22 +lon_1=0 +lat_2=22 +lon_2=45 +lat_3=-22 +lon_3=22.5 +datum=WGS84 +type=crs +units=km"))
plot(precip_df)

# calculate fixed effect for each pixel
fixed_afr <- fmodel(precip_df@data[["Africa_bgs_LTA_AnnPrecip"]])

# get kriging prediction for a residual
rand_afr <- krige(formula=res ~ 1, dt, precip_df, model=FittedModel)

spplot(rand_afr["var1.pred"], main = "ordinary kriging predictions")
spplot(rand_afr["var1.var"],  main = "ordinary kriging variance")

# combine and back-transform
rand_afr@data$var1.pred <- exp(fixed_afr + rand_afr@data$var1.pred)

spplot(rand_afr["var1.pred"],  main = "final recharge prediction")

# reproject again! doesnt work
recharge_df <- spTransform(rand_afr, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
