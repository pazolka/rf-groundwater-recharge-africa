######################################################################################################
# R codes to add residual kriging to predicted LTA groundwater recharge (Random Forest in Python)
# authors: Anna Pazola, Mohammad Shamsudduha
######################################################################################################

library(gstat); library(sp); library(rgdal); library(RColorBrewer); library(ggplot2);
library(raster); library(rasterize)

# parent directory is the working directory
af.shape <- readOGR("Data/Africa_continent_shape.shp")

# read residuals file - calculated in python
res <- read.csv("Data/residuals.csv", header=T)

# cast residuals to numeric values
res$Res_05 <- sapply(res$Res_05, as.numeric)
res$Res_01 <- sapply(res$Res_01, as.numeric)

######## residuals at 0.5 deg
res5 <- res[,c(1,2,3)]
names(res5) <- c("X", "Y", "RES05")

# Setting existing coordinate as lat-long system
coordinates(res5) <- ~X+Y
proj4string(res5) <- CRS("+proj=longlat +datum=WGS84")
head(res5)
plot(res5)

# Senegal - two points for the same location with nearly the same values - take only one
res5 <- res5[-zerodist(res5)[,1],]
# remove na (on the edges of the continent)
sel <- !is.infinite(res5$RES05) & !is.na(res5$RES05)
res5 <- res5[sel,]

# variogram
vario5 <- variogram(RES05 ~ 1, res5[sel,], width=100, cutoff=950)
vario.mod5 <- vgm(psill=0.05, model="Sph", nugget=0.05, range=200)      # model variogram
vario.mod5 <- fit.variogram(vario5, model=vario.mod5)       # fit model variogram
plot(vario5, model=vario.mod5)
print(vario.mod5)
attr(vario.mod5, "SSErr")

############ residuals at 0.1 deg

res1 <- res[,c(1,2,4)]
names(res1) <- c("X", "Y", "RES01")

# Setting existing coordinate as lat-long system
coordinates(res1) <- ~X+Y
proj4string(res1) <- CRS("+proj=longlat +datum=WGS84")
head(res1)
plot(res1)

# Senegal - two points for the same location with nearly the same values - take only one
res1 <- res1[-zerodist(res1)[,1],]
# remove na (on the edges of the continent)
sel <- !is.infinite(res1$RES01) & !is.na(res1$RES01)
res1 <- res1[sel,]

# variogram
vario1 <- variogram(RES01 ~ 1, res1[sel,], width=75, cutoff=900)
vario.mod1 <- vgm(psill=0.04, model="Sph", nugget=0.04, range=200)    # model variogram
vario.mod1 <- fit.variogram(vario1, model=vario.mod1, fit.kappa=T)       # fit model variogram
plot(vario1, model=vario.mod1)   
print(vario.mod1)
attr(vario.mod1, "SSErr")

res.df <- as.data.frame(res1)

###### # make grid 0.5 deg - alternative to minimize the information loss when reprojecting
rf5.recharge <- raster("Data/Low_res_data_05/RF_recharge.tif")
rf5.recharge.afr <- crop(rf5.recharge, af.shape)
rf5.pixels <- as(rf5.recharge, 'SpatialPixels')
plot(rf5.pixels)

####### krige at 0.5 deg
res5.ok <- gstat::krige(RES05 ~ 1, res5, newdata=rf5.pixels, vario.mod5)

###### # make grid 0.1 deg - alternative to minimize the information loss when reprojecting
rf1.recharge <- raster("Data/High_res_data_01/RF_recharge_01.tif")
rf1.recharge.afr <- crop(rf1.recharge, af.shape)
rf1.pixels <- as(rf1.recharge, 'SpatialPixels')
plot(rf1.pixels)

####### krige at 0.1 deg
res1.ok <- gstat::krige(RES01 ~ 1, res1, newdata=rf1.pixels, vario.mod1)

# transform to decimal degrees and write rasters
res5.ok.ras <- raster(res5.ok)
res1.ok.ras <- raster(res1.ok)

res5.ok.ras <- projectRaster(res5.ok.ras, crs="+proj=longlat +datum=WGS84")
res1.ok.ras <- projectRaster(res1.ok.ras, crs="+proj=longlat +datum=WGS84")

####### combine RF predictions (python) and residual kriging (R)

# resample bc resolutions and extents are slightly off
res5.ok.ras.resampled = resample(res5.ok.ras, rf5.recharge, method='bilinear')
res1.ok.ras.resampled = resample(res1.ok.ras, rf1.recharge, method='bilinear')

plot(res5.ok.ras)
plot(res5.ok.ras.resampled)

plot(res1.ok.ras)
plot(res1.ok.ras.resampled)

# combine RF prediction (python) with residual kriging
rkrf5 <- 10**sum(log10(rf5.recharge), res5.ok.ras.resampled)
proj4string(rkrf5) <- CRS("+proj=longlat +datum=WGS84")
rkrf1 <- 10**sum(log10(rf1.recharge), res1.ok.ras.resampled)
proj4string(rkrf1) <- CRS("+proj=longlat +datum=WGS84")

# write rasters
writeRaster(rkrf5, filename="Data/Low_res_data_05/R_RKRF_recharge_05.tif", format="GTiff", overwrite=T)
writeRaster(rkrf1, filename="Data/High_res_data_01/R_RKRF_recharge_01.tif", format="GTiff", overwrite=T)
