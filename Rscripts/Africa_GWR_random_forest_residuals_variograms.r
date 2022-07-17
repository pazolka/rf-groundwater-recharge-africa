# source("Africa_GWR_random_forest_residuals_variograms.r")

af.shape <- readOGR("CRU_precip_CGIAR_AI_data/Africa_continent_shape.shp")
af.shp.utm <- spTransform(af.shape, CRS("+proj=utm +zone=35 +south +datum=WGS84"))

library(gstat); library(sp); library(rgdal); library(RColorBrewer); library(ggplot2); 
library(raster); library(rasterize)

res <- read.csv("residuals.csv", header=T)
res5.wgs <- res[,c(1,2,3)]
names(res5.wgs) <- c("X", "Y", "RES05")

# Setting existing coordinate as lat-long system

coordinates(res5.wgs) <- ~X+Y
proj4string(res5.wgs) <- CRS("+proj=longlat +datum=WGS84")
head(res5.wgs)
plot(res5.wgs)

res5.utm <- spTransform(res5.wgs, CRS("+proj=utm +zone=35 +south +datum=WGS84"))
# Senegal - two points for the same location with nearly the same values - take one only
res5.utm <- res5.utm[-zerodist(res5.utm)[,1],]

head(res5.utm)
plot(res5.utm)

# spplot(res5.utm["RES05"], col.regions=brewer.pal(9,"RdBu"), scales=list(draw=TRUE))

# par(mfrow=c(2,2), oma=c(0,0,0,0))

sel <- !is.infinite(res5.utm$RES05)
vario5 <- variogram(RES05 ~ 1, res5.utm[sel,], width=100000, cutoff=950000)
plot(vario5, pch=19, ylim=c(0.0, 0.12), cex=2, col=4)

vario.mod5 <- vgm(psill=0.2, model="Exp", nugget=0.45, range=200000)      # model variogram
# plot(vario5, pch=19, ylim=c(0.0, 0.12), cex=2, col=4)

fit.vario5 <- fit.variogram(vario5, model=vario.mod5,  fit.kappa=T)       # fit model variogram
plot(vario5, model=fit.vario5, lwd=2, pch=20, cex=2, col=4)   

attr(fit.vario5, "SSErr")

###############################################################################################

# res <- read.csv("RF_residuals.csv", header=T)
res1.wgs <- res[,c(1,2,4)]
names(res1.wgs) <- c("X", "Y", "RES01")

# Setting existing coordinate as lat-long system

coordinates(res1.wgs) <- ~X+Y
proj4string(res1.wgs) <- CRS("+proj=longlat +datum=WGS84")
head(res1.wgs)
plot(res1.wgs)

res1.utm <- spTransform(res1.wgs, CRS("+proj=utm +zone=35 +south +datum=WGS84"))
# Senegal - two points for the same location with nearly the same values - take one only
res1.utm <- res1.utm[-zerodist(res1.utm)[,1],]
head(res1.utm)
plot(res1.utm)

# spplot(res1.utm["RES01"], col.regions=brewer.pal(9,"RdBu"), scales=list(draw=TRUE))

# par(mfrow=c(2,2), oma=c(0,0,0,0))

sel <- !is.infinite(res1.utm$RES01)
vario1 <- variogram(RES01 ~ 1, res1.utm[sel,], width=75000, cutoff=900000)
plot(vario1, pch=19, ylim=c(0.0, 0.12), cex=2, col=4)

vario.mod1 <- vgm(psill=0.02, model="Exp", nugget=0.02, range=200000)    # model variogram
# plot(vario1, pch=19, ylim=c(0.0, 0.12), cex=2, col=4)

fit.vario1 <- fit.variogram(vario1, model=vario.mod1, fit.kappa=T)       # fit model variogram
plot(vario1, model=fit.vario1, lwd=2, pch=20, cex=2, col=4)   

attr(fit.vario1, "SSErr")

res.df <- as.data.frame(res1.utm)

###### # make grid 0.5 deg - alternative to minimize the information loss when reprojecting
rf5.recharge <- raster("Low_res_data_05/RF_recharge.tif")
rf5.recharge.afr <- crop(rf5.recharge, af.shape)
rf5.recharge.utm <- projectRaster(rf5.recharge.afr, crs='+proj=utm +zone=35 +south +datum=WGS84')
rf5.pixels <- as(rf5.recharge.utm, 'SpatialPixels')
plot(rf5.pixels)

####### krige at 0.5 deg
res5.ok <- gstat::krige(RES05 ~ 1, res5.utm, newdata=rf5.pixels, fit.vario5)

###### # make grid 0.1 deg - alternative to minimize the information loss when reprojecting
rf1.recharge <- raster("High_res_data_01/RF_recharge_01.tif")
rf1.recharge.afr <- crop(rf1.recharge, af.shape)
rf1.recharge.utm <- projectRaster(rf1.recharge.afr, crs='+proj=utm +zone=35 +south +datum=WGS84')
rf1.pixels <- as(rf1.recharge.utm, 'SpatialPixels')
plot(rf1.pixels)

####### krige at 0.1 deg
res1.ok <- gstat::krige(RES01 ~ 1, res1.utm, newdata=rf1.pixels, fit.vario1)

# transform to decimal degrees and write rasters
res5.ok.ras <- raster(res5.ok)
res1.ok.ras <- raster(res1.ok)

res5.ok.ras <- projectRaster(res5.ok.ras, crs="+proj=longlat +datum=WGS84")
res1.ok.ras <- projectRaster(res1.ok.ras, crs="+proj=longlat +datum=WGS84")

# plot(res.ok.ras, col=tim.colors(50))
# plot(af.shp.utm, add=T)

####### combine RF predictions (python) and residual kriging (R)

# resample bc resolutions and extents are slightly off
res5.ok.ras.resampled = resample(res5.ok.ras, rf5.recharge, method='bilinear')
res1.ok.ras.resampled = resample(res1.ok.ras, rf1.recharge, method='bilinear')

plot(res1.ok.ras)
plot(res1.ok.ras.resampled)
# add
rkrf5 <- 10**sum(log10(rf5.recharge), res5.ok.ras.resampled)
rkrf1 <- 10**sum(log10(rf1.recharge), res1.ok.ras.resampled)

# write rasters
writeRaster(rkrf5, filename="Low_res_data_05/R_RKRF_recharge.tif", format="GTiff", overwrite=T)
writeRaster(rkrf1, filename="High_res_data_01/R_RKRF_recharge.tif", format="GTiff", overwrite=T)

writeRaster(res5.ok.ras, filename="Low_res_data_05/RK_orig.tif", format="GTiff", overwrite=T)
writeRaster(res5.ok.ras.resampled, filename="Low_res_data_05/RK_res.tif", format="GTiff", overwrite=T)
