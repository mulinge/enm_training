
#18/08/2021
# calculating and visualizing uncertainty between models


#=============================================================================================================
#INSTALLING PACKAGES-----------------

#install.packages ('devtools')
#install.packages('raster')
#install.packages('sp')
#install.packages('rgeos')
#install.packages ('rgdal')
#install.packages ('maptools')
#install.packages('mapdata')
#install.packages ('dismo')
#install.packages('spThin')
#install.packages ('ENMeval')

#if(!require(kuenm)){
#  devtools::install_github("marlonecobos/kuenm")
#}
#=============================================================================================================
#LIBRARIES----------------------------
library (devtools)
library(raster)
library(sp)
library(rgeos)
library (rgdal)
library (maptools)
library(mapdata)
library (dismo) #SDM algorithms plus other functions
library(spThin) #thinning of occurrences
library (ENMeval) #maxnet algorithm 
library (kuenm) #maxent algorithm dissected 
library (virtualspecies) #create virtual species
library(corrplot) #create correlation matrices
library (gatepoints) #for advance functions
library (rgl)
library (Rcpp)
library (hypervolume) #hypervolume based models 
#=============================================================================================================
setwd("E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana")
#=============================================================================================================
#VISUALIZATION OF MODELS---------------------------------
#median of bootstrap
#Read all the calculated medians from the Median RCPs folder output
rudgeana_current_median = raster ('E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/MedianRCPs/MedianRCPsrudgeanaoutput/rudgeana_current_median.asc')
rudgeana_4550_median = raster ('E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/MedianRCPs/MedianRCPsrudgeanaoutput/rudgeana_4550_median.asc')
rudgeana_4570_median = raster ('E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/MedianRCPs/MedianRCPsrudgeanaoutput/rudgeana_4570_median.asc')
rudgeana_8550_median = raster ('E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/MedianRCPs/MedianRCPsrudgeanaoutput/rudgeana_8550_median.asc')
rudgeana_8570_median = raster ('E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/MedianRCPs/MedianRCPsrudgeanaoutput/rudgeana_8570_median.asc')

#=============================================================================================================
#uncertainty depicted as range: 
#Read median of the max layers from their folder
rudgeana_currentmax_median = raster ('E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/Maxmedians/Maxmediansrudgeanaoutput/rudgeana_currentmax_median.asc')
rudgeana_futuremax4550_median = raster ('E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/Maxmedians/Maxmediansrudgeanaoutput/rudgeana_futuremax4550_median.asc')
rudgeana_futuremax4570_median = raster ('E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/Maxmedians/Maxmediansrudgeanaoutput/rudgeana_futuremax4570_median.asc')
rudgeana_futuremax8550_median = raster ('E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/Maxmedians/Maxmediansrudgeanaoutput/rudgeana_futuremax8550_median.asc')
rudgeana_futuremax8570_median = raster ('E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/Maxmedians/Maxmediansrudgeanaoutput/rudgeana_futuremax8570_median.asc')

#Read median of the min layers from their folder
rudgeana_currentMin_median = raster ('E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/Minmedians/Minmediansrudgeanaoutput/rudgeana_currentMin_median.asc')
rudgeana_futureMin4550_median = raster ('E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/Minmedians/Minmediansrudgeanaoutput/rudgeana_futureMin4550_median.asc')
rudgeana_futureMin4570_median = raster ('E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/Minmedians/Minmediansrudgeanaoutput/rudgeana_futureMin4570_median.asc')
rudgeana_futureMin8550_median = raster ('E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/Minmedians/Minmediansrudgeanaoutput/rudgeana_futureMin8550_median.asc')
rudgeana_futureMin8570_median = raster ('E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/Minmedians/Minmediansrudgeanaoutput/rudgeana_futureMin8570_median.asc')

# calculate the range
rudgeanacurrent_uncer_mod = rudgeana_currentmax_median - rudgeana_currentMin_median 
rudgeana4550_uncer_mod = rudgeana_futuremax4550_median - rudgeana_futureMin4550_median
rudgeana4570_uncer_mod = rudgeana_futuremax4570_median - rudgeana_futureMin4570_median
rudgeana8550_uncer_mod = rudgeana_futuremax8550_median - rudgeana_futureMin8550_median
rudgeana8570_uncer_mod = rudgeana_futuremax8570_median - rudgeana_futureMin8570_median
#=============================================================================================================
#color ramp for uncertainty: 
color = colorRampPalette(c('#c90020', 'azure2', '#0571b1'))
#=============================================================================================================

rudgeana_pts <- read.csv("E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/Uncertainities/rudgeana.csv")

#=============================================================================================================
#visualization
#current
dev.new()
par (mfrow = c(1,1))

plot (rudgeana_current_median, main = 'rudgeanacurrent_KUENM model')
plot (rudgeanacurrent_uncer_mod, col= color(100), main = 'rudgeanacurrent_uncertainty (range)')

par (mfrow = c(1,1))
plot (rudgeana_current_median, main = 'rudgeanacurrent_KUENM model')
points (rudgeana_pts[,2:3], col = 'red', pch = 4, cex = 0.3)
#points (ind_p2[,2:3],col = 'blue', pch = 3, cex = 0.3)

#zoom(rudgeana_current_median, ext=drawExtent(), new=TRUE, useRaster=T) #zooming in a particular part of the world

#writeRaster(rudgeanacurrent_uncer_mod, filename = 'rudgeana_current_uncertainity_range', format = 'ascii', 
#            bylayer = T, suffix = 'rudgeana_current_uncertainity_range')

writeRaster(rudgeanacurrent_uncer_mod, file="E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/Uncertainities/rudgeanacurrent_uncer_mod", format="ascii", overwrite=F)
#=============================================================================================================
#4550
dev.new()
par (mfrow = c(1,1))

plot (rudgeana_4550_median, main = 'rudgeana_4550_KUENM model')
plot (rudgeana4550_uncer_mod, col= color(100), main = 'rudgeana_4550_uncertainty (range)')

par (mfrow = c(1,1))
plot (rudgeana_4550_median, main = 'rudgeana_4550_KUENM model')
points (rudgeana_pts[,2:3], col = 'red', pch = 4, cex = 0.3)
#points (ind_p2[,2:3],col = 'blue', pch = 3, cex = 0.3)

#zoom(rudgeana_4550_median, ext=drawExtent(), new=TRUE, useRaster=T) #zooming in a particular part of the world

#writeRaster(rudgeana4550_uncer_mod, filename = 'rudgeana_4550_uncertainity_range', format = 'ascii', 
#            bylayer = T, suffix = 'rudgeana_4550_uncertainity_range')

writeRaster(rudgeana4550_uncer_mod, file="E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/Uncertainities/rudgeana4550_uncer_mod", format="ascii", overwrite=F)
#=============================================================================================================
#4570
dev.new()
par (mfrow = c(1,1))

plot (rudgeana_4570_median, main = 'rudgeana_4570_KUENM model')
plot (rudgeana4570_uncer_mod, col= color(100), main = 'rudgeana_4570_uncertainty (range)')

par (mfrow = c(1,1))
plot (rudgeana_4570_median, main = 'rudgeana_4570_KUENM model')
points (rudgeana_pts[,2:3], col = 'red', pch = 4, cex = 0.3)
#points (ind_p2[,2:3],col = 'blue', pch = 3, cex = 0.3)

#zoom(rudgeana_4570_median, ext=drawExtent(), new=TRUE, useRaster=T) #zooming in a particular part of the world

#writeRaster(rudgeana4570_uncer_mod, filename = 'rudgeana_4570_uncertainity_range', format = 'ascii', 
#            bylayer = T, suffix = 'rudgeana_4570_uncertainity_range')

writeRaster(rudgeana4570_uncer_mod, file="E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/Uncertainities/mazonum4570_uncer_mod", format="ascii", overwrite=F)
#=============================================================================================================
#8550
dev.new()
par (mfrow = c(1,1))

plot (rudgeana_8550_median, main = 'rudgeana_8550_KUENM model')
plot (rudgeana8550_uncer_mod, col= color(100), main = 'rudgeana_8550_uncertainty (range)')

par (mfrow = c(1,1))
plot (rudgeana_8550_median, main = 'rudgeana_8550_KUENM model')
points (rudgeana_pts[,2:3], col = 'red', pch = 4, cex = 0.3)
#points (ind_p2[,2:3],col = 'blue', pch = 3, cex = 0.3)

#zoom(rudgeana_8550_median, ext=drawExtent(), new=TRUE, useRaster=T) #zooming in a particular part of the world

#writeRaster(rudgeana8550_uncer_mod, filename = 'rudgeana_8550_uncertainity_range', format = 'ascii', 
#            bylayer = T, suffix = 'rudgeana_8550_uncertainity_range')

writeRaster(rudgeana8550_uncer_mod, file="E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/Uncertainities/rudgeana8550_uncer_mod", format="ascii", overwrite=F)
#=============================================================================================================
#8570
dev.new()
par (mfrow = c(1,1))

plot (rudgeana_8570_median, main = 'KUENM model')
plot (rudgeana8570_uncer_mod, col= color(100), main = 'rudgeana_8570_uncertainty (range)')

par (mfrow = c(1,1))
plot (rudgeana_8570_median, main = 'rudgeana_8570_KUENM model')
points (rudgeana_pts[,2:3], col = 'red', pch = 4, cex = 0.3)
#points (ind_p2[,2:3],col = 'blue', pch = 3, cex = 0.3)

#zoom(rudgeana_8570_median, ext=drawExtent(), new=TRUE, useRaster=T) #zooming in a particular part of the world

#writeRaster(rudgeana8570_uncer_mod, filename = 'rudgeana_8570_uncertainity_range', format = 'ascii', 
#            bylayer = T, suffix = 'rudgeana_8570_uncertainity_range')

writeRaster(rudgeana8570_uncer_mod, file="E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/Uncertainities/rudgeana8570_uncer_mod", format="ascii", overwrite=F)
#=============================================================================================================

#MODEL THRESHOLDS---------------------------------------
rudgeana2_pts <- read.csv("E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/Uncertainities/rudgeana.csv", header = T)

dim(rudgeana2_pts)

crs <- CRS("+proj=longlat +datum=WGS84")
rudgeana3_pts <- SpatialPointsDataFrame(rudgeana2_pts[,2:3], proj4string = crs, rudgeana2_pts)

#extracting the values from the model using all the calibration points
suit_vals_current = extract (rudgeana_current_median, rudgeana3_pts[,2:3]) 

#using the value of 5% as threshold 
th1 = quantile(suit_vals_current, 0.05) 

#collecting only those values above that threshold, automatically transforms in 1 presence 0 absence 
rudgeana_current_median_th1 = rudgeana_current_median >= th1

plot (rudgeana_current_median_th1)
points (rudgeana_pts[,2:3], col = 'red', pch = 4, cex = 0.3) #calibration points
#points (ind_p2[,2:3],col = 'blue', pch = 3, cex = 0.3) #independent points 

#writeRaster(rudgeana_current_median_th1, filename = 'rudgeana_current_median_th5', format = 'ascii', 
#            bylayer = T, suffix = 'rudgeana_current_median_th5') # output layer after setting a suitability scale equivalent to the omission threshold

writeRaster(rudgeana_current_median_th1, file="E:/Jmaster/Nymp2/Nymp_analyse/Analysis2/Analysis/Analysis/Nrudgeana/Uncertainities/rudgeana_current_median_th1", format="ascii", overwrite=F)
#=============================================================================================================
#=============================================================================================================


