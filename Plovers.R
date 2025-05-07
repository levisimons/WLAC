#Background on code is here: https://docs.google.com/document/d/1khr6h0livaJWeVYUU2aVJ1YYye_mJeBsD6sJNt3_Z8U/edit?usp=sharing
rm(list=ls())
require(data.table)
require(sf)
require(raster)
require(stars)
require(terra)
require(dplyr)
require(plyr)
require(dismo)
require(randomForest)
require(DescTools)
require(geodata)
require(ggplot2)

#Set working directory
wd <- "/Users/levisimons/Desktop/WLAC/Plovers"
setwd(wd)

#Set random number string
set.seed(1)

#Get English county boundaries from https://simplemaps.com/gis/country/gb
#Read in county shapefiles
England <- st_read("gb.shp")
#Retain only Cumbria
Cumbria <- England %>% dplyr::filter(name == "Cumbria")
#Reproject to EPSG:4326
Cumbria <- st_transform(Cumbria,crs=st_crs(4326))

#Read in soil map layers downloaded from https://soilgrids.org/
soil_vars <- list.files(path=".",pattern = "^Soil_")
soil_layer <- c()
j=1
for(soil_var in soil_vars){
  tmp <- raster(soil_var)
  #Transform the coordinates to a CRS of 4326
  tmp <- projectRaster(tmp,crs = "+init=epsg:4326")
  #Clip soil layers to the boundaries of Cumbria
  tmp <- crop(tmp,Cumbria)
  tmp <- mask(tmp,Cumbria)
  if(j>1){tmp <- resample(tmp,soil_layer[[1]], method='ngb')}
  soil_layer[[j]] <- tmp
  j=j+1
}
#Stack soil layers
soil_layers <- stack(soil_layer)
