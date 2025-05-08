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
wd <- ""
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

#Note: all Bioclim variables have a CRS of EPSG:4326
#Get all environmental variables.
#srad	incident solar radiation	kJ m-2 day-1
#wind	wind speed (2 m above the ground)	m s-1
#bio for bioclimatic variables
env_vars <- c("srad","wind","bio")
#Loop through each variable category and clip out climate data for Cumbria.
i=1
environmental_layer <- c()
for(env_var in env_vars){
  #Read in environmental layer for Great Britain
  tmp <- worldclim_country("GBR",res=0.5,var=env_var,path=tempdir())
  #Clip climate layers to the boundaries of Cumbria
  tmp <- crop(tmp,Cumbria)
  tmp <- mask(tmp,Cumbria)
  #Loop through each variable group and store it in a list
  for(j in 1:length(names(tmp))){
    environmental_layer[[i]] <- resample(raster(tmp[[j]]), soil_layers[[1]], method = "bilinear")
    i=i+1
  }
}
#Stack environmental layers
environmental_layers <- stack(environmental_layer)
environmental_layers <- stack(environmental_layers,soil_layers)

#Filter collinear environmental variables
#https://onlinelibrary.wiley.com/doi/pdf/10.1002/ece3.10901
env_retain <- removeCollinearity(environmental_layers,method="spearman",
                                 multicollinearity.cutoff = 0.75,sample.points = TRUE,
                                 nb.points = 1000,select.variables = TRUE)

#Build a raster stack of all environmental rasters with filtered layers.
environmental_layers <- subset(environmental_layers, env_retain)
environmental_layers <- stack(environmental_layers)

#Read in plover occurrences from GBIF.
#GBIF.org (13 April 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.sk6enc
presence_points <- read.table("Plovers.csv", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

#Convert species data to a spatial points object
presence_points <- st_as_sf(presence_points, coords = c("decimalLongitude","decimalLatitude"),  crs = 4326) 

#Filter points closer than 9.26km
#Initialize output
keep <- rep(TRUE, nrow(presence_points))
#Loop through points and remove points closer than 9260 meters
for (i in seq_len(nrow(presence_points))) {
  if (!keep[i]) next
  
  dists <- st_distance(presence_points[i, ], presence_points)
  close_points <- which(as.numeric(dists) < 1000 & as.numeric(dists) > 0)
  keep[close_points] <- FALSE
}
#Retain only distance filtered points
presence_points <- presence_points[keep, ]

