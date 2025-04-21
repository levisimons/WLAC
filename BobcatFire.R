#Background on code is here: https://docs.google.com/document/d/1HYc_Lw2yV8T63o-ixHWLCDJ-BGCQbdnVac_Xjep4voA/edit?usp=sharing
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
require(viridis)

#Set working directory
wd <- ""
setwd(wd)

#Set random number string
set.seed(1)

#Read in Los Angeles National Forest boundaries from https://egis-lacounty.hub.arcgis.com/datasets/lacounty::national-forest/about 
#Read in county shapefiles
LANF <- st_read("National_Forest.shp")
#Reproject to EPSG:4326
LANF <- st_transform(LANF,crs=st_crs(4326))

#Get bounding box for LANF
LANF_bbox <- st_bbox(LANF)

#Select year
year_selected <- 2019

#Read in nicotiana_glauca occurrences from GBIF.
#GBIF.org (14 April 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.4733vf 
presence_points <- read.table("nicotiana_glauca.csv", header=TRUE, quote = "",sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

#Filter occurrences by year.
presence_points <- presence_points[presence_points$year==year_selected,]

#Convert species data to a spatial points object
presence_points <- st_as_sf(presence_points, coords = c("decimalLongitude","decimalLatitude"),  crs = 4326) 

#Clip species data to Los Angeles National Forest
presence_points <- st_intersection(presence_points, LANF)

#Read in burn severity layer
#Google Earth Engine export code: https://code.earthengine.google.com/4106294b83776bc9cf8e37aa0fc8d4c3
burn_severity <- raster(paste("mean_nbr_",year_selected,"_la.tif",sep=""))

#Clip burn severity to Los Angeles National Forest boundaries
burn_severity <- crop(burn_severity,LANF)
burn_severity <- mask(burn_severity,LANF)

#Read in NDVI layer
#Google Earth Engine export code: https://code.earthengine.google.com/d15d0aba4b09e7dca63307576d2e7bcb
ndvi <- raster(paste("mean_ndvi_",year_selected,"_la.tif",sep=""))

#Resample NDVI layer to match spatial resolution of burn severity layer
ndvi <- resample(ndvi, burn_severity, method = "bilinear")

#Clip burn severity to Los Angeles National Forest boundaries
ndvi <- crop(ndvi,LANF)
ndvi <- mask(ndvi,LANF)

#Note: all Bioclim variables have a CRS of EPSG:4326
#Get all environmental variables.
#srad	incident solar radiation	kJ m-2 day-1
#wind	wind speed (2 m above the ground)	m s-1
#bio for bioclimatic variables
env_vars <- c("wind","bio")
#Loop through each variable category and clip out climate data for Cumbria.
i=1
environmental_layer <- c()
for(env_var in env_vars){
  #Read in environmental layer for the USA
  tmp <- worldclim_tile(var = env_var, lon = mean(c(LANF_bbox["xmin"], LANF_bbox["xmax"])), 
                        lat = mean(c(LANF_bbox["ymin"], LANF_bbox["ymax"])),
                        res = 0.5, path = tempdir())
  #Clip climate layers to the boundaries of Cumbria
  tmp <- crop(tmp,LANF)
  tmp <- mask(tmp,LANF)
  #Store environmental layer in list.
  if(env_var!="bio"){
    environmental_layer[[i]] <- resample(raster(tmp), burn_severity, method = "bilinear")
    i=i+1
  }
  #Store all of the bioclimatic variables as well.
  if(env_var=="bio"){
    #Loop through each bioclimatic variable
    for(j in 1:length(names(tmp))){
      environmental_layer[[i]] <- resample(raster(tmp[[j]]), burn_severity, method = "bilinear")
      i=i+1
    }
  }
  
}
#Stack environmental layers
environmental_layers <- stack(environmental_layer)
environmental_layers <- stack(environmental_layers,burn_severity)
environmental_layers <- stack(environmental_layers,ndvi)
