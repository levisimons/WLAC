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
