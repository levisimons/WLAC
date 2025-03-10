rm(list=ls())
require(data.table)
require(sf)

#Set working directory
wd <- ""
setwd(wd)

#Intake species occurrence data from GBIF
input_species_data <- fread(input="macrocystis_pyrifera.csv",sep="\t")

#Convert species data to a spatial points object
input_species_points <- st_as_sf(input_species_data, coords = c("decimalLongitude","decimalLatitude"),  crs = 4326) 

#Read in 3-Nautical Mile Coastal Boundary, California, 2001
#Source: https://earthworks.stanford.edu/catalog/stanford-sg211gq3741 
#Relevant shapefile data is initially in the directory sg211gq3741/data_EPSG_4326 , move it to the same directory as the script.
california_coast <- st_read("CA_cst3nm.shp")

#Crop the California coastal boundary to a bounding box covering the Southern California Bight
#Northwestern point of study area is Point Conception (34.4486° N, 120.4716° W), southeastern point is Tijuana Beach Promenade (32.5340° N, 117.1235° W)
scb_boundaries <-  st_crop(california_coast,st_bbox(c(xmin=-120.4716,xmax=-117.1235,ymin=32.5340,ymax=34.4486)))
