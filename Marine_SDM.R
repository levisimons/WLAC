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

#Generate a set of random background points, equal in number to actual occurrences.
num_occurrences <- nrow(input_species_data)
background_points = sf::st_sample(scb_boundaries, size=num_occurrences)

#Convert single column coordinates to standard longitude/latitude columns
background_points <- sf::st_coordinates(background_points)

#Retain species occurrence data from within the SCB
scb_species_points <- st_intersection(input_species_points, scb_boundaries) 

#Convert background points object to a data table
background_points <- as.data.table(background_points)

#Convert from longitude (X) and latitude (Y) columns to sf object.
background_points <- sf::st_as_sf(background_points,coords = c("X","Y"),  crs = 4326)
