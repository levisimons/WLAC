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
