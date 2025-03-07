rm(list=ls())
require(data.table)
require(sf)

#Set working directory
wd <- "" #Path to file on your own machine
setwd(wd)

#Intake species occurrence data from GBIF
#Available from querying GBIF: https://doi.org/10.15468/dl.j46jh7
#Copy of sample file here: https://github.com/levisimons/WLAC/blob/main/macrocystis_pyrifera.csv
input_species_data <- fread(input="macrocystis_pyrifera.csv",sep="\t")

#Convert species data to a spatial points object
input_species_points <- st_as_sf(input_species_data, coords = c("decimalLongitude","decimalLatitude"),  crs = 4326) 
