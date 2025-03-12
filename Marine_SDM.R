rm(list=ls())
require(data.table)
require(sf)
require(raster)

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

#Convert background points object to a data table
background_points <- as.data.table(background_points)

#Convert from longitude (X) and latitude (Y) columns to sf object.
background_points <- sf::st_as_sf(background_points,coords = c("X","Y"),  crs = 4326)

#Retain species occurrence data from within the SCB
scb_species_points <- st_intersection(input_species_points, scb_boundaries) 

#Get a list of all environmental map layers with the file extension .nc (NetCDF)
map_layers <- list.files(path="MapLayers",pattern = "\\.nc$")

#Convert all of the NetCDF map layers from Bio-Oracle into clipped rasters in tif format.
#Check if list of tif layers has already been generated first.
if(length(list.files(path="MapLayers",pattern = "\\.tif$"))<length(map_layers)){
  for(map_layer in map_layers){
    #Read in marine layer
    marine_layer <- brick(paste("MapLayers/",map_layer,sep=""))
    #Get base file name
    map_layer_name <- gsub("\\.nc$","", map_layer)
    #Convert from NetCDF to tif format.
    writeRaster(marine_layer, paste("MapLayers/",map_layer_name,".tif",sep=""), bylayer=FALSE)
    #Read in tiff formatted raster.
    marine_raster <- raster(paste("MapLayers/",map_layer_name,".tif",sep=""))
    #Crop the raster to the study extent.
    scb_raster <- crop(marine_raster,scb_boundaries)
    #Export the cropped raster.
    writeRaster(scb_raster,paste("MapLayers/",map_layer_name,".tif",sep=""), bylayer=FALSE,overwrite=TRUE)
  }
}

#Get a list of all environmental rasters in tif format.
scb_layers <- list.files(path="MapLayers",pattern = "\\.tif$")

#Build a raster stack of all environmental rasters
scb_rasters <- stack(paste("MapLayers/",scb_layers,sep=""))

#Update column names so the column names match the environmental raster file names
names(scb_rasters) <- scb_layers

#Remove empty rows
scb_extracted <- as.data.frame(scb_extracted[complete.cases(scb_extracted),])

#Add presence column
scb_extracted$presence <- 1

#Set presence variable to factor for modeling.
scb_extracted$presence <- as.factor(scb_extracted$presence)

#Extract raster values at background points
background_extracted <- raster::extract(scb_rasters, background_points)

#Update column names so the column names match the environmental raster file names
colnames(background_extracted) <- scb_layers

#Remove empty rows
background_extracted <- as.data.frame(background_extracted[complete.cases(background_extracted),])

#Add presence column
background_extracted$presence <- 0

#Set presence variable to factor for modeling.
background_extracted$presence <- as.factor(background_extracted$presence)

#Set a predictable random number generator seed for reproducibility.
set.seed(1)

#Create a subset of the presence/background data with the following properties:
#1. Composed of a randomly selected 80% of rows from scb_extracted.
#2. Composed of rows randomly selected from background_extracted. The number of rows will also be 80% of rows found in scb_extracted.
#3. Merged these two subsets together.
subset_extracted <- rbind(scb_extracted[sample(nrow(scb_extracted),0.8*nrow(scb_extracted)),],background_extracted[sample(nrow(background_extracted),0.8*nrow(scb_extracted)),])

#Run a random forest model over this data subset.
rf1 <- suppressWarnings(tuneRF(x=subset_extracted[,!(colnames(subset_extracted) %in% "presence")],y=subset_extracted$presence,stepFactor=1,plot=FALSE,doBest=TRUE))

#Make a prediction raster from the random forest model and store it in a list.
raster_predict_list[[i]] <- dismo::predict(scb_rasters,rf1,progress='text')

#Plot predicted raster
plot(raster_predict_list[[i]])
