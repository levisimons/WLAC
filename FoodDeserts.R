#Background on code is here: https://docs.google.com/document/d/11_zQAQ2co13Ac9W0ZCp0K-YhFZ6vb1DocPAndwitIQU/edit?usp=sharing
rm(list=ls())
require(data.table)
require(sf)
require(raster)
require(stars)
require(terra)
require(dplyr)

#Set working directory
wd <- ""
setwd(wd)

#Read in community garden coordinates
CommunityGardens <- fread(input="CommunityGardens.csv",sep=",")
#Convert community garden coordinates to spatial coordinates with a CRS of 4326
CommunityGardens <- st_as_sf(CommunityGardens,coords=c("Longitude","Latitude"),crs = st_crs(4326))
#Transform the coordinates to a CRS of 2229 to match the map layers
CommunityGardens <- st_transform(CommunityGardens,crs=st_crs(2229))

#Read in food deserts map shapefiles
Food_Deserts_Input <- st_read("Food_Deserts.shp")

#Choose variables to use in modeling
selected_variables <- c("PovertyRat","LA1and10")
#Initialize index and list of map layers object
i=1
selected_layer <- c()
#Loop through variables to rasterize from the food desert map shapefiles. Store outputs in a list of map rasters.
for(selected_variable in selected_variables){
  selected_layer[[i]] <- raster(rast(st_rasterize(Food_Deserts_Input %>% dplyr::select(!!selected_variable, geometry))))
  i=i+1
}

#Stack all of the food desert variables so far into a single object
food_desert_variables <- stack(selected_layer)

#Set a blank raster which matches the shape and extent of the map layers
r <- rast(selected_layer[[1]])
values(r) <- NA
#Create a map layer which calculates distances from community gardens
distance_community_gardens <- distanceFromPoints(raster(r), st_coordinates(CommunityGardens))
#Clip distance to community gardens layer to match the other map layers (boundaries of LA county)
distance_community_gardens <- crop(distance_community_gardens, Food_Deserts_Input)
distance_community_gardens <- mask(distance_community_gardens,Food_Deserts_Input)

#Add in distance to community gardens layer to the model input stack
food_desert_variables <- stack(food_desert_variables,distance_community_gardens)

#Read in CAL Enviro Screen shapefiles
CAL_Enviro_Screen <- st_read("CES4 Final Shapefile.shp")
#Transform the coordinates to a CRS of 2229 to match the map layers
CAL_Enviro_Screen <- st_transform(CAL_Enviro_Screen,crs=st_crs(2229))

#Set environmental variables to rasterize
environmental_variables <- c("CIscore","Educatn")
#Initialize index and list of map layers object
j=1
environmental_layer <- c()
#Loop through variables to rasterize from the food desert map shapefiles. Store outputs in a list of map rasters.
for(environmental_variable in environmental_variables){
  #Set raster values of -999 to NA before saving them.
  tmp <-raster(rast(st_rasterize(CAL_Enviro_Screen %>% dplyr::select(!!environmental_variable, geometry))))
  values(tmp)[values(tmp) <= -999] = NA
  #Clip CAL Enviro Screen layers to match the other map layers (boundaries of LA county)
  tmp <- crop(tmp, Food_Deserts_Input)
  tmp <- mask(tmp,Food_Deserts_Input)
  #Resample CAL Enviro Screen layers to match to match those of the other food desert layers
  tmp <- resample(tmp, food_desert_variables[[1]], method = "bilinear")  # Use "near" for categorical data
  environmental_layer[[j]] <- tmp
  j=j+1
}

#Stack all of the CAL Enviro Screen variables into a single object
CAL_Enviro_Screen_layers <- stack(environmental_layer)

#Add in CAL Enviro Screen layers to the model input stack
food_desert_variables <- stack(food_desert_variables,CAL_Enviro_Screen_layers)

#Rename map layers in raster stack
names(food_desert_variables) <- c(selected_variables,"Distance_From_Community_Gardens",environmental_variables)
