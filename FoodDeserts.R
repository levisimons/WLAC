rm(list=ls())
require(data.table)
require(sf)
require(raster)
require(stars)
require(terra)

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

selected_variables <- c("PovertyRat","LA1and10")
i=1
selected_layer <- c()
for(selected_variable in selected_variables){
  selected_layer[[i]] <- st_rasterize(Food_Deserts_Input %>% dplyr::select(!!selected_variable, geometry))
  i=i+1
}

#Set a blank raster which matches the shape and extent of the map layers
r <- rast(selected_layer[[1]])
values(r) <- NA
#Create a map layer which calculates distances from community gardens
distance_community_gardens <- distanceFromPoints(raster(r), st_coordinates(CommunityGardens))
#Clip distance to community gardens layer to match the other map layers (boundaries of LA county)
distance_community_gardens <- crop(distance_community_gardens, Food_Deserts_Input)
distance_community_gardens <- mask(distance_community_gardens,Food_Deserts_Input)
