rm(list=ls())
require(data.table)
require(sf)
require(raster)

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
