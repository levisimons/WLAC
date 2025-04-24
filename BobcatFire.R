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
#wind	wind speed (2 m above the ground)	m s-1
#bio for bioclimatic variables
env_vars <- c("wind","bio")
#Loop through each variable category and clip out climate data for LA National Forest.
i=1
environmental_layer <- c()
for(env_var in env_vars){
  #Read in environmental layer for the USA
  tmp <- worldclim_tile(var = env_var, lon = mean(c(LANF_bbox["xmin"], LANF_bbox["xmax"])), 
                        lat = mean(c(LANF_bbox["ymin"], LANF_bbox["ymax"])),
                        res = 0.5, path = tempdir())
  #Clip climate layers to the boundaries of LANF
  tmp <- crop(tmp,LANF)
  tmp <- mask(tmp,LANF)
  #Store environmental layer in list.
  if(env_var=="wind"){
    #Loop through each wind speed month variable
    for(j in 1:length(names(tmp))){
      environmental_layer[[i]] <- resample(raster(tmp[[j]]), burn_severity, method = "bilinear")
      i=i+1
    }
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

#Read in soil map layers downloaded from https://soilgrids.org/
soil_vars <- list.files(path=".",pattern = "^Soil_.*\\.tif$")
soil_layer <- c()
j=1
for(soil_var in soil_vars){
  tmp <- raster(soil_var)
  #Transform the coordinates to a CRS of 4326
  tmp <- projectRaster(tmp,crs = "+init=epsg:4326")
  #Clip soil layers to the boundaries of LA National Forest
  tmp <- crop(tmp,LANF)
  tmp <- mask(tmp,LANF)
  #Resample layers to align with other map layers
  if(names(tmp)=="Soil_Groups"){
    tmp <- resample(tmp, burn_severity, method = "ngb")
  } else{
    tmp <- resample(tmp, burn_severity, method = "bilinear")
  }
  soil_layer[[j]] <- tmp
  j=j+1
}
#Stack more environmental layers
soil_layers <- stack(soil_layer)
environmental_layers <- stack(environmental_layers,soil_layers)

#Extract raster values at occurrence points
presence_extracted <- as.data.frame(raster::extract(environmental_layers, presence_points))

#Remove empty rows
presence_extracted <- as.data.frame(presence_extracted[complete.cases(presence_extracted),])

#Add presence/absence column
presence_extracted$presence <- 1

#Set presence variable to factor for modeling.
presence_extracted$presence <- as.factor(presence_extracted$presence)

#Set soil groups variable to factor for modeling.
presence_extracted$Soil_Groups <- as.factor(as.integer(presence_extracted$Soil_Groups))

#Generate a set of random background points, equal in number to actual occurrences.
num_occurrences <- nrow(presence_extracted)
background_points = sf::st_sample(LANF, size=num_occurrences)

#Convert single column coordinates to standard longitude/latitude columns
background_points <- sf::st_coordinates(background_points)

#Convert background points object to a data table
background_points <- as.data.table(background_points)

#Convert from longitude (X) and latitude (Y) columns to sf object.
background_points <- sf::st_as_sf(background_points,coords = c("X","Y"),  crs = 4326)

#Extract raster values at background points
background_extracted <- as.data.frame(raster::extract(environmental_layers, background_points))

#Remove empty rows
background_extracted <- as.data.frame(background_extracted[complete.cases(background_extracted),])

#Add presence/absence column
background_extracted$presence <- 0

#Set presence variable to factor for modeling.
background_extracted$presence <- as.factor(background_extracted$presence)

#Set soil groups variable to factor for modeling.
background_extracted$Soil_Groups <- as.factor(as.integer(background_extracted$Soil_Groups))

#Create an empty list to store prediction rasters.
raster_predict_list <- c()
#Create an empty list to store relative importance outputs.
importance_list <- c()
#Create an empty list to store accuracy outputs.
accuracy_list <- c()
#Create an empty list to store partial plot outputs.
partial_plot_list <- c()
j <- 1
i_max <- 25
for(i in 1:i_max){
  #Create a subset of the presence/background data with the following properties:
  #1. Composed of a randomly selected 80% of rows from scb_extracted.
  #2. Composed of rows randomly selected from background_extracted. The number of rows will also be 80% of rows found in scb_extracted.
  #3. Merged these two subsets together.
  subset_extracted <- rbind(presence_extracted[sample(nrow(presence_extracted),0.8*nrow(presence_extracted)),],background_extracted[sample(nrow(background_extracted),0.8*nrow(presence_extracted)),])
  
  #Run a random forest model over this data subset.
  rf1 <- suppressWarnings(tuneRF(x=subset_extracted[,!(colnames(subset_extracted) %in% "presence")],y=subset_extracted$presence,stepFactor=1,plot=FALSE,doBest=TRUE))
  
  #Make a prediction raster from the random forest model and store it in a list.
  raster_predict_list[[i]] <- dismo::predict(environmental_layers,rf1,progress='text')
  #Plot predicted raster
  plot(raster_predict_list[[i]])
  
  #Store relative importance of variable outputs as a temporary data frame.
  tmp <- as.data.frame(rf1$importance)
  #Set one column to store the variable names from the row names.
  tmp$VariableName <- rownames(tmp)
  #Store this importance data frame in the importance list.
  importance_list[[i]] <- tmp
  
  #Calculate the true skill statistic TSS to evaluate model accuracy.
  sensitivity <- rf1$confusion[[1]] / (rf1$confusion[[1]]+rf1$confusion[[2]])
  specificity <- rf1$confusion[[4]] / (rf1$confusion[[4]]+rf1$confusion[[3]])
  TSS <- sensitivity+specificity-1
  #Store TSS results
  accuracy_list[i] <- TSS
  
  #Loop through each environmental variable and store the partial response outputs in a temporary data frame.
  for(environmental_layer in names(environmental_layers)){
    #Store partial plot chart data in a temporary data frame.
    tmp <- as.data.frame(partialPlot(rf1,subset_extracted[,!(colnames(subset_extracted) %in% "presence")],x.var=c(environmental_layer),plot=F))
    #Transform logistic probabilities to regular probabilities.
    tmp$y <- exp(tmp$y) / (1+exp(tmp$y))
    #Rename probability column
    colnames(tmp) <- c(environmental_layer,"Detection_Probability")
    #Store partial plot data in a list of data frames.
    partial_plot_list[[j]] <- tmp
    j <- j+1
  }
  print(paste(i,Sys.time()))
}
#Stack the list of prediction rasters.
raster_predict <- brick(raster_predict_list)
#Sum the list of prediction rasters.
raster_predict <- calc(raster_predict, sum)
#Save raster output.
writeRaster(raster_predict,paste("nicotiana_glauca_prediction_",year_selected,".tif",sep=""),overwrite=T)
#Plot raster output.
plot(raster_predict,col=viridis(i_max))

# Convert full raster to data frame
raster_df <- as.data.frame(raster_predict, xy = TRUE)
names(raster_df)[3] <- "value"

# Extract non-zero points
raster_points <- subset(raster_df, value > 0.5*i_max)

#Plot prediction raster
ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = value)) +
  scale_fill_gradient(low = "white", high = "white", na.value = "grey") +
  guides(fill = "none") +
  geom_point(data = raster_points, aes(x = x, y = y, color = value), size = 1) +
  scale_color_viridis_c() +
  coord_fixed() +
  theme_minimal() +
  labs(title = paste("Predicted occurrences of tree tobacco",sep=""),
       x="Longitude degrees East",y = "Latitude degrees North",
       color = paste("Predicted frequency\nout of ",i_max," models",sep=""))

#Convert list of importance data frames to a single data frame.
importance_total <- rbind.fill(importance_list)
#Calculate the mean relative importance for each variable.
importance_total <- aggregate(x=importance_total$MeanDecreaseGini,by = list(importance_total$VariableName),FUN = mean)
#Rename columns.
colnames(importance_total) <- c("VariableName","Importance")
#Convert importance to rank importance.
importance_total$Importance <- rank(desc(importance_total$Importance))
#Save rank importance table.
write.table(importance_total,paste("nicotiana_glauca_rank_importance",year_selected,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)

#Calculate the mean TSS for the models
mean(accuracy_list)
sd(accuracy_list)

#Collapse partial plot outputs into single data frame.
partial_plots <- rbind.fill(partial_plot_list)
partial_plots <- as.data.frame(partial_plots)
write.table(partial_plots,paste("nicotiana_glauca_partial_plots",year_selected,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
partial_plots <- read.table(paste("nicotiana_glauca_partial_plots",year_selected,".txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
#Set soil types to categorical for plotting
partial_plots$Soil_Groups <- as.factor(partial_plots$Soil_Groups)

k <- 29
#Plot heat maps for continuous data.
if(!is.factor(partial_plots[,names(environmental_layers[[k]])])){
  ggplot(partial_plots, aes(x=!!sym(names(environmental_layers[[k]])), y=Detection_Probability) )+
    xlab(names(environmental_layers[[k]]))+ylab("Detection\nProbability")+
    geom_bin2d(bins = 50)+
    scale_fill_continuous(type = "viridis",name=paste("Frequency\n(Out of ",i_max," models)",sep=""))+
    stat_smooth(aes(y = Detection_Probability, fill=Detection_Probability),method="auto",formula=y~x,color="violet",fill="red",n=0.1*sum(!is.na(partial_plots[,names(environmental_layers[[k]])])))+
    theme_bw(base_size=25)
}
#Plot violin plots for categorical data.
if(is.factor(partial_plots[,names(environmental_layers[[k]])])){
  ggplot(partial_plots, aes(x=!!sym(names(environmental_layers[[k]])), y=Detection_Probability) )+
    xlab(names(environmental_layers[[k]]))+ylab("Detection\nProbability")+
    geom_violin(aes(x=!!sym(names(environmental_layers[[k]])), y=Detection_Probability))+
    theme_bw(base_size=25)
}
