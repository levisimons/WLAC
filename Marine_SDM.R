rm(list=ls())
require(data.table)
require(sf)
require(raster)
require(dismo)
require(randomForest)
require(viridis)
require(ggplot2)

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

#Extract raster values at occurrence points
scb_extracted <- raster::extract(scb_rasters, scb_species_points) 

#Update column names so the column names match the environmental raster file names
colnames(scb_extracted) <- scb_layers

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

#Create an empty list to store prediction rasters.
raster_predict_list <- c()
#Create an empty list to store relative importance outputs.
importance_list <- c()
#Create an empty list to store accuracy outputs.
accuracy_list <- c()
#Create an empty list to store partial plot outputs.
partial_plot_list <- c()
j <- 1
i_max <- 5
for(i in 1:i_max){
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
  for(scb_layer in scb_layers){
    #Store partial plot chart data in a temporary data frame.
    tmp <- as.data.frame(partialPlot(rf1,subset_extracted[,!(colnames(subset_extracted) %in% "presence")],x.var=c(scb_layer),plot=F))
    #Transform logistic probabilities to regular probabilities.
    tmp$y <- exp(tmp$y) / (1+exp(tmp$y))
    #Rename probability column
    colnames(tmp) <- c(scb_layer,"Detection_Probability")
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
writeRaster(raster_predict,"macrocystis_pyrifera_prediction.tif",overwrite=T)

# Convert sum prediction raster into data frame for mapping.
raster_predict_df <- as.data.frame(rasterToPoints(raster_predict, xy=TRUE))
# Map the frequency with which a species shows up across the study area.
ggplot() + 
  geom_sf(data = scb_boundaries, size = 1.5, color = "white", fill = "white") + 
  ggtitle("Prediction frequency of Macrocystis pyrifera\nIn the Southern California Bight") + 
  geom_raster(data = raster_predict_df, aes(x=x, y=y, fill=layer))+
  scale_fill_gradientn(colours=c("blue","orange"),name=paste("Frequency\n(Out of ",i_max," models)",sep=""),transform="log1p")

#Calculate the mean TSS for the models
mean(accuracy_list)
sd(accuracy_list)

#Convert list of importance data frames to a single data frame.
importance_total <- rbind.fill(importance_list)
#Calculate the mean relative importance for each variable.
importance_total <- aggregate(x=importance_total$MeanDecreaseGini,by = list(importance_total$VariableName),FUN = mean)
#Rename columns.
colnames(importance_total) <- c("VariableName","Importance")
#Convert importance to rank importance.
importance_total$Importance <- rank(importance_total$Importance)
#Save rank importance table.
write.table(importance_total,"macrocystis_pyrifera_rank_importance.txt",quote=FALSE,sep="\t",row.names = FALSE)

#Collapse partial plot outputs into single data frame.
partial_plots <- rbind.fill(partial_plot_list)
partial_plots <- as.data.frame(partial_plots)
write.table(partial_plots,"macrocystis_pyrifera_partial_plots.txt",quote=FALSE,sep="\t",row.names = FALSE)
partial_plots <- read.table("macrocystis_pyrifera_partial_plots.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
k <- 10
ggplot(partial_plots, aes(x=!!sym(scb_layers[k]), y=Detection_Probability) )+
  xlab(scb_layers[k])+ylab("Detection\nProbability")+
  geom_bin2d(bins = 50)+
  scale_fill_continuous(type = "viridis",name=paste("Frequency\n(Out of ",i_max," models)",sep=""))+
  stat_smooth(aes(y = Detection_Probability, fill=Detection_Probability),method="auto",formula=y~x,color="violet",fill="red",n=0.1*sum(!is.na(partial_plots[,scb_layers[k]])))+
  theme_bw(base_size=25)
