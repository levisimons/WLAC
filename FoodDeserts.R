#Background on code is here: https://docs.google.com/document/d/11_zQAQ2co13Ac9W0ZCp0K-YhFZ6vb1DocPAndwitIQU/edit?usp=sharing
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

#Set working directory
wd <- ""
setwd(wd)

#Set random number string
set.seed(1)

#Read in community garden coordinates
CommunityGardens <- fread(input="CommunityGardens.csv",sep=",")
#Convert community garden coordinates to spatial coordinates with a CRS of 4326
CommunityGardens <- st_as_sf(CommunityGardens,coords=c("Longitude","Latitude"),crs = st_crs(4326))
#Transform the coordinates to a CRS of 2229 to match the map layers
CommunityGardens <- st_transform(CommunityGardens,crs=st_crs(2229))

#Read in food deserts map shapefiles
Food_Deserts_Input <- st_read("Food_Deserts.shp")

#Choose variables to use in modeling
selected_variables <- c("PovertyRat","LILATrac_1")
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
environmental_variables <- c("CIscore")
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

#Read in Healthy Places Index shapefile
Healthy_Places_Index <- st_read("Healthy_Places_Index_(3.0).shp")
#Transform the coordinates to a CRS of 2229 to match the map layers
Healthy_Places_Index <- st_transform(Healthy_Places_Index,crs=st_crs(2229))
#Choose variables to use in modeling
HPI_variables <- c("LEB","bachelor_1","insured_pc","transpor_1","automobi_1","parkacce_1")
#Initialize index and list of map layers object
k=1
HPI_layer <- c()
#Loop through variables to rasterize from the food desert map shapefiles. Store outputs in a list of map rasters.
for(HPI_variable in HPI_variables){
  #Set raster values of -999 to NA before saving them.
  tmp <-raster(rast(st_rasterize(Healthy_Places_Index %>% dplyr::select(!!HPI_variable, geometry))))
  values(tmp)[values(tmp) <= -999] = NA
  #Clip HPI layers to match the other map layers (boundaries of LA county)
  tmp <- crop(tmp,Food_Deserts_Input)
  tmp <- mask(tmp,Food_Deserts_Input)
  #Resample HPI layers to match to match those of the other food desert layers
  tmp <- resample(tmp, food_desert_variables[[1]], method = "bilinear")  # Use "near" for categorical data
  HPI_layer[[k]] <- tmp
  k=k+1
}

#Stack all of the CAL Enviro Screen variables into a single object
HPI_layers <- stack(HPI_layer)

#Add in CAL Enviro Screen layers to the model input stack
food_desert_variables <- stack(food_desert_variables,HPI_layers)

#Rename map layers in raster stack
names(food_desert_variables) <- c(selected_variables,"Distance_From_Community_Gardens",environmental_variables,HPI_variables)
food_desert_layers <- names(food_desert_variables)[!(names(food_desert_variables) %in% "LEB")]

#Generate a set of random background points, equal in number to actual occurrences.
num_occurrences <- 1000
background_points = sf::st_sample(Food_Deserts_Input, size=num_occurrences)

#Convert single column coordinates to standard longitude/latitude columns
background_points <- sf::st_coordinates(background_points)

#Convert background points object to a data table
background_points <- as.data.table(background_points)

#Convert from longitude (X) and latitude (Y) columns to sf object.
background_points <- sf::st_as_sf(background_points,coords = c("X","Y"),  crs = 2229)

#Extract raster values at occurrence points
food_desert_extracted <- as.data.frame(raster::extract(food_desert_variables, background_points))

#Set food desert variable to categorical
food_desert_extracted$LILATrac_1 <- as.factor(food_desert_extracted$LILATrac_1)

#Remove empty rows
food_desert_extracted <- as.data.frame(food_desert_extracted[complete.cases(food_desert_extracted),])

#Create an empty list to store prediction rasters.
raster_predict_list <- c()
#Create an empty list to store relative importance outputs.
importance_list <- c()
#Create an empty list to store accuracy outputs.
accuracy_list <- c()
#Create an empty list to store partial plot outputs.
partial_plot_list <- c()
j <- 1
i_max <- 100
for(i in 1:i_max){
  #Create a subset of the food desert data with the following properties:
  #1. Composed of a randomly selected 80% of rows from food_desert_extracted.
  #Construct a training and testing set for the streams subset data.
  group <- kfold(food_desert_extracted,5)
  food_desert_subset <- food_desert_extracted[group!=1,]
  food_desert_test <- food_desert_extracted[group==1,]
  
  #Run a random forest model over this data subset.
  rf1 <- suppressWarnings(tuneRF(x=food_desert_subset[,!(colnames(food_desert_subset) %in% "LEB")],y=food_desert_subset$LEB,stepFactor=1,plot=FALSE,doBest=TRUE))
  
  #Make a prediction raster from the random forest model and store it in a list.
  raster_predict_list[[i]] <- dismo::predict(food_desert_variables,rf1,progress='text')
  #Plot predicted raster
  plot(raster_predict_list[[i]])
  
  #Store relative importance of variable outputs as a temporary data frame.
  tmp <- as.data.frame(rf1$importance)
  #Set one column to store the variable names from the row names.
  tmp$VariableName <- rownames(tmp)
  #Store this importance data frame in the importance list.
  importance_list[[i]] <- tmp
  
  #Calculate the correlation between actual and predicted life expectancy to evaluate model accuracy.
  #Calculate the predicted life expectancy
  prediction <- predict(rf1,food_desert_test)
  #Store correlation results
  accuracy_list[i] <- cor.test(prediction,food_desert_test$LEB)$estimate
  
  #Loop through each environmental variable and store the partial response outputs in a temporary data frame.
  for(food_desert_layer in food_desert_layers){
    #Store partial plot chart data in a temporary data frame.
    tmp <- as.data.frame(partialPlot(rf1,food_desert_subset[,!(colnames(food_desert_subset) %in% "LEB")],x.var=c(food_desert_layer),plot=F))
    #Rename probability column
    colnames(tmp) <- c(food_desert_layer,"Life Expectancy at Birth")
    #Store partial plot data in a list of data frames.
    partial_plot_list[[j]] <- tmp
    j <- j+1
  }
  print(paste(i,Sys.time()))
}

#Stack the list of prediction rasters.
raster_predict <- brick(raster_predict_list)
#Sum the list of prediction rasters.
raster_predict <- calc(raster_predict, mean)

#Get the mean Pearson correlation coefficient between predicted and actual life expectancy
FisherZInv(mean(FisherZ(accuracy_list)))
FisherZInv(sd(FisherZ(accuracy_list)))

#Convert list of importance data frames to a single data frame.
importance_total <- rbind.fill(importance_list)
#Calculate the mean relative importance for each variable.
importance_total <- aggregate(x=importance_total$IncNodePurity,by = list(importance_total$VariableName),FUN = mean)
#Rename columns.
colnames(importance_total) <- c("VariableName","Importance")
#Convert importance to rank importance.
importance_total$Importance <- rank(importance_total$Importance)

#Collapse partial plot outputs into single data frame.
partial_plots <- rbind.fill(partial_plot_list)
partial_plots <- as.data.frame(partial_plots)

#Make heat map partial dependence plots for life expectancy
k <- 1
ggplot(partial_plots, aes(x=!!sym(food_desert_layers[k]), y=`Life Expectancy at Birth`) )+
  xlab(food_desert_layers[k])+ylab("`Life Expectancy\nat Birth`")+
  geom_bin2d(bins = 50)+
  scale_fill_continuous(type = "viridis",name=paste("Frequency\n(Out of ",i_max," models)",sep=""))+
  stat_smooth(aes(y = `Life Expectancy at Birth`, fill=`Life Expectancy at Birth`),method="auto",formula=y~x,color="violet",fill="red",n=0.1*sum(!is.na(partial_plots[,food_desert_layers[k]])))+
  theme_bw(base_size=25)
