rm(list=ls())
require(data.table)
require(dplyr)
require(lubridate)
require(ggplot2)
require(jsonlite)
detach("package:tigris", unload = TRUE) #Run in case of issues with httr
require(httr)

#Set working directory
#Put in the file path to where you want to work on your computer
wd <- ""
setwd(wd)

#Read in USGS earthquake data for earthquakes with a greater magnitude than 4.0.
#Data obtained via https://earthquake.usgs.gov/earthquakes/search/
earthquakes <- fread(input="USGS_2025.csv",sep=",")

#Filter earthquakes based on magnitude.
magnitude_minimum <- 7
earthquakes <- earthquakes[earthquakes$mag >= magnitude_minimum]

#Set date column to a standard format
earthquakes$time <- as.Date(format(ymd_hms(earthquakes$time), "%Y-%m-%d"),format="%Y-%m-%d")

#Returns a data frame with the average daily temperatures in Celsius
#computed from 24 hourly readings at 2 m above ground.
fetch_avg_temp <- function(lat, lon, date, units = "celsius") {
  
  #Process dates and return errors in case of missing data or improper formatting.
  date_parsed <- tryCatch(
    as.Date(date, format = "%Y-%m-%d"),
    error = function(e) stop("date must be in YYYY-MM-DD format")
  )
  if (is.na(date_parsed)) stop("date must be in YYYY-MM-DD format")
  if (date_parsed < as.Date("1940-01-01"))
    stop("Open-Meteo archive starts 1940-01-01")
  if (date_parsed > Sys.Date() - 5)
    stop("Date must be at least 5 days in the past")
  
  #Define acceptable temperature units.
  units <- match.arg(units, c("celsius", "fahrenheit"))
  
  #Build API request to Open Meteo
  base_url <- "https://archive-api.open-meteo.com/v1/archive"
  
  #Run API call against Open Meteo.
  resp <- GET(
    url   = base_url,
    query = list(
      latitude         = lat,
      longitude        = lon,
      start_date       = as.character(date_parsed),
      end_date         = as.character(date_parsed),
      hourly           = "temperature_2m",
      temperature_unit = units,
      timezone         = "auto"          # local time for the coordinates
    )
  )
  
  #If there's an error with connecting with the API return it and stop the function.
  if (http_error(resp)) {
    body <- tryCatch(
      fromJSON(content(resp, as = "text", encoding = "UTF-8"))$reason,
      error = function(e) http_status(resp)$message
    )
    stop(sprintf("API error %d: %s", status_code(resp), body))
  }
  
  #Parse the JSON returned from the API call. Temperatures will be extracted from this.
  payload <- fromJSON(content(resp, as = "text", encoding = "UTF-8"))
  
  #Extract air temperature and time from JSON export.
  temps <- payload$hourly$temperature_2m
  times <- payload$hourly$time
  
  if (all(is.na(temps)))
    stop("No temperature data returned — check coordinates and date.")
  
  #Compute average temperature for location and date.
  avg_temp <- mean(temps, na.rm = TRUE)
  
  #Return average temperature for earthquake location
  result <- data.frame(
    date       = as.Date(date),
    latitude   = lat,
    longitude  = lon,
    avg_temp   = round(avg_temp, 2),
    stringsAsFactors = FALSE
  )
  result
}

#Initialize and empty list to store temperature results.
results <- c()
#Run a loop finding the average temperature associated with each earthquake.
for(i in 1:nrow(earthquakes)){
  #Supply earthquake coordinates
  lat  <- earthquakes[i,]$latitude
  lon  <- earthquakes[i,]$longitude
  date <- earthquakes[i,]$time
  
  #Get average temperature per earthquake
  result <- fetch_avg_temp(lat, lon, date, units = "celsius")
  
  #Print how far the loop has gone through all data.
  print(paste(i,nrow(earthquakes)))
  results[[i]] <- result
}

#Combine average temperature outputs
results <- rbindlist(results)

#Merge in average temperature data with earthquake data
earthquakes <- dplyr::left_join(earthquakes,results,by=c("time"="date","latitude"="latitude","longitude"="longitude"))

#Plot the magnitude of earthquakes versus average temperature
ggplot(data=earthquakes,aes(x=avg_temp, y=mag))+
  xlab("Average temperature in degrees C")+ylab("Earthquake magnitude")+
  geom_point(aes(x=avg_temp, y=mag))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

#Test if variations in earthquake magnitudes are normally distributed
#This test is needed to see which correlation test is appropriate
#for checking for significant correlations between earthquake magnitudes and average temperature
#The test used is a Kolmogorov-Smirnov test
ks.test(earthquakes$mag,"pnorm")

#Check for correlations between earthquake magnitudes and average temperature.
cor.test(earthquakes$avg_temp,earthquakes$mag,method="spearman")
