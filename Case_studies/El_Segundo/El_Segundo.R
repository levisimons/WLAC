rm(list=ls())
require(data.table)
require(sf)
require(leaflet)
require(dplyr)
require(tidyr)
require(terra)

#Set working directory. This is the path to where you want to work on your computer.
wd <- ""
setwd(wd)

#Set sf settings to reduce risk of spatial errors
sf_use_s2(FALSE)

#Read in project boundary
project_boundary <- st_read("Project_Boundary.shp")

#Read in butterfly counts by grid for 2019 - 2021
#Make sure all of the other files with the prefix ESBB_GridTotals_2019_2021 are also in the working directory.
grid_counts_1 <- st_read("ESBB_GridTotals_2019_2021.shp")
#Read in butterfly counts by grid for 2022 - 2025
#Make sure all of the other files with the prefix Block_Counts_within_Grids are also in the working directory.
grid_counts_2 <- st_read("Block_Counts_within_Grids.shp")

#Reproject butterfly the spatial object for the counts by grid for 2022 - 2025 to match
#those of 2019 - 2021. This is needed to eventually combine these data.
st_crs(grid_counts_2) <- st_crs(grid_counts_1)

#Aggregate butterfly counts by year and check for aggregate trends over time
#Make it seven rows, one for each year.
#The two columns are to store the year and the total count of butterflies from across the study area.
grid_count_summary <- data.frame(matrix(nrow=7,ncol=2))
colnames(grid_count_summary) <- c("year","total count")
grid_count_summary$year <- c(2019,2020,2021,2022,2023,2024,2025)
grid_count_summary[grid_count_summary$year==2019,"total count"] <- sum(grid_counts_1$ESBB_2019)
grid_count_summary[grid_count_summary$year==2020,"total count"] <- sum(grid_counts_1$ESBB_2020)
grid_count_summary[grid_count_summary$year==2021,"total count"] <- sum(grid_counts_1$ESBB_2021)
grid_count_summary[grid_count_summary$year==2022,"total count"] <-sum(grid_counts_2$ESBB_2022)
grid_count_summary[grid_count_summary$year==2023,"total count"] <-sum(grid_counts_2$ESBB_2023)
grid_count_summary[grid_count_summary$year==2025,"total count"] <-sum(grid_counts_2$ESBB_2025)

#What does the total butterfly count over time look like?
ggplot(data=grid_count_summary,aes(x=year, y=`total count`))+
  xlab("Year")+ylab("Total butterfly count")+
  geom_point(aes(x=year, y=`total count`))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

#Test if variations in total butterfly counts are normally distributed
#This test is needed to see which correlation test is appropriate
#for checking for significant correlations between total butterfly counts over time
#The test used is a Kolmogorov-Smirnov test
ks.test(grid_count_summary$`total count`,"pnorm")

#Does the total butterfly count vary significantly over time?
cor.test(grid_count_summary$year,grid_count_summary$`total count`,method="spearman")

#Make a unified data frame of El Segundo Blue Butterfly grid counts
tmp1 <- grid_counts_1[,c("ESBB_2019","Cell_ID")]
tmp1$year <- 2019
colnames(tmp1) <- c("count","Cell_ID","geometry","year")
tmp2 <- grid_counts_1[,c("ESBB_2020","Cell_ID")]
tmp2$year <- 2020
colnames(tmp2) <- c("count","Cell_ID","geometry","year")
tmp3 <- grid_counts_1[,c("ESBB_2021","Cell_ID")]
tmp3$year <- 2021
colnames(tmp3) <- c("count","Cell_ID","geometry","year")
tmp4 <- grid_counts_2[,c("ESBB_2022","Cell_ID")]
tmp4$year <- 2022
colnames(tmp4) <- c("count","Cell_ID","geometry","year")
tmp5 <- grid_counts_2[,c("ESBB_2023","Cell_ID")]
tmp5$year <- 2023
colnames(tmp5) <- c("count","Cell_ID","geometry","year")
tmp6 <- grid_counts_2[,c("ESBB_2025","Cell_ID")]
tmp6$year <- 2025
colnames(tmp6) <- c("count","Cell_ID","geometry","year")
grid_counts <- rbind(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
grid_counts <- grid_counts[!duplicated(grid_counts),]

#Map El Segundo Blue Butterfly grid counts for a selected year
pal <- colorNumeric(palette = "plasma", domain = grid_counts$count)
year_selected <- 2025
butterfly_plot <- grid_counts[grid_counts$year==year_selected,] |>
  filter(!st_is_empty(geometry)) |>
  st_make_valid() |>
  filter(st_geometry_type(geometry) %in% c("POLYGON", "MULTIPOLYGON")) |>
  st_cast("MULTIPOLYGON") |>
  st_transform(4326)
leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addPolygons(data=butterfly_plot,color = "transparent",fillColor = pal(st_drop_geometry(grid_counts$count))) |>
  addLegend(pal = pal, values =st_drop_geometry(grid_counts$count) , opacity = 1, title = paste(year_selected,"<br>El Segundo Blue<br>butterfly counts",sep=""),
            position = "bottomright")

#Calculate the percent count change per grid cell for 2019 - 2025
tmp1 <- (grid_counts[grid_counts$year==2019,c("count","Cell_ID")])
colnames(tmp1) <- c("count_2019","Cell_ID","geometry")
tmp2 <- (grid_counts[grid_counts$year==2025,c("count","Cell_ID")])
colnames(tmp2) <- c("count_2025","Cell_ID","geometry")
grid_change <- st_join(tmp1,tmp2)
grid_change$change <- ifelse(grid_change$count_2019==0,NA,(grid_change$count_2025-grid_change$count_2019)/grid_change$count_2019)

#Map changes in El Segundo Blue Butterfly grid counts 2019-2025
pal <- colorNumeric(palette = "plasma", domain = grid_change$change)
butterfly_plot <- grid_change |>
  filter(!st_is_empty(geometry)) |>
  st_make_valid() |>
  filter(st_geometry_type(geometry) %in% c("POLYGON", "MULTIPOLYGON")) |>
  st_cast("MULTIPOLYGON") |>
  st_transform(4326)
leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addPolygons(data=butterfly_plot,color = "transparent",fillColor = pal(st_drop_geometry(grid_change$change))) |>
  addLegend(pal = pal, values =st_drop_geometry(grid_change$change) , opacity = 1, title = "Percent change in<br>El Segundo Blue<br>butterfly counts<br>2019 - 2025",
            position = "bottomright")
