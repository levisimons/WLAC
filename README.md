# WLAC
Ecological code repository for West Los Angeles College.

The following instructions are for getting the data needed to run the species distribution modeling script [Marine_SDM.R](https://github.com/levisimons/WLAC/blob/main/Marine_SDM.R).

# To install R and R Studio
Following the tutorial here: [R and R Studio tutorial](https://alexd106.github.io/intro2R/setup.html).

# Getting biological data
## Go to [https://www.gbif.org/occurrence/search](https://www.gbif.org/occurrence/search)
## Set query filters
### Occurrence status: present
### Scientific name: Macrocystis pyrifera .  Note: Giant kelp, test example
### Location: including coordinates
### Administrative area: California
### Time: 2010 - 2020
## A download link will be generated.  The output file will be a tab-delimited CSV file.
### Example output file for Macrocystis pyrifera: [macrocystis_pyrifera.csv](https://github.com/levisimons/WLAC/blob/main/macrocystis_pyrifera.csv)

# Defining the study area
## Got to download the shapefile for [3-Nautical Mile Coastal Boundary, California, 2001](https://earthworks.stanford.edu/catalog/stanford-sg211gq3741).
### Relevant shapefile data is initially in the directory sg211gq3741/data_EPSG_4326 , move it to the same directory as the [script](https://github.com/levisimons/WLAC/blob/main/Marine_SDM.R).
### Within [Marine_SDM.R](https://github.com/levisimons/WLAC/blob/main/Marine_SDM.R) clip the coastal boundary zone of California to the extent of the Southnern California Bight. Northwestern point of study area is Point Conception (34.4486° N, 120.4716° W), southeastern point is Tijuana Beach Promenade (32.5340° N, 117.1235° W)
