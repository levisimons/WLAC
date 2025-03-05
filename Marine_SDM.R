rm(list=ls())
require(data.table)

#Set working directory
wd <- "" #Path to file on your own machine
setwd(wd)

#Intake species occurrence data from GBIF
#Available from querying GBIF: https://doi.org/10.15468/dl.j46jh7
input_species_data <- fread(input="macrocystis_pyrifera.csv",sep="\t")
