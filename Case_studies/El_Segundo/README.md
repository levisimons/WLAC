# Central question: Are the ongoing efforts at habitat restoration helping the El Segundo Blue butterfly?

## How to consider the question.

Since 1986 the Los Angeles World Airports (LAWA) has been focused on restoring the habitat of the endangered El Segundo Blue butterfly. These efforts have involved demolishing infrastructure in the LAX dunes, removing invasive plants, and replanting Seacliff Buckwheat (Eriogonum parvifolium), the butterfly's sole food source and life-cycle host. As their habitat is being restored, are the number of butterflies significantly increasing?

## Data analysis

### How to install R and RStudio

In order to do our data analysis we will need to set up R and RStudio on our machines. R is a computer programming language commonly used in the biological sciences, but also in a number of projects involving statistical analysis. RStudio is a Integrated Development Environment (IDE), a tool for building and running code written in R. Both R and RStudio are free and open source, and can be run on a number of different operating systems. Instructions on how to install both R and RStudio can be found [here](https://rstudio-education.github.io/hopr/starting.html).

### Where to get data?

Data were provided by [LAWA](https://www.lawa.org/), which oversees operations of LAX airport. They are also involved in restoring dune habitat for the El Segundo Blue butterfly. In their restoration operations they have also be tracking the number of butterflies by counting them in a set of 30 m wide grid cells covering their 203 acre restoration area. These counts were made and mapped annually from 2019 to 2025, with the count map data available [here](Case_studies/El_Segundo/ESBB_GridTotals_2019_2021.zip) for the period 2019 - 2021, and [here](Case_studies/El_Segundo/Block_Counts_within_Grids.zip) for the period 2022 - 2025.

### Running analysis code

The script we can use to get started on analyzing our question statistically can be found [here](Case_studies/El_Segundo/El_Segundo.R). But what do the individual parts mean?

Line 1: Clear memory. This is good coding practice to make sure that there's nothing in memory before you run your current script.

Lines 2-7: These are packages you'll want to install to run this script. There are a number of functions which R comes pre-installed with, but for many tasks you'll need to install packages to run other functions.

Lines 10-11: Define a path to your working directory, then tell your computer to set this as your working directory. A working directory is where you'll be running your script, and it's where you'll be telling your computer where to look for information.
