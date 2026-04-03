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

Line 14: This command switches off the default s2 engine (developed by Google), which requires data to be in a projected coordinate system for operations like area, length, and distance if spherical geometry is disabled.This is done to disable spherical geometry in R and reduce the risk of some potential errors in working with spatial data.

Line 19: Read in mapped count data for 2019 - 2021.

Line 23: Read in mapped count data for 2022 - 2025.

Line 27: Make sure that the mapped count data for 2022 - 2025 has the same coordinate reference system as for the mapped count data collected in 2019 - 2021. A coordinate reference system is a framework that defines how coordinates relate to real locations on Earth, bridging the gap between the planet's 3D spherical surface and a 2D flat map.

Lines 32 - 40: Aggregate the total butterfly counts by year and check for aggregate trends over time. Make this table be seven rows long, one for each year. The two columns in the table are to store the year and the total count of butterflies from across the study area.

Line 43 - 46: Plot the total annual butterfly counts by year. You should get a plot like this:

<img width="583" height="407" alt="image" src="https://github.com/user-attachments/assets/1d58f7c1-8c6a-40e3-a7c6-08c93147c0f2" />

There definitely appears to be a trend over time, but let's test for how significant it is.

Line 52: Now we want to test if the number of butterflies counted over time is significantly increasing over time. However, to know which test to run we'll need to first determine if the distribution of our butterfly counts varies normally or not. A [normal distribution](https://en.wikipedia.org/wiki/Normal_distribution), otherwise known as a bell-curve, occurs with a lot of different data sets and whether or not our data follows it will determine which statistical test is appropriate to use. Here we will use a [Kolmogorov-Smirnov test](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test). If the significance output, that is the value of p, is less that 0.05 then we can assume that our data are not normally distributed.

Line 55: Here we run a [Spearman correlation](https://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient) between the number of butterflies counted, and the year of the count. We're using a Spearman correlation because it is used for data which are not normally distributed. Sometimes you will see this type of test being referred to as being non-parametric. A parametric test assumes normality with the data, and an example of such a test is known as a [Pearson correlation](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient).

Lines 58 - 76: Take all of the butterfly map count data and build a table out of it with the following columns: the number of butterflies counted per 30 m wide map cell, the ID number of that map cell, the location of that map cell, and the year.

Line 80: Create a color scale for visualizing the number of butterflies counted per map cell for a selected year.

Line 81: Select a year you want to map the butterfly count data by.

Lines 82 - 92: Create a map, with color scale, for the number of butterflies counted per map cell for a selected year. the map should look something like:

<img width="1162" height="810" alt="image" src="https://github.com/user-attachments/assets/a2f68181-7f15-424d-89c0-e3a2d32debad" />


Line 77: Remove any duplicated rows from this butterfly count table.
