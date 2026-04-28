# Central question: Is there such a thing as earthquake weather?

## How to consider the question.

The USGS keeps global records of earthquake magnitudes, locations, and dates. By combining this information with records of average temperatures recorded at each of those earthquake locations and dates we can test for potential relationships between the two.

## Data analysis

### How to install R and RStudio

In order to do our data analysis we will need to set up R and RStudio on our machines. R is a computer programming language commonly used in the biological sciences, but also in a number of projects involving statistical analysis. RStudio is a Integrated Development Environment (IDE), a tool for building and running code written in R. Both R and RStudio are free and open source, and can be run on a number of different operating systems. Instructions on how to install both R and RStudio can be found <a href="https://rstudio-education.github.io/hopr/starting.html">here</a>.

### Where to get data

USGS earthquake data can be downloaded via their site <a href="https://earthquake.usgs.gov/earthquakes/search/">here</a>. For this case study we downloaded the global data set for 2025 for all earthquakes with a magnitude of at least 4.0. The result of this query can also be found <a href="https://github.com/levisimons/WLAC/blob/main/Case_studies/Earthquake_Weather/USGS_2025.csv">here</a>.

Temperature data for each location and date will be obtained programmatically using an Application Programming Interface (API) to connect to <a href="https://archive-api.open-meteo.com/v1/archive">Open Meteo site</a>. Using their API will allow for us to directly connect to the site and download data using specific inputs.

### Running the analysis code

Line 1: Clear memory. This is good coding practice to make sure that there's nothing in memory before you run your current script.

Lines 2-8: These are packages you'll want to install to run this script. There are a number of functions which R comes pre-installed with, but for many tasks you'll need to install packages to run other functions.

Lines 10-13: Define a path to your working directory, then tell your computer to set this as your working directory. A working directory is where you'll be running your script, and it's where you'll be telling your computer where to look for information.

Lines 15-17: Read in downloaded earthquake data.

Lines 19-21: Filter earthquake data to only contain entries where the magnitude of the earthquakes exceeds some set threshold.

Lines 23-24: Reformat the collection date column so that R knows to deal with it as an actual calendar date. By default, R will read in those data values as a string of characters and not know to deal with it as a date.

---

Lines 26-92: Create a function which takes in the latitude, longitude, and date of each earthquake, as well as specifying the units of temperature, and returns the average temperature associated with those earthquakes. Within this function the code functions as follows:
