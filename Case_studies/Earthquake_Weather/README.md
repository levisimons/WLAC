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

Lines 26-92: Create a function which takes in the latitude, longitude, and date of each earthquake, as well as specifying the units of temperature, and returns the average temperature associated with those earthquakes. Within this function the code functions as follows:

---

Lines 30-39: Process dates and return errors in case of missing data or improper formatting.

Lines 41-42: Define acceptable temperature units.

Lines 44-45: Set the url to put requests to for data using an API.

Lines 47-59: Run the API against the Open Meteo website to download temperature data using locations and dates for each earthquake.

Lines 61-68: Check if there's an error with connecting with the API. If there is then stop the function.

Lines 70-71: Store data returned from the Open Meteo website as a JSON object. JSON (JavaScript Object Notation) is a lightweight data-interchange format. 

Lines 73-75: Extract the temperature and time data from the JSON object.

Lines 77-78: Check if there's no temperature data in the JSON object. If there is then give an error message.

Lines 80-81: Calculate the average temperature from the temperature data.

Lines 83-91: Store date, location, and average temperature as a data frame to export from the function.

---

Lines 94-95: Initialize an empty list to store temperature results.

Lines 96-109: Run a loop to call the function to get temperature data via an API for each earthquake in your data. Within this loop the code functions as follows:

---

Lines 98-101: Supply coordinates in space and time for each earthquake.

Lines 103-104: Call the function to get the average temperature for each earthquake.

Lines 106-108: Print the status of the loop and store the average temperature results in a list.

---

Lines 111-112: Take the list of temperature values and convert them to a single data frame.

Lines 114-115: Join the average temperature data for each earthquake into the overall earthquake data.

Lines 117-121: Plot the magnitude of each earthquake against the average temperature associated with it. The plot should look like: 

<img width="583" height="407" alt="image" src="https://github.com/user-attachments/assets/75713f1c-7841-4f75-95ab-38158943108c" />

Lines 123-127: Now we want to test if the magnitude of earthquakes are significantly associated with average temperature. However, to know which test to run we'll need to first determine if the distribution of our earthquake magnitudes varies normally or not. A <a href="https://en.wikipedia.org/wiki/Normal_distribution">normal distribution</a>, otherwise known as a bell-curve, occurs with a lot of different data sets and whether or not our data follows it will determine which statistical test is appropriate to use. Here we will use a <a href="https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test">Kolmogorov-Smirnov test</a>. If the significance output, that is the value of p, is less that 0.05 then we can assume that our data are not normally distributed.

Lines 129-130: Here we run a Spearman correlation between the magnitude of earthquakes, and the average temperature associated with each of them. We're using a Spearman correlation because it is used for data which are not normally distributed. Sometimes you will see this type of test being referred to as being non-parametric. A parametric test assumes normality with the data, and an example of such a test is known as a Pearson correlation.
