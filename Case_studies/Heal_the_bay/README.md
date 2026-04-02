# Central question: Was the Los Angeles county plastic bag ban effective?


# How to consider the question.

By 1 January 2012, Los Angeles county had put in a ban for single-use plastic bans. Here we will consider various facets for addressing the question of how effective this policy has been.

## Data analysis

### How to install R and RStudio

In order to do our data analysis we will need to set up R and RStudio on our machines. R is a computer programming language commonly used in the biological sciences, but also in a number of projects involving statistical analysis. RStudio is a Integrated Development Environment (IDE), a tool for building and running code written in R. Both R and RStudio are free and open source, and can be run on a number of different operating systems. Instructions on how to install both R and RStudio can be found [here](https://rstudio-education.github.io/hopr/starting.html).

### Where to get data?

Data was provided by [Heal The Bay](https://healthebay.org/), a non-profit dedicated to the health of the Santa Monica bay. The particular type of data consisted of [beach cleanup](https://healthebay.org/take-part/) records collected by volunteers over the period 2001 - 2025. These records contained information on the type of trash collected, as well as when and where it was collected. A zipped version of this file can be found [here](Case_studies/Heal_the_bay/HTB.csv.zip).

### Running analysis code

The script we can use to get started on analyzing our question statistically can be found [here](Case_studies/Heal_the_bay/HTB.R). But what do the individual parts mean?

Line 1: Clear memory. This is good coding practice to make sure that there's nothing in memory before you run your current script.

Lines 2-4: These are packages you'll want to install to run this script. There are a number of functions which R comes pre-installed with, but for many tasks you'll need to install packages to run other functions.

Lines 8-9: Define a path to your working directory, then tell your computer to set this as your working directory. A working directory is where you'll be running your script, and it's where you'll be telling your computer where to look for information.

Line 14: Read in our data. Once you do this you'll end up storing all of this information in an object called a data table. You can explore what's inside this data table by going to your console in the upper right-hand panel of your RStudio and clicking on the name of the data table. You'll see that the information is stored as a table of values, much like a spreadsheet.



## Policy and economics

Once we are able to test our question with data the next step is to figure out the implications. Plastic bags can be banned, but what policies were involved, what were consequences. Here we can begin to investigate the political and economic context of our findings.

### Examples of effects of plastic bag bans
[A Survey on the Economic Effects of Los Angeles County’s Plastic Bag Ban](https://www.ncpathinktank.org/w18/st340/).

[Considerations, benefits and unintended consequences of banning plastic shopping bags for environmental sustainability: A systematic literature review](https://pmc.ncbi.nlm.nih.gov/articles/PMC8847762/).

[Plastic bag bans in the US reduced plastic bag use by billions, study finds](https://www.weforum.org/stories/2024/01/plastic-bag-bans-reduce-waste/).
### Examples of other plastic bag bans

In [California](https://www.cawrecycles.org/list-of-local-bag-bans).

## Art and science communication
