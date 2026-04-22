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

Lines 2-5: These are packages you'll want to install to run this script. There are a number of functions which R comes pre-installed with, but for many tasks you'll need to install packages to run other functions.

Lines 9-10: Define a path to your working directory, then tell your computer to set this as your working directory. A working directory is where you'll be running your script, and it's where you'll be telling your computer where to look for information.

Line 15: Read in our data. Once you do this you'll end up storing all of this information in an object called a data table. You can explore what's inside this data table by going to your console in the upper right-hand panel of your RStudio and clicking on the name of the data table. You'll see that the information is stored as a table of values, much like a spreadsheet.

Line 18: Remove entries from our data without count values. Since we'll ultimately be using count information to track the number of plastic bags found in beach cleanups we can only use rows in our data table which have count values.

Line 21: Create a plastic bag-specific data table which only contains rows from our data table where the category of garbage is designated as 'Plastic Bags' and there is a specified site location.

Line 24: Reformat the collection date column so that R knows to deal with it as an actual calendar date. By default, R will read in those data values as a string of characters and not know to deal with it as a date.

Line 27: Determine the earliest date in our plastic bag data set. This will be used to calculate the number of days from the start of data collection.

Line 31: Calculate the number of days from the earliest plastic bag collection date. Treat this value as a number.

Lines 34-37: Make a plot of the number of plastic bags collected versus the number of days since plastic bag data was being collected. Your plot should look something like this:

<img width="583" height="407" alt="image" src="https://github.com/user-attachments/assets/0571db6b-e889-4362-8b4c-8c335daff80d" />

Just looking at this plot it's difficult to tell if there's any trend, so we'll need to investigate further.

Line 43: Now we want to test if the number of plastic bags collected over time is significantly declining over time. However, to know which test to run we'll need to first determine if the distribution of our plastic bag counts varies normally or not. A [normal distribution](https://en.wikipedia.org/wiki/Normal_distribution), otherwise known as a bell-curve, occurs with a lot of different data sets and whether or not our data follows it will determine which statistical test is appropriate to use. Here we will use a [Kolmogorov-Smirnov test](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test). If the significance output, that is the value of p, is less that 0.05 then we can assume that our data are not normally distributed.

Line 50: Here we run a [Spearman correlation](https://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient) between the number of bags collected, and the number of days since our data collection began. We're using a Spearman correlation because it is used for data which are not normally distributed. Sometimes you will see this type of test being referred to as being non-parametric. A parametric test assumes normality with the data, and an example of such a test is known as a [Pearson correlation](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient).

Lines 56 - 61: Within our plastic bag data table we're going to create a new column whose value depends on whether or not plastic bags were collected and counted before or after the bag date of 1 January 2012. We're going to use this new columns to group our count data into before and after ban groups.

Lines 64 - 67: Here we'll plot the counts of plastic bags collected across all beach cleanups before and after the plastic bag ban. We'll use a [violin plot](https://en.wikipedia.org/wiki/Violin_plot) to visualize this data. We'll also be using the log of the count values to help stretch the plot which makes the plot easier to read. A violin plot is similar to a bar chart, except that the width of the bar depends on the number of measurements taken with that value. A violin plot gets wide at values which are frequently recorded, and thinner and ones which are not. Your plot should look something like this:

<img width="583" height="407" alt="image" src="https://github.com/user-attachments/assets/64e369d2-5cca-4279-a8b0-cbd6b74124d7" />

Line 72: Just looking at the plot it's hard to tell if there's a significant difference between our before and after picture. To really test this question we'll need to use some statistical tests. The first test we'll use is called a [Kruskal-Wallis test](https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_test). This is another non-parametric test, and it allows for us to test if the ban has a significant influence on the number of bags collected and counted. That is, if we treat the plastic bag ban as a binary variable does it have a significant effect on our count data? If the output significance value, p, is less than 0.05 we can say that it does.

Line 75: Now, can we test if the number of plastic bags collected and counted are significantly lower after the ban. This is done using a non-parametric test known as a [Wilcoxon rank-sum test](https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test). This is a test to see if, on average, a randomly selected value from one distribution of values tends to be larger or smaller than a second distribution. If the output significance value, p, is less than 0.05 we can say that the number of plastic bags found following the ban is significantly lower than before.

Line 78: Now we want to re-check if the plastic ban bag was effective by normalizing the number of plastic bags collected per beach cleanup event, and then seeing if there is a significant decline. We first do this by calculating the average number of bags collected per beach cleanup event.

Line 84: Now we want to test if the average number of plastic bags collected per beach cleanup over time is significantly declining over time. However, to know which test to run we'll need to first determine if the distribution of our average plastic bag counts varies normally or not. A [normal distribution](https://en.wikipedia.org/wiki/Normal_distribution), otherwise known as a bell-curve, occurs with a lot of different data sets and whether or not our data follows it will determine which statistical test is appropriate to use. Here we will use a [Kolmogorov-Smirnov test](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test). If the significance output, that is the value of p, is less that 0.05 then we can assume that our data are not normally distributed.

Lines 88 - 91: Here we'll plot the average counts of plastic bags collected per beach cleanup across all beach cleanups before and after the plastic bag ban. We'll use a [violin plot](https://en.wikipedia.org/wiki/Violin_plot) to visualize this data. We'll also be using the log of the count values to help stretch the plot which makes the plot easier to read. A violin plot is similar to a bar chart, except that the width of the bar depends on the number of measurements taken with that value. A violin plot gets wide at values which are frequently recorded, and thinner and ones which are not. Your plot should look something like this:

<img width="583" height="407" alt="image" src="https://github.com/user-attachments/assets/9037756b-3aa1-4dad-afa4-204e30314762" />

Line 95: Just looking at the plot it's hard to tell if there's a significant difference between our before and after picture. To really test this question we'll need to use some statistical tests. The first test we'll use is called a [Kruskal-Wallis test](https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_test). This is another non-parametric test, and it allows for us to test if the ban has a significant influence on the number of bags collected and counted. That is, if we treat the plastic bag ban as a binary variable does it have a significant effect on our count data? If the output significance value, p, is less than 0.05 we can say that it does.

Line 98: Now, can we test if the average number of plastic bags collected per beach cleanup are significantly lower after the ban. This is done using a non-parametric test known as a [Wilcoxon rank-sum test](https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test). This is a test to see if, on average, a randomly selected value from one distribution of values tends to be larger or smaller than a second distribution. If the output significance value, p, is less than 0.05 we can say that the average number of plastic bags per beach cleanup found following the ban is significantly lower than before.

###Extension: Now we want to check how the proportions of foodware and packaging shifting over time.

Lines 102 - 110: We want to again designate a numerical variable which is the number of days since data collection began.

Lines 112 - 114: Here we designate lists of trash categories associated with either packaging or foodware.

Lines 116 - 117: Here we designate a binary variable which state if a piece of trash is packaging or not.

Lines 119 - 123: Calculate the fraction of trash items which are packaging per beach cleanup collection event.

Lines 125 - 129: Plot the fraction of trash items which are packaging per beach cleanup collection event versus day. The plot should look like: 

<img width="583" height="407" alt="image" src="https://github.com/user-attachments/assets/f10489a1-2041-423c-a6fb-bb5d67bb1d11" />

Line 135: Now we want to test if the fraction of trash composed of packaging per beach cleanup over time is significantly declining over time. However, to know which test to run we'll need to first determine if the distribution of our fraction of trash composed of packaging per beach cleanup varies normally or not. A [normal distribution](https://en.wikipedia.org/wiki/Normal_distribution), otherwise known as a bell-curve, occurs with a lot of different data sets and whether or not our data follows it will determine which statistical test is appropriate to use. Here we will use a [Kolmogorov-Smirnov test](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test). If the significance output, that is the value of p, is less that 0.05 then we can assume that our data are not normally distributed.

Line 143: Here we run a [Spearman correlation](https://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient) between the fraction of trash composed of packaging per beach cleanup, and the number of days since our data collection began. We're using a Spearman correlation because it is used for data which are not normally distributed. Sometimes you will see this type of test being referred to as being non-parametric. A parametric test assumes normality with the data, and an example of such a test is known as a [Pearson correlation](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient).

Lines 145 - 146: Here we designate a binary variable which state if a piece of trash is foodware or not.

Lines 148 - 152: Calculate the fraction of trash items which are foodware per beach cleanup collection event.

Lines 154 - 158: Plot the fraction of trash items which are foodware per beach cleanup collection event versus day. The plot should look like: 

<img width="583" height="407" alt="image" src="https://github.com/user-attachments/assets/ebee80c6-4720-4df1-a161-b2933d7583d2" />

Line 164: Now we want to test if the fraction of trash composed of foodware per beach cleanup over time is significantly declining over time. However, to know which test to run we'll need to first determine if the distribution of our fraction of trash composed of foodware per beach cleanup varies normally or not. A [normal distribution](https://en.wikipedia.org/wiki/Normal_distribution), otherwise known as a bell-curve, occurs with a lot of different data sets and whether or not our data follows it will determine which statistical test is appropriate to use. Here we will use a [Kolmogorov-Smirnov test](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test). If the significance output, that is the value of p, is less that 0.05 then we can assume that our data are not normally distributed.

Line 172: Here we run a [Spearman correlation](https://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient) between the fraction of trash composed of foodware per beach cleanup, and the number of days since our data collection began. We're using a Spearman correlation because it is used for data which are not normally distributed. Sometimes you will see this type of test being referred to as being non-parametric. A parametric test assumes normality with the data, and an example of such a test is known as a [Pearson correlation](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient).

## Policy and economics

Once we are able to test our question with data the next step is to figure out the implications. Plastic bags can be banned, but what policies were involved, what were consequences. Here we can begin to investigate the political and economic context of our findings.

### Examples of effects of plastic bag bans
[A Survey on the Economic Effects of Los Angeles County’s Plastic Bag Ban](https://www.ncpathinktank.org/w18/st340/).

[Considerations, benefits and unintended consequences of banning plastic shopping bags for environmental sustainability: A systematic literature review](https://pmc.ncbi.nlm.nih.gov/articles/PMC8847762/).

[Plastic bag bans in the US reduced plastic bag use by billions, study finds](https://www.weforum.org/stories/2024/01/plastic-bag-bans-reduce-waste/).
### Examples of other plastic bag bans

In [California](https://www.cawrecycles.org/list-of-local-bag-bans).

## Art and science communication

Now that we have started to dig into the political and economic implications of our plastic bag ban findings, how can we communicate the results in a meaningful way?

### Examples of using art to communicate the science related to plastic bags

[UCLA’s Center for the Art of Performance Presents Robin Frohardt’s The Plastic Bag Store](https://hyperallergic.com/cap-ucla-presents-robin-frohardt-the-plastic-bag-store/)

Art project ideas from the [Because Turtles Eat Plastic Bags](https://becauseturtleseatplasticbags.com/resources/plastic-art-projects/) blog.

### Examples of public campaigns to communicate and/or advocate on the plastic bag ban issue

[How we finally banned plastic grocery bags in California](https://pirg.org/california/articles/how-we-finally-banned-plastic-grocery-bags-in-california/)

[Artivist Series - Dianna Cohen](https://womenmindthewater.com/artivist-series/artivist-series-dianna-cohen)
