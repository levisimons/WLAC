rm(list=ls())
require(data.table)
require(lubridate)
require(dplyr)
require(ggplot2)

#Set working directory
#Put in the file path to where you want to work on your computer
wd <- ""
setwd(wd)

#Read in Heal The Bay data
#Data can be downloaded and unzipped from https://github.com/levisimons/WLAC/blob/main/Case_studies/Heal_the_bay/HTB.csv.zip
#Make sure data is in your working directory
HTB_input <- fread(input="HTB.csv",sep=",")

#Remove entries without plastic bag counts
HTB_input <- HTB_input[!is.na(HTB_input$count),]

#Subset data to only include plastic bag entries
HTB_plastic_bags <- HTB_input[HTB_input$subcategory=='Plastic Bags' & HTB_input$Site!="",]

#Set date column to a standard format
HTB_plastic_bags$`Collected Date` <- as.Date(format(mdy_hm(HTB_plastic_bags$`Collected Date`), "%m/%d/%Y"),format="%m/%d/%Y")

#Find the earliest sampling date
min_date <- min(HTB_plastic_bags$`Collected Date`)

#Create a variable which is the number of days following the earliest sampling date
#Make it numeric for plotting purposes.
HTB_plastic_bags$day <- as.numeric(HTB_plastic_bags$`Collected Date`-min_date)

#Plot plastic bag counts versus day
ggplot(data=HTB_plastic_bags,aes(x=day, y=count))+
  xlab("Days from 21 February 2001")+ylab("Plastic bag counts")+
  geom_point(aes(x=day, y=count))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

#Test if variations in plastic bag counts are normally distributed
#This test is needed to see which correlation test is appropriate
#for checking for significant correlations between plastic bag counts over time
#The test used is a Kolmogorov-Smirnov test
ks.test(HTB_plastic_bags$count,"pnorm")

#The output for the Kolmogorov-Smirnov test is a p value less than 0.05.
#This indicates that the distribution of plastic bag counts is not normally distributed.
#This then means that a non-parameteric test, such as Spearman, will be needed to 
#check for significant correlations between plastic bag counts over time
#Test, using a Spearman correlation, if there are any significant trends over time
#for the number of plastic bags counted.
cor.test(HTB_plastic_bags$count,HTB_plastic_bags$day,method="spearman")

#Was there a significant drop in plastic bags after the bag ban on 1 January 2012?

#Create a grouping variable for when samples were collected
HTB_plastic_bags <- HTB_plastic_bags %>%
  dplyr::mutate(`Collected Date` = as.Date(`Collected Date`, format = "%m/%d/%Y"),
    ban_status = case_when(
    `Collected Date` < as.Date("01/01/2012",format="%m/%d/%Y") ~ "before plastic bag ban",
    `Collected Date` >= as.Date("01/01/2012",format="%m/%d/%Y") ~ "after plastic bag ban"
  ))

#Violin plot of plastic bag counts whether or not they were collected before or after the plastic bag ban
#Use a log-scale on plastic bag counts to help visualize the distributions
ggplot(HTB_plastic_bags, aes(x=ban_status, y=count) )+
  xlab("Plastic ban status")+ylab("log(Plastic bag counts)")+
  geom_violin(aes(x=ban_status, y=log10(count)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Test if ban status is significant in influencing plastic bag counts.
#Use a Kruskal-Wallis test.
kruskal.test(HTB_plastic_bags$count,HTB_plastic_bags$ban_status)

#Test if plastic bag counts are significantly lower following plastic bag ban.
wilcox.test(HTB_plastic_bags[HTB_plastic_bags$ban_status=="after plastic bag ban",]$count,HTB_plastic_bags[HTB_plastic_bags$ban_status=="before plastic bag ban",]$count,alternative="less")

#Calculate the average number of plastic bags per collection event.
HTB_plastic_bags[, average_count := mean(count), by = .(Site, day)]

#Test if variations in the average plastic bag counts per collection event are normally distributed
#This test is needed to see which correlation test is appropriate
#for checking for significant correlations between plastic bag counts over time
#The test used is a Kolmogorov-Smirnov test
ks.test(HTB_plastic_bags$average_count,"pnorm")

#Violin plot of the average plastic bag counts per collection event on whether or not they were collected before or after the plastic bag ban
#Use a log-scale on plastic bag counts to help visualize the distributions
ggplot(HTB_plastic_bags, aes(x=ban_status, y=average_count) )+
  xlab("Plastic ban status")+ylab("log(Average plastic bag counts\nper collection event)")+
  geom_violin(aes(x=ban_status, y=log10(count)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Test if ban status is significant in influencing the collection event average of plastic bag counts.
#Use a Kruskal-Wallis test.
kruskal.test(HTB_plastic_bags$average_count,HTB_plastic_bags$ban_status)

#Test if collection event average of plastic bag counts are significantly lower following plastic bag ban.
wilcox.test(HTB_plastic_bags[HTB_plastic_bags$ban_status=="after plastic bag ban",]$average_count,HTB_plastic_bags[HTB_plastic_bags$ban_status=="before plastic bag ban",]$average_count,alternative="less")

#How are the proportions of foodware and packaging shifting over time?

#Set date column to a standard format
HTB_input$`Collected Date` <- as.Date(format(mdy_hm(HTB_input$`Collected Date`), "%m/%d/%Y"),format="%m/%d/%Y")

#Find the earliest sampling date
min_date <- min(HTB_input$`Collected Date`)

#Create a variable which is the number of days following the earliest sampling date
#Make it numeric for plotting purposes.
HTB_input$day <- as.numeric(HTB_input$`Collected Date`-min_date)

#Define lists of foodware and packaging terms
packaging_list <- c("Plastic Beverage Bottles","Plastic Snack Bags / Wrappers","Plastic Bottle Caps / Rings","6-Pack Rings","Liquid Bottles / Large Containers","Tobacco Package","Glass Bottles")
foodware_list <- c("Plastic Cups / Lids","Plastic Utensils","Plastic Plates","Plastic Straws / Stirrers","Foam Take-Out Containers","Foam Cups / Plates","Paper Containers","Paper Cups","Paper Plates")

#Designate a packaging category
HTB_input$packaging <- ifelse(HTB_input$subcategory %in% packaging_list,1,0)

#Calculate the fraction of trash items which are packaging per collection event
HTB_input <- HTB_input %>%
  group_by(Site, day) %>%
  mutate(packaging_fraction = sum(count[packaging == 1])/sum(count[packaging %in% c(1,0)])) %>%
  ungroup()

#Plot the fraction of trash items which are packaging per collection event versus day
ggplot(data=HTB_input,aes(x=day, y=packaging_fraction))+
  xlab("Days from 21 February 2001")+ylab("Trash fraction as packaging\nper collection event")+
  geom_point(aes(x=day, y=packaging_fraction))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

#Test if variations in plastic bag counts are normally distributed
#This test is needed to see which correlation test is appropriate
#for checking for significant correlations between plastic bag counts over time
#The test used is a Kolmogorov-Smirnov test
ks.test(HTB_input$packaging_fraction,"pnorm")

#The output for the Kolmogorov-Smirnov test is a p value less than 0.05.
#This indicates that the distribution of the fraction of packaging trash is not normally distributed.
#This then means that a non-parameteric test, such as Spearman, will be needed to 
#check for significant correlations between the fraction of packaging trash over time
#Test, using a Spearman correlation, if there are any significant trends over time
#for the number of plastic bags counted.
cor.test(HTB_input$packaging_fraction,HTB_input$day,method="spearman")

#Designate a foodware category
HTB_input$foodware <- ifelse(HTB_input$subcategory %in% foodware_list,1,0)

#Calculate the fraction of trash items which are foodware per collection event
HTB_input <- HTB_input %>%
  group_by(Site, day) %>%
  mutate(foodware_fraction = sum(count[foodware == 1])/sum(count[foodware %in% c(1,0)])) %>%
  ungroup()

#Plot the fraction of trash items which are foodware per collection event versus day
ggplot(data=HTB_input,aes(x=day, y=foodware_fraction))+
  xlab("Days from 21 February 2001")+ylab("Trash fraction as foodware\nper collection event")+
  geom_point(aes(x=day, y=foodware_fraction))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

#Test if variations in the fraction of trash items which are foodware per collection event are normally distributed
#This test is needed to see which correlation test is appropriate
#for checking for significant correlations between the fraction of trash items which are foodware per collection event over time
#The test used is a Kolmogorov-Smirnov test
ks.test(HTB_input$foodware_fraction,"pnorm")

#The output for the Kolmogorov-Smirnov test is a p value less than 0.05.
#This indicates that the distribution of the fraction of packaging trash is not normally distributed.
#This then means that a non-parameteric test, such as Spearman, will be needed to 
#check for significant correlations between the fraction of trash items which are foodware per collection event over time
#Test, using a Spearman correlation, if there are any significant trends over time
#for the fraction of trash items which are foodware per collection event.
cor.test(HTB_input$foodware_fraction,HTB_input$day,method="spearman")

#Read in beach site grouping information to filter out sites in state beaches.
site_groups <- fread(input="CleanupSites_Groups.csv",sep=",")

#Designate state beaches
state_beaches <- c("Manhattan State Beach","Dockweiler State Beach","Santa Monica State Beach","Will Rogers State Beach","Topanga State Beach","Malibu Lagoon State Beach / Surfrider Beach","Point Dume State Beach","Robert H. Meyer Memorial State Beach","Leo Carrillo State Beach")

#Identify cleanup sites within state beaches
state_beach_sites <- unique(site_groups[site_groups$`Parent Beach Category` %in% state_beaches,]$`Cleanup Site Name`)

#Get cigarette butt count data for cleanup sites within state beaches
HTB_cigarettes <- HTB_input[HTB_input$Site %in% state_beach_sites & HTB_input$subcategory=="Cigarette Butts",]

#Determine if samples are taken before or after 1 January 2020, the date of Senate Bill 8
#the law which banned smoking in state beaches.
HTB_cigarettes <- HTB_cigarettes %>%
  dplyr::mutate(`Collected Date` = as.Date(`Collected Date`, format = "%m/%d/%Y"),
                ban_status = case_when(
                  `Collected Date` < as.Date("01/01/2020",format="%m/%d/%Y") ~ "before smoking ban",
                  `Collected Date` >= as.Date("01/01/2020",format="%m/%d/%Y") ~ "after smoking ban"
                ))

#Calculate the average number of cigarette butts per collection event.
HTB_cigarettes <- as.data.table(HTB_cigarettes)
HTB_cigarettes[, average_count := mean(count), by = .(Site, day)]

#Violin plot of the average cigarette counts per collection event on whether or not they were collected before or after the plastic bag ban
#Use a log-scale on cigarette counts to help visualize the distributions
ggplot(HTB_cigarettes, aes(x=ban_status, y=average_count) )+
  xlab("Smoking ban status")+ylab("log(Average cigarette counts\nper collection event)")+
  geom_violin(aes(x=ban_status, y=log10(average_count)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Test if ban status is significant in influencing the collection event average of cigarette butt counts.
#Use a Kruskal-Wallis test.
kruskal.test(HTB_cigarettes$average_count,HTB_cigarettes$ban_status)

#Test if collection event average of cigarette counts are significantly lower following smoking ban.
wilcox.test(HTB_cigarettes[HTB_cigarettes$ban_status=="after smoking ban",]$average_count,HTB_cigarettes[HTB_cigarettes$ban_status=="before smoking ban",]$average_count,alternative="less")
