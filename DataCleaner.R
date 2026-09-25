#!/usr/bin/env Rscript
#Above line allows code to be run using ./DataCleaner.R in terminal

#Libraries and what they are used for commented next to them
library(dplyr)#as_tibble and many other dataframe manipulation shortcuts
library(data.table)#setnames function
library(insight)#print_color function

#Collect data
fcidf <- read.csv(paste0('realdata/FCI-ARK.csv'))
fcitbbl <- as_tibble(fcidf)

#Define columns of interest
fcipost <- paste0('AQ',1:30,'.y')
fcipostgender <- 'FinalGender.y'
fcipre <- paste0('AQ',1:30,'.x')
fcipregender <- 'FinalGender.x'

#New column names
fciItem <- paste0('Item',1:30)

#Printing out the full tibble so one can see column names and data types
print_color('============================================================================\n','bold')
print_color('==============================Cleaned Data Set==============================\n','bold')
print_color('============================================================================\n','bold')

print_color('====================================FCI=====================================\n','bgreen')
fcipostdata <- fcitbbl %>%
	select(all_of(c(fcipost,fcipostgender)))
setnames(fcipostdata, old = c(fcipost,fcipostgender), new = c(fciItem,'Gender'))
print(fcipostdata)
write.csv(fcipostdata, 'realdata/FCI-post.csv', row.names = FALSE)

fcipredata <- fcitbbl %>%
	select(all_of(c(fcipre,fcipregender)))
setnames(fcipredata, old = c(fcipre,fcipregender), new = c(fciItem,'Gender'))
print(fcipredata)
write.csv(fcipredata, 'realdata/FCI-pre.csv', row.names = FALSE)

