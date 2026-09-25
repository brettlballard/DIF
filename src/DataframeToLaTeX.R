#!/usr/bin/env Rscript
#Contains function to sort item names that is used in several places

#Libraries and what they are used for commented next to them
library(dplyr)#as_tibble and many other dataframe manipulation shortcuts
library(insight)#print_color function

#Converting a dataframe to LaTeX table format
dftoLaTeX <- function(data, suppress=NA, filename='Test'){
	df <- data
	fl <- file(paste0(filename,'.txt'),'w')
	print(df)
	print(colnames(df))
	print(rownames(df))
	#Suppressing table values if desired 
	if (!is.na(suppress)){
		print_color(paste0('=====================Suppressing Table Values Below ',suppress,'===================\n'),'byellow')
		df[abs(df) < suppress] <- NA
		df <- df %>%
			mutate(across(where(is.numeric), \(x) round(x,3))) %>%
			print()
	}

	#Convert dataframe entries besides title into LaTeX format 
	for (i in 1:nrow(df)){
		row <- df[i,]
		row[is.na(row)] <- ''
		writeLines(paste0(rownames(df)[i], ' & ', paste(row, collapse=' & '), ' \\\\','\n'), con = fl)
		writeLines(paste0('\\hline\n'), con = fl)
	}
	close(fl)
}
