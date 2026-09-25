#!/usr/bin/env Rscript
#Contains functions to generate values with multiple probability density functions
library(insight)#print_color function

#Function for multiple random normal distributions
multirnorm <- function(ns, mean, sd, w=c('eq')){
	if ('eq' %in% w){
		if (length(mean) == length(sd)){
			num <- length(mean)
			w <- rep(1/num, num)
		}else {
			print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!UNEQUAL NUMBER OF ARGUMENTS!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bred')
        		break
		}
	}else {
		w <- as.numeric(w)
		if (!all.equal(sum(w),1)){
			print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!!WEIGHTS DO NOT EQUAL ONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bred')
        		break
		}
		if ((length(mean) == length(sd)) & (length(mean) == length(w))){
			num <- length(mean)
		}else {
			print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!UNEQUAL NUMBER OF ARGUMENTS!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bred')
        		break
		}
	}

	#Sort weights and make sure the means and sd are properly sorted as well
	templab <- order(w, decreasing  = TRUE)
	w <- w[templab]
	mean <- as.numeric(mean[templab])
	sd <- as.numeric(sd[templab])

	#Iterate over number of samples
	values <- c()
	for (i in 1:ns){
		set.seed(NULL)#makes sure parameters end up in new bin each time 
		r <- runif(1, min=0, max=1)

		#Iterate over possible pdfs
		check <- w[1]
		for (j in 1:num){
			if (r < check){
				val <- rnorm(1, mean=mean[j], sd=sd[j])
				break
			}else {
				if (j == num){
					val <- rnorm(1, mean=mean[j], sd=sd[j])
					break
				}else {
					check <- check + w[j+1]
				}
			}
		}
		values <- c(values,val)
	}	
	return(values)
}#end of multirnorm

#Function for multiple random uniform distributions
multirunif <- function(ns, min, max, w=c('eq')){
	if ('eq' %in% w){
		if (length(min) == length(max)){
			num <- length(min)
			w <- rep(1/num, num)
		}else {
			print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!UNEQUAL NUMBER OF ARGUMENTS!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bred')
        		break
		}
	}else {
		w <- as.numeric(w)
		if (!all.equal(sum(w),1)){
			print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!!WEIGHTS DO NOT EQUAL ONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bred')
        		break
		}
		if ((length(min) == length(max)) & (length(min) == length(w))){
			num <- length(min)
		}else {
			print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!UNEQUAL NUMBER OF ARGUMENTS!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bred')
        		break
		}
	}

	#Sort weights and make sure the means and sd are properly sorted as well
	templab <- order(w, decreasing  = TRUE)
	w <- w[templab]
	min <- as.numeric(min[templab])
	max <- as.numeric(max[templab])

	#Iterate over number of samples
	values <- c()
	for (i in 1:ns){
		set.seed(NULL)#makes sure parameters end up in new bin each time 
		r <- runif(1, min=0, max=1)
	
		#Iterate over possible pdfs
		check <- w[1]
		for (j in 1:num){
			if (r < check){
				val <- runif(1, min=min[j], max=max[j])
				break
			}else {
				if (j == num){
					val <- runif(1, min=min[j], max=max[j])
					break
				}else {
					check <- check + w[j+1]
				}
			}
		}
		values <- c(values,val)
	}	
	return(values)
}#end of multirunif
