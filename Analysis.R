#!/usr/bin/env Rscript
#Above line allows code to be run using ./Analysis.R in terminal

#Curious about runtime
start <- Sys.time()

#Libraries and what they are used for commented next to them
library(dplyr)#as_tibble and many other dataframe manipulation shortcuts
library(data.table)#setnames function
library(insight)#print_color function
library(argparser)#anything parser related
library(reshape2)#melt function
library(difR)#MH Test
library(ggplot2)#plot related
ggshapes <- c(0:14,32:127)
library(geomtextpath)#geom_text_segment
library(ggrepel)#geom_text_repel
library(treemapify)#treemap plots
library(ggmosaic)#mosaic plots
library(cowplot)#combining plots

#Adding argument parsers so that I can vary the simulated data from the command line
parser <- arg_parser('Options for varying the simulated data generated')
parser <- add_argument(parser, "--name", help = 'name of output when in flex/run mode if other name desired',nargs='*',default='TEST')
parser <- add_argument(parser, "--DIFtype", help = 'uniform, nonuniform, or mixed: default is uniform',nargs='*',default='uniform')
parser <- add_argument(parser, "--DIFsimtype", help = 'simple or variable: default is simple',nargs='*',default='simple')

#Test arguments
parser <- add_argument(parser, "--test", help = 'MH, binMH, SIBU, SIBNU, Lord, RajuU, RajuS: default is MH',nargs='*',default='MH')
parser <- add_argument(parser, '--nbins', help = 'number of bins for the Binned MH Test: default is 5',nargs='*',default=5)
parser <- add_argument(parser, '--niter', help = 'number of iterations for the Purified MH Tests: default is 10',nargs='*',default=10)
parser <- add_argument(parser, "--padj", help = 'options are none, BH, or holm: default is holm',nargs='*',default='holm')
parser <- add_argument(parser, "--purify", help = 'purify when TRUE: default is FALSE',nargs='*',default=FALSE)

#additional arguments
parser <- add_argument(parser, "--nitems", help = 'number of items when in flex mode: format input as begin,end,increment',nargs='*',default=c(20,20,0))
parser <- add_argument(parser, "--ns", help = 'number of students when in flex mode: format input as begin,end,increment',nargs='*',default=c(1000,1000,0))
parser <- add_argument(parser, "--nfiles", help = 'number of files in each folder: default is 1',nargs='*',default=c(1,1,0))
arg <- parse_args(parser)

#Turning multiple input arguments into vectors
nitems <- seq(from = arg$nitems[1], to = arg$nitems[2], by = arg$nitems[3])
numst <- seq(from = arg$ns[1], to = arg$ns[2], by = arg$ns[3])
nfiles <- seq(from = arg$nfiles[1], to = arg$nfiles[2], by = arg$nfiles[3])

if (arg$test == 'MH'){
	print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!!!!!!RUNNING MH ANALYSIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bcyan')
}else if (arg$test == 'binMH'){
	print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!!!!RUNNING binMH ANALYSIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bcyan')
}else if (arg$test == 'SIBU'){
	print_color(paste0('!!!!!!!!!!!!!!!!!!!!!RUNNING Uniform SIBTEST ANALYSIS!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bcyan')
}else if (arg$test == 'SIBNU'){
	print_color(paste0('!!!!!!!!!!!!!!!!!!!RUNNING Nonuniform SIBTEST ANALYSIS!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bcyan')
}else if (arg$test == 'Lord'){
	print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!!!!RUNNING Lord\'s ANALYSIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bcyan')
}else if (arg$test == 'RajuU'){
	print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!RUNNING Raju\'s Unsigned ANALYSIS!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bcyan')
}else if (arg$test == 'RajuS'){
	print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!RUNNING Raju\'s Signed ANALYSIS!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bcyan')
}else {
	print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!INCORRECT INPUT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bred')
	break
}

if (arg$test %in% c('MH','binMH','SIBU','SIBNU','Lord','RajuU','RajuS')){
	if (arg$purify){
		print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!!!!RUNNING Purified ANALYSIS!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bcyan')
	}else {
		print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!!!RUNNING Unpurified ANALYSIS!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bcyan')
	}

	if (arg$padj == 'none'){
		print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!!!!!NO p-value adjustment!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bcyan')
	}else if (arg$padj == 'BH'){
		print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!!!!!BH p-value adjustment!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bcyan')
	}else if (arg$padj == 'holm'){
		print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!!!!Holm p-value adjustment!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bcyan')
	}else {
		print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!INCORRECT INPUT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bred')
		break
	}
}

if (!dir.exists(paste0('analysisout/',arg$name,'/',arg$test))){dir.create(paste0('analysisout/',arg$name,'/',arg$test), recursive = TRUE)}

###############################################################################################################
##################################################FUNCTIONS####################################################
###############################################################################################################
#Function that removes NA from a vector so that it's length is 0
rmna <- function(vec){
	vec <- vec[!is.na(vec)]
	return(vec)
}

#Binned Mantel-Haenzsel Test
#Function that bins Mantel-Haenszel Test
binMH <- function(data, items, group, purify=FALSE, p.adjust.method='none'){
	padj <- p.adjust.method
	nbins <- as.numeric(arg$nbins)
	nitems <- length(items)
	if (purify){
		total <- 'AnchorTotal'
	}else {
		total <- 'Total'
	}
	#Automating bin selection based on user input
	bins <- unname(round(quantile(data[[total]], probs = seq(0,1, by = 1/nbins)[2:nbins])))
	print_color(paste0('Number of bins being used: ',nbins,'\n'),'bgreen')
	print_color(paste0('Scores being used to separate bins: ',paste(bins, collapse = ', '),'\n'),'bgreen')
	for (i in 1:nbins){
		if (i == 1){
			data[[total]] <- ifelse((data[[total]] < bins[i]),i,data[[total]])
		}else if (i == nbins){
			data[[total]] <- ifelse((data[[total]] >= bins[(i-1)]),i,data[[total]])
		}else {
			data[[total]] <- ifelse((data[[total]] >= bins[(i-1)] & data[[total]] < bins[i]),i,data[[total]])
		}
	}	
	difvalues <- c()
	for (i in colnames(data[,items])){
		data[[i]] <- factor(data[[i]], levels = c(1,0))
		data[[group]] <- factor(data[[group]], levels = c(0,1))
		alphanum <- c()
		alphaden <- c()
		testnum <- c()
		testden <- c()
		#Calculations for each bin
		for (b in 1:nbins){
			tab <- table(data[data[[total]] == b, c(group,i)])
			ni <- tab[1,1]+tab[1,2]+tab[2,1]+tab[2,2]
			alphanum <- append(alphanum, (tab[1,1]*tab[2,2]/ni))
			alphaden <- append(alphaden, (tab[1,2]*tab[2,1]/ni))
			testnum <- append(testnum, (tab[1,1] - (tab[1,1]+tab[1,2])*(tab[1,1]+tab[2,1])/ni))
			testden <- append(testden, ((tab[1,1]+tab[1,2])*(tab[2,1]+tab[2,2])*(tab[1,1]+tab[2,1])*((tab[1,2]+tab[2,2])/((ni**2)*(ni-1)))))
		}
		#Calculations for item
		alphaMH <- sum(alphanum, na.rm=TRUE)/sum(alphaden, na.rm=TRUE)
		deltaMH <- -2.35*log(alphaMH)
		if (abs(deltaMH) < 1){
			ETSDeltaCode <- 'A'
		}else if (abs(deltaMH) < 1.5){
			ETSDeltaCode <- 'B'
		}else if (abs(deltaMH) >= 1.5){
			ETSDeltaCode <- 'C'
		}

		MHtest <- ((abs(sum(testnum, na.rm=TRUE)) - .5)**2)/sum(testden, na.rm=TRUE)
		pvalue <- 1 - pchisq(MHtest, df=1)
		itemvalues <- c(round(alphaMH,4),round(deltaMH,4),ETSDeltaCode,round(MHtest,4),round(pvalue,4))
		difvalues <- append(difvalues, itemvalues)
	}
	rownam <- colnames(data[,items])
	colnam <- c('alphaMH','deltaMH','ETS.Code','MH.Stat', 'P.value')
	binmhmat <- matrix(difvalues, nrow = nitems, byrow = TRUE, dimnames = list(rownam,colnam))
	binmhdf <- as.data.frame(binmhmat)
	colnames(binmhdf) <- colnam
	if (padj == 'holm' | padj == 'BH'){
		binmhdf$Adj.Alpha.05 <- rep(0, times = nitems)
	}
	binmhdf$Sig <- rep(' ', times = nitems)
	index <- order(as.numeric(binmhdf$P.value), decreasing = FALSE)
	rank <- 1

	#Making either Holm or BH adjustments for multiple comparisons
	if (padj == 'holm'){
		#Making Holm's adjuments to control family-wise Type I error rate
		for (i in index){
			pval <- as.numeric(binmhdf$P.value[i])
			adjalpha05 <- .05/(nitems+1-rank)
			binmhdf$Adj.Alpha.05[i] <- round(adjalpha05,4)
			adjalpha1 <- .1/(nitems+1-rank)
			adjalpha01 <- .01/(nitems+1-rank)
			adjalpha001 <- .001/(nitems+1-rank)
			if (pval < adjalpha001){
				binmhdf$Sig[i] <- '***'
			}else if (pval < adjalpha01){
				binmhdf$Sig[i] <- '**'
			}else if (pval < adjalpha05){
				binmhdf$Sig[i] <- '*'
			}else if (pval < adjalpha1){
				binmhdf$Sig[i] <- '.'
			}
			rank <- rank + 1
		}
	}else if (padj == 'BH'){
		#Making Benjamini-Hochburg adjuments to control the false discovery rate
		for (i in index){
			pval <- as.numeric(binmhdf$P.value[i])
			adjalpha05 <- (.05*rank)/nitems
			binmhdf$Adj.Alpha.05[i] <- round(adjalpha05,4) 
			adjalpha1 <- (.1*rank)/nitems
			adjalpha01 <- (.01*rank)/nitems
			adjalpha001 <- (.001*rank)/nitems
			if (pval < adjalpha001){
				binmhdf$Sig[i] <- '***'
			}else if (pval < adjalpha01){
				binmhdf$Sig[i] <- '**'
			}else if (pval < adjalpha05){
				binmhdf$Sig[i] <- '*'
			}else if (pval < adjalpha1){
				binmhdf$Sig[i] <- '.'
			}
			rank <- rank + 1
		}
	}else if (padj == 'none'){
		#Making no adjuments to control the false discovery rate or family-wise Type I error rate
		for (i in index){
			pval <- as.numeric(binmhdf$P.value[i])
			if (pval < .001){
				binmhdf$Sig[i] <- '***'
			}else if (pval < .01){
				binmhdf$Sig[i] <- '**'
			}else if (pval < .05){
				binmhdf$Sig[i] <- '*'
			}else if (pval < .1){
				binmhdf$Sig[i] <- '.'
			}
		}
	}else {
		print_color('Value of p value adjustment not allowed!','bred')
	}
	sigcodes <- 'Significance Codes For Global Alphas: 0 \'***\' .001 \'**\' .01 \'*\' .05 \'.\' .1 \' \' 1\n'
	
	#Finding dif items through significance tests
	if (padj == 'none'){
		binmhdif <- binmhdf[as.numeric(binmhdf$P.value) < .05,]
	}else {
		binmhdif <- binmhdf[as.numeric(binmhdf$P.value) < as.numeric(binmhdf$Adj.Alpha.05),]
	}
	binmhdifitems <- row.names(binmhdif)
	binmhdifitemsout <- ifelse(length(row.names(binmhdif)) == 0,NA,paste(gsub('s','',row.names(binmhdif)),collapse=', '))
	
	#Finding large effect size dif items
	binmheffdif <- binmhdf[binmhdf$ETS.Code == 'C',]
	binmheffdifitems <- row.names(binmheffdif)
	binmheffdifitemsout <- ifelse(length(row.names(binmheffdif)) == 0,NA,paste(gsub('s','',row.names(binmheffdif)),collapse=', '))
	
	#Ouputting stuff
	binmhout <- list('binMHdf' = binmhdf, 'sig' = sigcodes, 'DIFitems' = binmhdifitemsout, 'EffDIFitems' = binmheffdifitemsout)
	return(binmhout)
}

#Function that outputs binMH results  
printbinMH <- function(binMHobj){
	print(binMHobj$binMHdf)
	cat(binMHobj$sig)
	cat(paste0('Items detected as DIF items: ', binMHobj$DIFitems,'\n'))
	return(invisible(NULL))
}

#Function that removes DIF positive items from matching criteria to be used in purification
rmitems <- function(data, difitems){
	rmitems <- strsplit(difitems, ', ')[[1]]
	print_color('Removing Items: ', 'bold')
	#Removing DIFitems from Total
	for (i in rmitems){
		rmitem <- i
		for (j in rmitem){
			print_color(paste0(j, ' '), 'bred')
			data$Total <- data$Total-data[[j]]
		}	
	}
	cat('\n')
	return(data$Total)	

}

#Function that purifies Binned Mantel-Haenszel Test
#Purification removes DIF items from matching criteria calculation iteratively
pbinMH <- function(data, items, group, niter, p.adjust.method='none', criteria='eff'){
	print_color(paste0('Number of maximum iterations in the binning process: ',niter,'\n'),'bgreen')
	initbinmh <- binMH(data = data, items = items, group = group, p.adjust.method = p.adjust.method)
	if (criteria == 'eff'){
		crit <- 'EffDIFitems'
	}else if (criteria == 'pval'){
		crit <- 'DIFitems'
	}
	#Output intial if no DIF items
	if (length(rmna(initbinmh[[crit]])) == 0){
		pbinmhout <- list('binMHdf' = initbinmh$binMHdf, 'sig' = initbinmh$sig, 'DIFitems' = initbinmh$DIFitems, 'EffDIFitems' = initbinmh$EffDIFitems)
		return(pbinmhout)
	#Output intial if all are DIF items
	}else if (length(rmna(initbinmh[[crit]])) == length(items)){
		pbinmhout <- list('binMHdf' = initbinmh$binMHdf, 'sig' = initbinmh$sig, 'DIFitems' = initbinmh$DIFitems, 'EffDIFitems' = initbinmh$EffDIFitems)
		return(pbinmhout)
	}else {
		difitems <- initbinmh[[crit]]
		data$AnchorTotal <- rmitems(data = data, difitems = difitems)
		#print(data$AnchorTotal)
		for (count in 2:niter){
			tempbinmh <- binMH(data = data, items = items, group = group, purify = TRUE, p.adjust.method = p.adjust.method)
			#Checking to see if no items dif positive once dif items removed from matching criteria
			if (length(rmna(tempbinmh[[crit]])) == 0){
				pbinmhout <- list('binMHdf' = tempbinmh$binMHdf, 'sig' = tempbinmh$sig, 'DIFitems' = tempbinmh$DIFitems, 'EffDIFitems' = tempbinmh$EffDIFitems)
				print_color(paste0('Purification ended after ',count,' iterations\n'),'bgreen')
				return(pbinmhout)
				break
			#Checking to see if same items dif positive once dif items removed from matching criteria
			}else if(tempbinmh[[crit]] == difitems){
				pbinmhout <- list('binMHdf' = tempbinmh$binMHdf, 'sig' = tempbinmh$sig, 'DIFitems' = tempbinmh$DIFitems, 'EffDIFitems' = tempbinmh$EffDIFitems)
				print_color(paste0('Purification ended after ',count,' iterations\n'),'bgreen')
				return(pbinmhout)
				break
			#Time to remove the new dif items and try again
			}else {
				if (count == niter){
					pbinmhout <- list('binMHdf' = tempbinmh$binMHdf, 'sig' = tempbinmh$sig, 'DIFitems' = tempbinmh$DIFitems, 'EffDIFitems' = tempbinmh$EffDIFitems)
					print_color(paste0('Purification did not converge after maximum number of iterations and results from final iteration are given\n'),'bred')
					return(pbinmhout)
				}else {
					difitems <- tempbinmh[[crit]]
					data$AnchorTotal <- rmitems(data = data, difitems = difitems)
				}
			}
		}
	}
}



###############################################################################################################
###################################################ANALYSIS####################################################
###############################################################################################################
#Analyzing simulated data 
namevec <- c()
itemsvec <- c()
stvec <- c()
fvec <- c()
typevec <- c()
simtypevec <- c()
magvec <- c()
focalstpropvec <- c()
consistvec <- c()
numdifitemsvec <- c()
numrefdifitemsvec <- c()
testvec <- c()
palphavec <- c()
pbetavec <- c()
pprevvec <- c()
if (grepl('MH',arg$test) | grepl('SIB',arg$test)){
	Balphavec <- c()
	Bbetavec <- c()
	Bprevvec <- c()
	Calphavec <- c()
	Cbetavec <- c()
	Cprevvec <- c()
	if (grepl('MH',arg$test)){
		ptpdeltavec <- c()
		pfpdeltavec <- c()
		pfndeltavec <- c()
		ptndeltavec <- c()
	}
}

#Collecting all data for test to compare
for (nit in nitems){
	for (nst in numst){
		for (f in nfiles){
			print_color(paste0('==============================================================================\n'),'bgreen')
			print_color(paste0('===============================DIF Analysis===================================\n'),'bgreen')
			print_color(paste0('==============================================================================\n'),'bgreen')
			print_color(paste0('====================================',f,'=======================================\n'),'bcyan')
			#Saving variables
			namevec <- c(namevec, arg$name)
			itemsvec <- c(itemsvec, nit)
			stvec <- c(stvec, nst)
			fvec <- c(fvec, f)
			typevec <- c(typevec, arg$DIFtype)
			simtypevec <- c(simtypevec, arg$DIFsimtype)
			testvec <- c(testvec, arg$test)
			
			#Retrieving and saving generators
			gen <- read.csv(paste0('simdata/',arg$name,'/',nit,'items/',nst,'students/',arg$DIFtype,'/',arg$DIFsimtype,'/',paste0(arg$name,f),'-Generators.csv'))
			magvec <- c(magvec, mean(gen$DIF.Magnitude))
			consistvec <- c(consistvec, mean(gen$DIF.Consistency))
			
			#Building dataframe for metric calculations
			DIFItems <- read.csv(paste0('simdata/',arg$name,'/',nit,'items/',nst,'students/',arg$DIFtype,'/',arg$DIFsimtype,'/',paste0(arg$name,f),'-Items.csv'))
			DIFItems <- DIFItems[,c('Items','DIF.Item')]
			origdifitems <- DIFItems$DIF.Item
			numdifitemsvec <- c(numdifitemsvec,length(origdifitems[origdifitems != 0]))
			numrefdifitemsvec <- c(numrefdifitemsvec,length(origdifitems[origdifitems == 2]))
			DIFItems$DIF.Item <- ifelse(DIFItems$DIF.Item > 0,1,0)
			DIFItems$DIF.Item <- factor(DIFItems$DIF.Item, levels = c(1,0))
		
			#Collecting data
			data <- read.csv(paste0('simdata/',arg$name,'/',nit,'items/',nst,'students/',arg$DIFtype,'/',arg$DIFsimtype,'/',paste0(arg$name,f),'-Data.csv'))
			Items <- paste0('Item',1:nit)
			data$Total <- apply(data[,Items],1,sum)
			focalstpropvec <- c(focalstpropvec, mean(data$Group))

			#Running DIF test
			if (arg$test == 'MH'){
				diftest <- difMH(Data = data[,c(Items,'Group')], group = 'Group', focal.name = 1, purify = arg$purify, p.adjust.method = arg$padj)
				delMH <- -2.35*log(diftest$alphaMH)
			}else if (arg$test == 'binMH'){
				if (arg$purify){
					bMHtest <- pbinMH(data = data, items = Items, group = 'Group', niter = arg$niter, p.adjust.method = arg$padj, criteria = 'pval')
				}else {
					bMHtest <- binMH(data = data, items = Items, group = 'Group', p.adjust.method = arg$padj)
				}
				diftest <- bMHtest$binMHdf
				delMH <- as.numeric(diftest$deltaMH)
				if (is.na(bMHtest$DIFitems)){
					bMHDIFitems <- c('NO DIF ITEMS') 
				}else {
					bMHDIFitems <- strsplit(bMHtest$DIFitems,', ')[[1]]
				}
			}else if (arg$test == 'SIBU'){
				testfail <- tryCatch(
					expr = {
						diftest <- difSIBTEST(Data = data[,c(Items,'Group')], group = 'Group', focal.name = 1, purify = arg$purify, p.adjust.method = arg$padj)
						betaSIB <- diftest$Beta
					},
					error = function(e){
						print_color(paste0('===========================UNIFORM SIBTEST FAILED===============================\n'),'bred')
						print(e)
					}
				)
			}else if (arg$test == 'SIBNU'){
				testfail <- tryCatch(
					expr = {
						diftest <- difSIBTEST(Data = data[,c(Items,'Group')], group = 'Group', focal.name = 1, type = 'nudif', purify = arg$purify, p.adjust.method = arg$padj)
						betaSIB <- diftest$Beta
					},
					error = function(e){
						print_color(paste0('==========================NONUNIFORM SIBTEST FAILED=============================\n'),'bred')
						print(e)
					}
				)
			}else if (arg$test == 'Lord'){
				testfail <- tryCatch(
					expr = {
						diftest <- difLord(Data = data[,c(Items,'Group')], group = 'Group', focal.name = 1, model = '2PL', purify = arg$purify, p.adjust.method = arg$padj)
					},
					error = function(e){
						print_color(paste0('=============================LORD\'S TEST FAILED================================\n'),'bred')
						print(e)
					}
				)
			}else if (arg$test == 'RajuU'){
				testfail <- tryCatch(
					expr = {
						diftest <- difRaju(Data = data[,c(Items,'Group')], group = 'Group', focal.name = 1, model = '2PL', signed = FALSE, purify = arg$purify, p.adjust.method = arg$padj)
					},
					error = function(e){
						print_color(paste0('========================Raju\'S Unsigned TEST FAILED============================\n'),'bred')
						print(e)
					}
				)
			}else if (arg$test == 'RajuS'){
				testfail <- tryCatch(
					expr = {
						diftest <- difRaju(Data = data[,c(Items,'Group')], group = 'Group', focal.name = 1, model = '2PL', signed = TRUE, purify = arg$purify, p.adjust.method = arg$padj)
					},
					error = function(e){
						print_color(paste0('=========================Raju\'S Signed TEST FAILED=============================\n'),'bred')
						print(e)
					}
				)
			}

			#Set fail condition to be false for all MH and SIBTEST type tests
			if (grepl('MH',arg$test)){
				testfail <- c('Pass') 
			}

			if (!any(class(testfail) == 'error')){
				print(diftest)
			}

			#Determine effect sizes for MH and SIBTEST
			if (grepl('MH',arg$test)){
				Bdif <- which(abs(delMH) > 1)
				Cdif <- which(abs(delMH) > 1.5)
			}
			if (!any(class(testfail) == 'error') & grepl('SIB',arg$test)){
				Bdif <- which(abs(betaSIB) >= .059)
				Cdif <- which(abs(betaSIB) >= .088)
			}

			#Build dataframe for metrics to be calculated
			if (!any(class(testfail) == 'error')){
				if (arg$test %in% c('MH','SIBU','SIBNU')){
					DIFItems$ptest <- ifelse(DIFItems$Items %in% diftest$names[diftest$DIFitems],1,0)
					DIFItems$Befftest <- ifelse(DIFItems$Items %in% diftest$names[Bdif],1,0)
					DIFItems$Cefftest <- ifelse(DIFItems$Items %in% diftest$names[Cdif],1,0)
				}else if (arg$test == 'binMH'){
					DIFItems$ptest <- ifelse(DIFItems$Items %in% strsplit(bMHDIFitems,', ')[[1]],1,0)
					DIFItems$Befftest <- ifelse(DIFItems$Items %in% rownames(diftest)[Bdif],1,0)
					DIFItems$Cefftest <- ifelse(DIFItems$Items %in% rownames(diftest)[Cdif],1,0)
				}else if (arg$test %in% c('Lord','RajuU','RajuS')){
					DIFItems$ptest <- ifelse(DIFItems$Items %in% diftest$names[diftest$DIFitems],1,0)
				}
				
				DIFItems$ptest <- factor(DIFItems$ptest, levels = c(1,0))
				if (grepl('MH',arg$test) | grepl('SIB',arg$test)){
					DIFItems$Befftest <- factor(DIFItems$Befftest, levels = c(1,0))
					DIFItems$Cefftest <- factor(DIFItems$Cefftest, levels = c(1,0))
				}
				print(DIFItems)
			}

			#Collect item values for comparisons
			if (grepl('MH',arg$test)){
				#Finding mean delta for items that fall in each cell
				tpitems <- which(DIFItems$DIF.Item == 1 & DIFItems$ptest == 1)
				temp <- delMH[tpitems]
				ptpdeltavec <- c(ptpdeltavec,mean(temp[is.finite(temp)]))
				fpitems <- which(DIFItems$DIF.Item == 0 & DIFItems$ptest == 1)
				temp <- delMH[fpitems]
				pfpdeltavec <- c(pfpdeltavec,mean(temp[is.finite(temp)]))
				fnitems <- which(DIFItems$DIF.Item == 1 & DIFItems$ptest == 0)
				temp <- delMH[fnitems]
				pfndeltavec <- c(pfndeltavec,mean(temp[is.finite(temp)]))
				tnitems <- which(DIFItems$DIF.Item == 0 & DIFItems$ptest == 0)
				temp <- delMH[tnitems]
				ptndeltavec <- c(ptndeltavec,mean(temp[is.finite(temp)]))
			}
				
			if (!any(class(testfail) == 'error')){
				#Calculating accuracy, true and false positive rate, and precision
				tab <- table(data = DIFItems[,c('DIF.Item','ptest')])
				print(tab)
				alpha <- tab[2,1]/(tab[2,1]+tab[2,2])
				beta <- tab[1,2]/(tab[1,2]+tab[1,1])
				prev <- (tab[1,1]+tab[1,2])/nit
				palphavec <- c(palphavec,alpha)
				pbetavec <- c(pbetavec,beta)
				pprevvec <- c(pprevvec,prev)
			}else {
				palphavec <- c(palphavec,NA)
				pbetavec <- c(pbetavec,NA)
				pprevvec <- c(pprevvec,NA)
			}

			if (grepl('MH',arg$test) | (grepl('SIB',arg$test) & !any(class(testfail) == 'error'))){
				#Calculating accuracy, true and false positive rate, and precision
				tab <- table(data = DIFItems[,c('DIF.Item','Befftest')])
				print(tab)
				alpha <- tab[2,1]/(tab[2,1]+tab[2,2])
				beta <- tab[1,2]/(tab[1,2]+tab[1,1])
				prev <- (tab[1,1]+tab[1,2])/nit
				Balphavec <- c(Balphavec,alpha)
				Bbetavec <- c(Bbetavec,beta)
				Bprevvec <- c(Bprevvec,prev)
				
				#Calculating accuracy, true and false positive rate, and precision
				tab <- table(data = DIFItems[,c('DIF.Item','Cefftest')])
				print(tab)
				alpha <- tab[2,1]/(tab[2,1]+tab[2,2])
				beta <- tab[1,2]/(tab[1,2]+tab[1,1])
				prev <- (tab[1,1]+tab[1,2])/nit
				Calphavec <- c(Calphavec,alpha)
				Cbetavec <- c(Cbetavec,beta)
				Cprevvec <- c(Cprevvec,prev)
			}else if ((grepl('SIB',arg$test) & any(class(testfail) == 'error'))){
				Balphavec <- c(Balphavec,NA)
				Bbetavec <- c(Bbetavec,NA)
				Bprevvec <- c(Bprevvec,NA)
				Calphavec <- c(Calphavec,NA)
				Cbetavec <- c(Cbetavec,NA)
				Cprevvec <- c(Cprevvec,NA)
			}
		}#end of nfiles loop
	}#end of nst loop
}#end of nitems loop

print_color(paste0('==============================================================================\n'),'bviolet')
print_color(paste0('=================================DIF Plots====================================\n'),'bviolet')
print_color(paste0('==============================================================================\n'),'bviolet')
if (grepl('MH',arg$test) | grepl('SIB',arg$test)){
	if (grepl('MH',arg$test)){
		difdata <- data.frame(Name = namevec, Num.Items = itemsvec, Num.Students = stvec, Num.File = fvec, DIF.Type = typevec, Sim.Type = simtypevec, Mag = magvec, Prop.Focal.St = focalstpropvec, Consistency = consistvec, Num.DIF.Items = numdifitemsvec, Num.Ref.DIF.Items = numrefdifitemsvec, Test = testvec, p.Alpha = palphavec, p.Beta = pbetavec, p.Prevalence = pprevvec, Beff.Alpha = Balphavec, Beff.Beta = Bbetavec, Beff.Prevalence = Bprevvec, Ceff.Alpha = Calphavec, Ceff.Beta = Cbetavec, Ceff.Prevalence = Cprevvec, TP.Delta = ptpdeltavec, TN.Delta = ptndeltavec, FP.Delta = pfpdeltavec, FN.Delta = pfndeltavec)
	}else {
		difdata <- data.frame(Name = namevec, Num.Items = itemsvec, Num.Students = stvec, Num.File = fvec, DIF.Type = typevec, Sim.Type = simtypevec, Mag = magvec, Prop.Focal.St = focalstpropvec, Consistency = consistvec, Num.DIF.Items = numdifitemsvec, Num.Ref.DIF.Items = numrefdifitemsvec, Test = testvec, p.Alpha = palphavec, p.Beta = pbetavec, p.Prevalence = pprevvec, Beff.Alpha = Balphavec, Beff.Beta = Bbetavec, Beff.Prevalence = Bprevvec, Ceff.Alpha = Calphavec, Ceff.Beta = Cbetavec, Ceff.Prevalence = Cprevvec)
	}
	difdata <- difdata %>%
		mutate(p.Sensitivity = 1 - p.Beta) %>%
		mutate(p.Specificity = 1 - p.Alpha) %>%
		mutate(p.Accuracy = 1 - (p.Prevalence*p.Beta + (1 - p.Prevalence)*p.Alpha)) %>%
		mutate(p.PPV = (p.Prevalence*(1 - p.Beta)/(p.Prevalence*(1 - p.Beta) + (1 - p.Prevalence)*p.Alpha))) %>%
		mutate(p.NPV = ((1 - p.Prevalence)*(1 - p.Alpha)/((1 - p.Prevalence)*(1 - p.Alpha) + p.Prevalence*p.Beta))) %>%
		mutate(Beff.Sensitivity = 1 - Beff.Beta) %>%
		mutate(Beff.Specificity = 1 - Beff.Alpha) %>%
		mutate(Beff.Accuracy = 1 - (Beff.Prevalence*Beff.Beta + (1 - Beff.Prevalence)*Beff.Alpha)) %>%
		mutate(Beff.PPV = (Beff.Prevalence*(1 - Beff.Beta)/(Beff.Prevalence*(1 - Beff.Beta) + (1 - Beff.Prevalence)*Beff.Alpha))) %>%
		mutate(Beff.NPV = ((1 - Beff.Prevalence)*(1 - Beff.Alpha)/((1 - Beff.Prevalence)*(1 - Beff.Alpha) + Beff.Prevalence*Beff.Beta))) %>%
		mutate(Ceff.Sensitivity = 1 - Ceff.Beta) %>%
		mutate(Ceff.Specificity = 1 - Ceff.Alpha) %>%
		mutate(Ceff.Accuracy = 1 - (Ceff.Prevalence*Ceff.Beta + (1 - Ceff.Prevalence)*Ceff.Alpha)) %>%
		mutate(Ceff.PPV = (Ceff.Prevalence*(1 - Ceff.Beta)/(Ceff.Prevalence*(1 - Ceff.Beta) + (1 - Ceff.Prevalence)*Ceff.Alpha))) %>%
		mutate(Ceff.NPV = ((1 - Ceff.Prevalence)*(1 - Ceff.Alpha)/((1 - Ceff.Prevalence)*(1 - Ceff.Alpha) + Ceff.Prevalence*Ceff.Beta)))
}else if (arg$test %in% c('Lord','RajuU','RajuS')){
	difdata <- data.frame(Name = namevec, Num.Items = itemsvec, Num.Students = stvec, Num.File = fvec, DIF.Type = typevec, Sim.Type = simtypevec, Mag = magvec, Prop.Focal.St = focalstpropvec, Consistency = consistvec, Num.DIF.Items = numdifitemsvec, Num.Ref.DIF.Items = numrefdifitemsvec, Test = testvec, Alpha = palphavec, Beta = pbetavec, Prevalence = pprevvec)
	difdata <- difdata %>%
		mutate(Sensitivity = 1 - Beta) %>%
		mutate(Specificity = 1 - Alpha) %>%
		mutate(Accuracy = 1 - (Prevalence*Beta + (1 - Prevalence)*Alpha)) %>%
		mutate(PPV = (Prevalence*(1 - Beta)/(Prevalence*(1 - Beta) + (1 - Prevalence)*Alpha))) %>%
		mutate(NPV = ((1 - Prevalence)*(1 - Alpha)/((1 - Prevalence)*(1 - Alpha) + Prevalence*Beta)))
}
print(as_tibble(difdata))
print(head(as.data.frame(difdata),20))

if (arg$test == 'MH'){
	title <- ggdraw()+draw_label(paste0('Mantel-Haenzel Test'), fontface='bold', hjust=.5)
}else if (arg$test == 'binMH'){
	title <- ggdraw()+draw_label(paste0('Binned Mantel-Haenzel Test with ',arg$nbins,' bins'), fontface='bold', hjust=.5)
}else if (arg$test == 'SIBU'){
	title <- ggdraw()+draw_label(paste0('Uniform SIBTEST'), fontface='bold', hjust=.5)
}else if (arg$test == 'SIBNU'){
	title <- ggdraw()+draw_label(paste0('Nonuniform SIBTEST'), fontface='bold', hjust=.5)
}else if (arg$test == 'Lord'){
	title <- ggdraw()+draw_label(paste0('Lord\'s Test'), fontface='bold', hjust=.5)
}else if (arg$test == 'RajuU'){
	title <- ggdraw()+draw_label(paste0('Unsigned Raju\'s Test'), fontface='bold', hjust=.5)
}else if (arg$test == 'RajuS'){
	title <- ggdraw()+draw_label(paste0('Signed Raju\'s Test'), fontface='bold', hjust=.5)
}

if (grepl('MH',arg$test)){
	#######################################
	#Plot delta means against DIF magnitude
	#######################################
	pldata <- melt(difdata, id=c('Name','Num.Items','Num.Students','DIF.Type','Sim.Type','Test','Mag','Num.DIF.Items'), measure=c('TP.Delta','FP.Delta'))
	pldata <- pldata %>% 
		rename(Measure = variable) %>%
		group_by(Mag,Num.DIF.Items,Measure) %>%
		summarize(value.mn = mean(value, na.rm=TRUE)) %>%
		print()

	plmin <- min(pldata$value.mn, na.rm=TRUE)
	plmax <- max(pldata$value.mn, na.rm=TRUE)

#Plot different number of dif items separately
	m2plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 2,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(a) 2 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Delta')+coord_cartesian(ylim=c(plmin,plmax))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
	m4plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 4,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(b) 4 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Delta')+coord_cartesian(ylim=c(plmin,plmax))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
	m6plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 6,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(c) 6 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Delta')+coord_cartesian(ylim=c(plmin,plmax))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
	m8plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 8,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(d) 8 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Delta')+coord_cartesian(ylim=c(plmin,plmax))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
	m10plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 10,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(f) 10 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Delta')+coord_cartesian(ylim=c(plmin,plmax))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')

	mplot <- plot_grid(m2plot+theme(legend.position='none'),m4plot+theme(legend.position='none'),m6plot+theme(legend.position='none'),m8plot+theme(legend.position='none'),m10plot+theme(legend.position='none'), ncol=2)
	legend <- get_legend(m2plot+guides(color = guide_legend(nrow=1))+theme(legend.position='bottom'))
	mplot <- plot_grid(title, mplot, legend, ncol=1, rel_heights=c(.05,.9,.05))
	if (arg$test == 'binMH'){
		ggsave(file=paste0('Delta-Plots-',arg$purify,'-',arg$padj,'-',arg$nbins,'.pdf'), path=paste0('analysisout/',arg$name,'/',arg$test,'/'), mplot, width=7.5, height=7.5, units='in')
	}else {
		ggsave(file=paste0('Delta-Plots-',arg$purify,'-',arg$padj,'.pdf'), path=paste0('analysisout/',arg$name,'/',arg$test,'/'), mplot, width=7.5, height=7.5, units='in')
	}
}

##################################
#Plot alphas against DIF magnitude
##################################
if (grepl('MH',arg$test) | grepl('SIB',arg$test)){
	pldata <- melt(difdata, id=c('Name','Num.Items','Num.Students','DIF.Type','Sim.Type','Test','Mag','Num.DIF.Items'), measure=c('p.Alpha','Beff.Alpha','Ceff.Alpha'))
}else {
	pldata <- melt(difdata, id=c('Name','Num.Items','Num.Students','DIF.Type','Sim.Type','Test','Mag','Num.DIF.Items'), measure=c('Alpha'))
}
pldata <- pldata %>% 
	rename(Measure = variable) %>%
	group_by(Mag,Num.DIF.Items,Measure) %>%
	summarize(value.mn = mean(value, na.rm=TRUE)) %>%
	print()

#Plot different number of dif items separately
m2plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 2,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(a) 2 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Type I Error Rate', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m4plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 4,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(b) 4 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Type I Error Rate', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m6plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 6,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(c) 6 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Type I Error Rate', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m8plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 8,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(d) 8 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Type I Error Rate', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m10plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 10,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(f) 10 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Type I Error Rate', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')

mplot <- plot_grid(m2plot+theme(legend.position='none'),m4plot+theme(legend.position='none'),m6plot+theme(legend.position='none'),m8plot+theme(legend.position='none'),m10plot+theme(legend.position='none'), ncol=2)
legend <- get_legend(m2plot+guides(color = guide_legend(nrow=1))+theme(legend.position='bottom'))
mplot <- plot_grid(title, mplot, legend, ncol=1, rel_heights=c(.05,.9,.05))
if (arg$test == 'binMH'){
	ggsave(file=paste0('TypeI-Plots-',arg$purify,'-',arg$padj,'-',arg$nbins,'.pdf'), path=paste0('analysisout/',arg$name,'/',arg$test,'/'), mplot, width=7.5, height=7.5, units='in')
}else {
	ggsave(file=paste0('TypeI-Plots-',arg$purify,'-',arg$padj,'.pdf'), path=paste0('analysisout/',arg$name,'/',arg$test,'/'), mplot, width=7.5, height=7.5, units='in')
}

#################################
#Plot betas against DIF magnitude
#################################
if (grepl('MH',arg$test) | grepl('SIB',arg$test)){
	pldata <- melt(difdata, id=c('Name','Num.Items','Num.Students','DIF.Type','Sim.Type','Test','Mag','Num.DIF.Items'), measure=c('p.Beta','Beff.Beta','Ceff.Beta'))
}else {
	pldata <- melt(difdata, id=c('Name','Num.Items','Num.Students','DIF.Type','Sim.Type','Test','Mag','Num.DIF.Items'), measure=c('Beta'))
}
pldata <- pldata %>% 
	rename(Measure = variable) %>%
	group_by(Mag,Num.DIF.Items,Measure) %>%
	summarize(value.mn = mean(value, na.rm=TRUE)) %>%
	print()

#Plot different number of dif items separately
m2plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 2,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(a) 2 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Type II Error Rate', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m4plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 4,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(b) 4 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Type II Error Rate', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m6plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 6,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(c) 6 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Type II Error Rate', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m8plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 8,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(d) 8 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Type II Error Rate', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m10plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 10,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(f) 10 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Type II Error Rate', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')

mplot <- plot_grid(m2plot+theme(legend.position='none'),m4plot+theme(legend.position='none'),m6plot+theme(legend.position='none'),m8plot+theme(legend.position='none'),m10plot+theme(legend.position='none'), ncol=2)
legend <- get_legend(m2plot+guides(color = guide_legend(nrow=1))+theme(legend.position='bottom'))
mplot <- plot_grid(title, mplot, legend, ncol=1, rel_heights=c(.05,.9,.05))
if (arg$test == 'binMH'){
	ggsave(file=paste0('TypeII-Plots-',arg$purify,'-',arg$padj,'-',arg$nbins,'.pdf'), path=paste0('analysisout/',arg$name,'/',arg$test,'/'), mplot, width=7.5, height=7.5, units='in')
}else {
	ggsave(file=paste0('TypeII-Plots-',arg$purify,'-',arg$padj,'.pdf'), path=paste0('analysisout/',arg$name,'/',arg$test,'/'), mplot, width=7.5, height=7.5, units='in')
}

#########################################
#Plot sensitivities against DIF magnitude
#########################################
if (grepl('MH',arg$test) | grepl('SIB',arg$test)){
	pldata <- melt(difdata, id=c('Name','Num.Items','Num.Students','DIF.Type','Sim.Type','Test','Mag','Num.DIF.Items'), measure=c('p.Sensitivity','Beff.Sensitivity','Ceff.Sensitivity'))
}else {
	pldata <- melt(difdata, id=c('Name','Num.Items','Num.Students','DIF.Type','Sim.Type','Test','Mag','Num.DIF.Items'), measure=c('Sensitivity'))
}
pldata <- pldata %>% 
	rename(Measure = variable) %>%
	group_by(Mag,Num.DIF.Items,Measure) %>%
	summarize(value.mn = mean(value, na.rm=TRUE)) %>%
	print()

#Plot different number of dif items separately
m2plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 2,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(a) 2 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Sensitivity', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m4plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 4,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(b) 4 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Sensitivity', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m6plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 6,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(c) 6 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Sensitivity', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m8plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 8,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(d) 8 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Sensitivity', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m10plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 10,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(f) 10 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Sensitivity', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')

mplot <- plot_grid(m2plot+theme(legend.position='none'),m4plot+theme(legend.position='none'),m6plot+theme(legend.position='none'),m8plot+theme(legend.position='none'),m10plot+theme(legend.position='none'), ncol=2)
legend <- get_legend(m2plot+guides(color = guide_legend(nrow=1))+theme(legend.position='bottom'))
mplot <- plot_grid(title, mplot, legend, ncol=1, rel_heights=c(.05,.9,.05))
if (arg$test == 'binMH'){
	ggsave(file=paste0('Sensitivity-Plots-',arg$purify,'-',arg$padj,'-',arg$nbins,'.pdf'), path=paste0('analysisout/',arg$name,'/',arg$test,'/'), mplot, width=7.5, height=7.5, units='in')
}else {
	ggsave(file=paste0('Sensitivity-Plots-',arg$purify,'-',arg$padj,'.pdf'), path=paste0('analysisout/',arg$name,'/',arg$test,'/'), mplot, width=7.5, height=7.5, units='in')
}

#########################################
#Plot specificities against DIF magnitude
#########################################
if (grepl('MH',arg$test) | grepl('SIB',arg$test)){
	pldata <- melt(difdata, id=c('Name','Num.Items','Num.Students','DIF.Type','Sim.Type','Test','Mag','Num.DIF.Items'), measure=c('p.Specificity','Beff.Specificity','Ceff.Specificity'))
}else {
	pldata <- melt(difdata, id=c('Name','Num.Items','Num.Students','DIF.Type','Sim.Type','Test','Mag','Num.DIF.Items'), measure=c('Specificity'))
}
pldata <- pldata %>% 
	rename(Measure = variable) %>%
	group_by(Mag,Num.DIF.Items,Measure) %>%
	summarize(value.mn = mean(value, na.rm=TRUE)) %>%
	print()

#Plot different number of dif items separately
m2plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 2,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(a) 2 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Specificity', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m4plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 4,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(b) 4 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Specificity', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m6plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 6,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(c) 6 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Specificity', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m8plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 8,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(d) 8 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Specificity', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m10plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 10,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(f) 10 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Specificity', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')

mplot <- plot_grid(m2plot+theme(legend.position='none'),m4plot+theme(legend.position='none'),m6plot+theme(legend.position='none'),m8plot+theme(legend.position='none'),m10plot+theme(legend.position='none'), ncol=2)
legend <- get_legend(m2plot+guides(color = guide_legend(nrow=1))+theme(legend.position='bottom'))
mplot <- plot_grid(title, mplot, legend, ncol=1, rel_heights=c(.05,.9,.05))
if (arg$test == 'binMH'){
	ggsave(file=paste0('Specificity-Plots-',arg$purify,'-',arg$padj,'-',arg$nbins,'.pdf'), path=paste0('analysisout/',arg$name,'/',arg$test,'/'), mplot, width=7.5, height=7.5, units='in')
}else {
	ggsave(file=paste0('Specificity-Plots-',arg$purify,'-',arg$padj,'.pdf'), path=paste0('analysisout/',arg$name,'/',arg$test,'/'), mplot, width=7.5, height=7.5, units='in')
}

######################################
#Plot accuracies against DIF magnitude
######################################
if (grepl('MH',arg$test) | grepl('SIB',arg$test)){
	pldata <- melt(difdata, id=c('Name','Num.Items','Num.Students','DIF.Type','Sim.Type','Test','Mag','Num.DIF.Items'), measure=c('p.Accuracy','Beff.Accuracy','Ceff.Accuracy'))
}else {
	pldata <- melt(difdata, id=c('Name','Num.Items','Num.Students','DIF.Type','Sim.Type','Test','Mag','Num.DIF.Items'), measure=c('Accuracy'))
}
pldata <- pldata %>% 
	rename(Measure = variable) %>%
	group_by(Mag,Num.DIF.Items,Measure) %>%
	summarize(value.mn = mean(value, na.rm=TRUE)) %>%
	print()

#Plot different number of dif items separately
m2plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 2,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(a) 2 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Accuracy', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m4plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 4,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(b) 4 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Accuracy', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m6plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 6,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(c) 6 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Accuracy', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m8plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 8,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(d) 8 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Accuracy', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m10plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 10,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(f) 10 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Accuracy', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')

mplot <- plot_grid(m2plot+theme(legend.position='none'),m4plot+theme(legend.position='none'),m6plot+theme(legend.position='none'),m8plot+theme(legend.position='none'),m10plot+theme(legend.position='none'), ncol=2)
legend <- get_legend(m2plot+guides(color = guide_legend(nrow=1))+theme(legend.position='bottom'))
mplot <- plot_grid(title, mplot, legend, ncol=1, rel_heights=c(.05,.9,.05))
if (arg$test == 'binMH'){
	ggsave(file=paste0('Accuracy-Plots-',arg$purify,'-',arg$padj,'-',arg$nbins,'.pdf'), path=paste0('analysisout/',arg$name,'/',arg$test,'/'), mplot, width=7.5, height=7.5, units='in')
}else {
	ggsave(file=paste0('Accuracy-Plots-',arg$purify,'-',arg$padj,'.pdf'), path=paste0('analysisout/',arg$name,'/',arg$test,'/'), mplot, width=7.5, height=7.5, units='in')
}

####################################################
#Plot positive predicted value against DIF magnitude
####################################################
if (grepl('MH',arg$test) | grepl('SIB',arg$test)){
	pldata <- melt(difdata, id=c('Name','Num.Items','Num.Students','DIF.Type','Sim.Type','Test','Mag','Num.DIF.Items'), measure=c('p.PPV','Beff.PPV','Ceff.PPV'))
}else {
	pldata <- melt(difdata, id=c('Name','Num.Items','Num.Students','DIF.Type','Sim.Type','Test','Mag','Num.DIF.Items'), measure=c('PPV'))
}
pldata <- pldata %>% 
	rename(Measure = variable) %>%
	group_by(Mag,Num.DIF.Items,Measure) %>%
	summarize(value.mn = mean(value, na.rm=TRUE)) %>%
	print()

#Plot different number of dif items separately
m2plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 2,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(a) 2 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Positive Predictive Value', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m4plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 4,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(b) 4 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Positive Predictive Value', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m6plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 6,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(c) 6 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Positive Predictive Value', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m8plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 8,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(d) 8 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Positive Predictive Value', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m10plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 10,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(f) 10 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Positive Predictive Value', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')

mplot <- plot_grid(m2plot+theme(legend.position='none'),m4plot+theme(legend.position='none'),m6plot+theme(legend.position='none'),m8plot+theme(legend.position='none'),m10plot+theme(legend.position='none'), ncol=2)
legend <- get_legend(m2plot+guides(color = guide_legend(nrow=1))+theme(legend.position='bottom'))
mplot <- plot_grid(title, mplot, legend, ncol=1, rel_heights=c(.05,.9,.05))
if (arg$test == 'binMH'){
	ggsave(file=paste0('PPV-Plots-',arg$purify,'-',arg$padj,'-',arg$nbins,'.pdf'), path=paste0('analysisout/',arg$name,'/',arg$test,'/'), mplot, width=7.5, height=7.5, units='in')
}else {
	ggsave(file=paste0('PPV-Plots-',arg$purify,'-',arg$padj,'.pdf'), path=paste0('analysisout/',arg$name,'/',arg$test,'/'), mplot, width=7.5, height=7.5, units='in')
}

####################################################
#Plot negative predicted value against DIF magnitude
####################################################
if (grepl('MH',arg$test) | grepl('SIB',arg$test)){
	pldata <- melt(difdata, id=c('Name','Num.Items','Num.Students','DIF.Type','Sim.Type','Test','Mag','Num.DIF.Items'), measure=c('p.NPV','Beff.NPV','Ceff.NPV'))
}else {
	pldata <- melt(difdata, id=c('Name','Num.Items','Num.Students','DIF.Type','Sim.Type','Test','Mag','Num.DIF.Items'), measure=c('NPV'))
}
pldata <- pldata %>% 
	rename(Measure = variable) %>%
	group_by(Mag,Num.DIF.Items,Measure) %>%
	summarize(value.mn = mean(value, na.rm=TRUE)) %>%
	print()

#Plot different number of dif items separately
m2plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 2,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(a) 2 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Negative Predictive Value', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m4plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 4,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(b) 4 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Negative Predictive Value', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m6plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 6,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(c) 6 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Negative Predictive Value', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m8plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 8,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(d) 8 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Negative Predictive Value', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')
m10plot <- ggplot()+geom_point(data=pldata[pldata$Num.DIF.Items == 10,], mapping=aes(x=Mag,y=value.mn,color=Measure,shape=Measure), size=3)+scale_shape_manual(values=ggshapes[1:length(unique(pldata$Measure))])+labs(title=paste0('(f) 10 DIF Items'))+scale_x_continuous(name='DIF Magnitude')+scale_y_continuous(name='Mean Negative Predictive Value', limits=c(0,1))+theme_bw()+theme(text=element_text(family='serif'))#, legend.position='none')

mplot <- plot_grid(m2plot+theme(legend.position='none'),m4plot+theme(legend.position='none'),m6plot+theme(legend.position='none'),m8plot+theme(legend.position='none'),m10plot+theme(legend.position='none'), ncol=2)
legend <- get_legend(m2plot+guides(color = guide_legend(nrow=1))+theme(legend.position='bottom'))
mplot <- plot_grid(title, mplot, legend, ncol=1, rel_heights=c(.05,.9,.05))
if (arg$test == 'binMH'){
	ggsave(file=paste0('NPV-Plots-',arg$purify,'-',arg$padj,'-',arg$nbins,'.pdf'), path=paste0('analysisout/',arg$name,'/',arg$test,'/'), mplot, width=7.5, height=7.5, units='in')
}else {
	ggsave(file=paste0('NPV-Plots-',arg$purify,'-',arg$padj,'.pdf'), path=paste0('analysisout/',arg$name,'/',arg$test,'/'), mplot, width=7.5, height=7.5, units='in')
}














#Curious about runtime 
end <- Sys.time()
hrdiff <- as.numeric(difftime(end, start, units = 'hours'))
mindiff <- as.numeric(difftime(end, start, units = 'mins'))
secdiff <- as.numeric(difftime(end, start, units = 'secs'))
print_color(paste0('Runtime: ',floor(hrdiff),' hours ',floor(mindiff %% 60),' mins ',round(secdiff %% 60),' seconds\n'),'bgreen')
