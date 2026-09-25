#!/usr/bin/env Rscript
#Above line allows code to be run using ./DIFSimulator.R in terminal

#Curious about runtime
start <- Sys.time()

#Libraries and what they are used for commented next to them
library(dplyr)#as_tibble and many other dataframe manipulation shortcuts
library(data.table)#setnames function
library(insight)#print_color function
library(argparser)#anything parser related
source('src/MultiPDF.R')#use of multiple random pdfs 

#Adding argument parsers so that I can vary the simulated data from the command line
parser <- arg_parser('Options for varying the simulated data generated')
parser <- add_argument(parser, "--name", help = 'name of output when in flex/run mode if other name desired',nargs='*',default='TEST')

#DIF arguments
parser <- add_argument(parser, "--DIFtype", help = 'uniform, nonuniform, or mixed: default is uniform',nargs='*',default='uniform')
parser <- add_argument(parser, "--DIFsimtype", help = 'simple or variable: default is simple',nargs='*',default='simple')
parser <- add_argument(parser, "--DIFmag", help = 'DIF magnitude: default is .5',nargs='*',default=c(.5,.5,0))
parser <- add_argument(parser, "--DIFfocstprop", help = 'proportion of students in focal group: default is .5',nargs='*',default=c(.5,.5,0))
parser <- add_argument(parser, "--DIFconsist", help = 'probability of seeing DIF in unfavored group of students: default is 1',nargs='*',default=c(1,1,0))
parser <- add_argument(parser, "--DIFitprop", help = 'proportion of items being DIF: default is .2',nargs='*',default=c(.2,.2,0))
parser <- add_argument(parser, "--DIFrefitprop", help = 'proportion of DIF items favoring reference group: default is .1',nargs='*',default=c(.1,.1,0))

#additional arguments
parser <- add_argument(parser, "--nitems", help = 'number of items when in flex mode: format input as begin,end,increment',nargs='*',default=c(20,20,0))
parser <- add_argument(parser, "--ns", help = 'number of students when in flex mode: format input as begin,end,increment',nargs='*',default=c(1000,1000,0))
parser <- add_argument(parser, "--nrun", help = 'number of runs if using variable: default is 10',nargs='*',default=c(10,10,0))
arg <- parse_args(parser)

#Turning multiple input arguments into vectors
DIFmag <- seq(from = arg$DIFmag[1], to = arg$DIFmag[2], by = arg$DIFmag[3])
DIFfocstprop <- seq(from = arg$DIFfocstprop[1], to = arg$DIFfocstprop[2], by = arg$DIFfocstprop[3])
DIFconsist <- seq(from = arg$DIFconsist[1], to = arg$DIFconsist[2], by = arg$DIFconsist[3])
DIFitprop <- seq(from = arg$DIFitprop[1], to = arg$DIFitprop[2], by = arg$DIFitprop[3])
DIFrefitprop <- seq(from = arg$DIFrefitprop[1], to = arg$DIFrefitprop[2], by = arg$DIFrefitprop[3])
nitems <- seq(from = arg$nitems[1], to = arg$nitems[2], by = arg$nitems[3])
numst <- seq(from = arg$ns[1], to = arg$ns[2], by = arg$ns[3])
if (arg$DIFsimtype == 'simple'){
	print_color(paste0('!!!!!!!!!!!!!!!!!!!!RUNNING SIMPLE DIF IRT SIMULATIONS!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bcyan')
	nrun <- c(1)
}else if (arg$DIFsimtype == 'variable'){
	print_color(paste0('!!!!!!!!!!!!!!!!!!!RUNNING VARIABLE DIF IRT SIMULATIONS!!!!!!!!!!!!!!!!!!!!!!!\n'),'bcyan')
	nrun <- seq(from = arg$nrun[1], to = arg$nrun[2], by = arg$nrun[3])
}else {
	print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!INCORRECT INPUT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bred')
	break
}
if (arg$DIFtype == 'uniform'){
	print_color(paste0('!!!!!!!!!!!!!!!!!!!!RUNNING UNIFORM DIF IRT SIMULATIONS!!!!!!!!!!!!!!!!!!!!!!!\n'),'bcyan')
}else if (arg$DIFtype == 'nonuniform'){
	print_color(paste0('!!!!!!!!!!!!!!!!!!RUNNING NONUNIFORM DIF IRT SIMULATIONS!!!!!!!!!!!!!!!!!!!!!!\n'),'bcyan')
}else if (arg$DIFtype == 'mixed'){
	print_color(paste0('!!!!!!!!!!!!!!!!NOT RUNNING MIXED DIF IRT SIMULATIONS YET!!!!!!!!!!!!!!!!!!!!!\n'),'bred')
	break
}else {
	print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!INCORRECT INPUT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bred')
	break
}
###############################################################################################################
##################################################FUNCTIONS####################################################
###############################################################################################################
#Function to generate student repsonse for item based on whether the item is DIF and the student's group
response <- function(itempar,th=0,g,consist){
	set.seed(NULL)#used to unset the seed so the responses will differ for each item
	resp <- runif(1, min=0, max=1)
	dif <- runif(1, min=0, max=1)
	
	#Set probabilities based on model
	if (g == 0){
		a <- itempar$Discrimination
		b <- itempar$Difficulty
	}else if (g == 1){
		if (dif < consist){
			a <- itempar$Focal.Discrimination
			b <- itempar$Focal.Difficulty
		}else {
			a <- itempar$Discrimination
			b <- itempar$Difficulty
		}
	}
	p <- exp(a*(th - b))/(1 + exp(a*(th - b)))
	
	#Set student response based on random uniform dice roll and probability
	if (resp < p){
		r <- 1
	}else {
		r <- 0
	}	

	return(r)
}

###############################################################################################################
##################################################SIMULATION###################################################
###############################################################################################################
#Simulating data 
#Flexible simulations based on user input

#Build datasets
for (nit in nitems){
	for (nst in numst){

		#Keep track of different files
		filecount <- 1
	
		for (mag in DIFmag){
			for (focstprop in DIFfocstprop){
				for (consist in DIFconsist){
					for (itprop in DIFitprop){
						for (refitprop in DIFrefitprop){
							for (r in nrun){
								#Generate data
								#Using while loop to ensure simulated data will converge in IRT fitting 
								DataCheck <- TRUE
								while (DataCheck){
									#Setting incremented values
									Item <- paste0('Item',1:nit)
									ns <- nst
									thmn <- c(0)
									thsd <- c(1)
									thw <- c('eq')
								
									#Setting DIF items
									difitems <- rep(0, nit)
									if (mag > 0){
										ndifitems <- round(itprop*nit)
										nrefdifitems <- round(refitprop*ndifitems)
										switch <- round(seq(from = 1, to = nit, length.out = ndifitems))
										refswitch <- sample(switch, nrefdifitems)
										for (i in switch){difitems[i] <- 1}
										for (i in refswitch){difitems[i] <- 2}
									}
									print(difitems)

									#Setting base item parameters
									if (arg$DIFsimtype == 'simple'){
										itemdiff <- rep(0, nit)
										itemdisc <- rep(2, nit)
									}else {
										itemdiff <- runif(nit, min=-1.5, max=1.5)
										itemdisc <- runif(nit, min=1.5, max=3.5)
									}

									#Saving item generators
									if (!dir.exists(paste0('simdata/',arg$name,'/',nit,'items/',ns,'students/',arg$DIFtype,'/',arg$DIFsimtype))){dir.create(paste0('simdata/',arg$name,'/',nit,'items/',ns,'students/',arg$DIFtype,'/',arg$DIFsimtype), recursive = TRUE)}
									if (arg$DIFsimtype == 'simple'){
										generators <- data.frame(Num.Items = nit, Num.Students = nst, DIF.Type = arg$DIFtype, DIF.Magnitude = mag, Prop.Focal.St = focstprop, DIF.Consistency = consist, Prop.DIF.Items = itprop, Prop.DIF.Ref.Items = refitprop)
									}else {
										generators <- data.frame(Num.Items = nit, Num.Students = nst, Num.Run = r, DIF.Type = arg$DIFtype, DIF.Magnitude = mag, Prop.Focal.St = focstprop, DIF.Consistency = consist, Prop.DIF.Items = itprop, Prop.DIF.Ref.Items = refitprop)
									}
									write.csv(generators, paste0('simdata/',arg$name,'/',nit,'items/',ns,'students/',arg$DIFtype,'/',arg$DIFsimtype,'/',paste0(arg$name,filecount),'-Generators.csv'), row.names = FALSE)	

									#True item parameters that will be used in the generated data
									par <- data.frame(Items = Item, DIF.Item = difitems, Difficulty = itemdiff, Discrimination = itemdisc)

									#Calculating focal item parameters
									if (arg$DIFtype == 'uniform'){
										par$Focal.Difficulty <- ifelse(par$DIF.Item == 1,(par$Difficulty + mag),par$Difficulty)
										par$Difficulty <- ifelse(par$DIF.Item == 2,(par$Difficulty + mag),par$Difficulty)
										par$Focal.Discrimination <- par$Discrimination
									}else if (arg$DIFtype == 'nonuniform'){
										par$Focal.Difficulty <- par$Difficulty
										par$Focal.Discrimination <- ifelse(par$DIF.Item == 1,(par$Discrimination - mag),par$Discrimination)
										par$Discrimination <- ifelse(par$DIF.Item == 2,(par$Discrimination - mag),par$Discrimination)
									}
									
									print_color(paste0('==============================================================================\n'),'bold')
									print_color(paste0('==============================Item Parameters=================================\n'),'bold')
									print_color(paste0('==============================================================================\n'),'bold')
									print(par)
									write.csv(par, paste0('simdata/',arg$name,'/',nit,'items/',ns,'students/',arg$DIFtype,'/',arg$DIFsimtype,'/',paste0(arg$name,filecount),'-Items.csv'), row.names = FALSE)	

									#Setting true proficiencies
									df <- data.frame(ID = 1:ns, Theta = multirnorm(ns, mean=0, sd=1), Group = ifelse(runif(ns, min=0, max=1) < focstprop,1,0))

									#Fill in student responses 
									print_color(paste0('==============================================================================\n'),'bcyan')
									print_color(paste0('========================Generating Student Responses==========================\n'),'bcyan')
									print_color(paste0('==============================================================================\n'),'bcyan')
									for (j in Item){
										temp <- c()
										
										for (i in 1:ns){
											resp <- response(itempar=par[par$Items == j,], th=df[df$ID == i,]$Theta, g=df[df$ID == i,]$Group, consist=consist)
											temp <- c(temp, resp)
										}
										df[[j]] <- temp
									}
									print(as_tibble(df))

									checkvec <- c()
									#Check data quality
									for (j in Item){
										freq <- table(df[[j]])
										#print(freq)#will suppress after testing
										#print(freq[1] > 50 & freq[2] > 50)#will suppress after testing
										checkvar <- freq[1] > 50 & freq[2] > 50
										checkvec <- c(checkvec,checkvar)
									}
									print(checkvec)#will suppress after testing
									print(all(checkvec))#will suppress after testing
									check <- all(checkvec)
									if (is.na(check)){
										check <- FALSE
									}	
									if (check){
										print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!!DATA PASSED QUALITY CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bgreen')
										DataCheck <- FALSE
									}else {
										print_color(paste0('!!!!!!!!!!!!!!!!!!!!!!!!DATA FAILED QUALITY CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'),'bred')
									}
									
								}#end of while loop to check data quality 
			
								#Saving flex datasets
								write.csv(df, paste0('simdata/',arg$name,'/',nit,'items/',ns,'students/',arg$DIFtype,'/',arg$DIFsimtype,'/',paste0(arg$name,filecount),'-Data.csv'), row.names = FALSE)
								filecount <- filecount + 1
							}#end of nrun loop
						}#end of refitprop loop
					}#end of itprop loop
				}#end of consist loop
			}#end of focstprop loop
		}#end of mag loop
	}#end of ns loop
}#end of nitems loop



#Curious about runtime 
end <- Sys.time()
hrdiff <- as.numeric(difftime(end, start, units = 'hours'))
mindiff <- as.numeric(difftime(end, start, units = 'mins'))
secdiff <- as.numeric(difftime(end, start, units = 'secs'))
print_color(paste0('Runtime: ',floor(hrdiff),' hours ',floor(mindiff %% 60),' mins ',round(secdiff %% 60),' seconds\n'),'bgreen')
