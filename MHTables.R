#!/usr/bin/env Rscript
#Above line allows code to be run using ./MHTables.R in terminal

#Curious about runtime
start <- Sys.time()

#Libraries and what they are used for commented next to them
library(dplyr)#as_tibble and many other dataframe manipulation shortcuts
library(tidyr)#unnest function
library(data.table)#setnames function
library(reshape2)#melt function
library(insight)#print_color function
library(ggplot2)#plot related
library(geomtextpath)#geom_text_segment
library(ggrepel)#geom_text_repel
library(treemapify)#treemap plots
library(ggmosaic)#mosaic plots
library(cowplot)#combining plots
library(argparser)#argument parser stuff

#Options for varying analysis in terminal
parser <- arg_parser('User Inputs for DIF Analyses')
parser <- add_argument(parser, '--post', help = 'analyzing posttest FCI if TRUE: default is TRUE',nargs='*',default=TRUE)
parser <- add_argument(parser, '--nbins', help = 'number of bins for the Binned MH Test: default is 5',nargs='*',default=5)
arg <- parse_args(parser)

#Redefining arguments that may have multiple inputs
nbins <- as.numeric(arg$nbins)

#Collecting data
print_color('============================================================================\n','bold')
print_color('==============================Collecting Data===============================\n','bold')
print_color('============================================================================\n','bold')

if (arg$post){
	df <- read.csv('realdata/FCI-post.csv')
}else {
	df <- read.csv('realdata/FCI-pre.csv')
}
Items <- paste0('Item',1:30)
df$Total <- apply(df[,Items],1,sum)
print(head(df))

#Making MH Contingency Tables
print_color('============================================================================\n','bcyan')
print_color('=============================Contingency Tables=============================\n','bcyan')
print_color('============================================================================\n','bcyan')

data <- df
data$Gender <- ifelse(data$Gender == 'M',0,1)
print(head(data))

#Automating bin selection based on user input
bins <- unname(round(quantile(data$Total, probs = seq(0,1, by = 1/nbins)[2:nbins])))
print_color(paste0('Number of bins being used: ',nbins,'\n'),'bgreen')
print_color(paste0('Scores being used to separate bins: ',paste(bins, collapse = ', '),'\n'),'bgreen')
for (b in 1:nbins){
	if (b == 1){
		data$Total <- ifelse((data$Total < bins[b]),b,data$Total)
	}else if (b == nbins){
		data$Total <- ifelse((data$Total >= bins[(b-1)]),b,data$Total)
	}else {
		data$Total <- ifelse((data$Total >= bins[(b-1)] & data$Total < bins[b]),b,data$Total)
	}
}	

#Storing values for barcharts
n11vec <- c()
n00vec <- c()
n10vec <- c()
n01vec <- c()
itemvec <- c()
binvec <- c()

difvalues <- c()
for (i in Items){
	data[[i]] <- factor(data[[i]], levels = c(1,0))
	data$Gender <- factor(data$Gender, levels = c(0,1))
	alphanum <- c()
	alphaden <- c()
	testnum <- c()
	testden <- c()
	#Calculations for each bin
	for (b in 1:nbins){
		tab <- table(data[data$Total == b, c('Gender',i)])
		#print(tab)
		n11vec <- c(n11vec, tab[2,1])
		n00vec <- c(n00vec, tab[1,2])
		n10vec <- c(n10vec, tab[1,1])
		n01vec <- c(n01vec, tab[2,2])
		itemvec <- c(itemvec, i)
		binvec <- c(binvec, b)
		ni <- tab[1,1]+tab[1,2]+tab[2,1]+tab[2,2]
		alphanum <- append(alphanum, (tab[1,1]*tab[2,2]/ni))
		alphaden <- append(alphaden, (tab[1,2]*tab[2,1]/ni))
		testnum <- append(testnum, (tab[1,1] - (tab[1,1]+tab[1,2])*(tab[1,1]+tab[2,1])/ni))
		testden <- append(testden, ((tab[1,1]+tab[1,2])*(tab[2,1]+tab[2,2])*(tab[1,1]+tab[2,1])*((tab[1,2]+tab[2,2])/((ni**2)*(ni-1)))))
	}
	#Calculations for item
	alphaMH <- sum(alphanum)/sum(alphaden)
	deltaMH <- -2.35*log(alphaMH)
	if (abs(deltaMH) < 1){
		ETSDeltaCode <- 'A'
	}else if (abs(deltaMH) < 1.5){
		ETSDeltaCode <- 'B'
	}else if (abs(deltaMH) >= 1.5){
		ETSDeltaCode <- 'C'
	}

	MHtest <- ((abs(sum(testnum)) - .5)**2)/sum(testden)
	pvalue <- 1 - pchisq(MHtest, df=1)
	itemvalues <- c(round(alphaMH,4),round(deltaMH,4),ETSDeltaCode,round(MHtest,4),round(pvalue,4))
	difvalues <- append(difvalues, itemvalues)
}
rownam <- Items
colnam <- c('alphaMH','deltaMH','ETS.Code','MH.Stat', 'P.value')
binmhmat <- matrix(difvalues, nrow = length(Items), byrow = TRUE, dimnames = list(rownam,colnam))
binmhdf <- as.data.frame(binmhmat)
colnames(binmhdf) <- colnam
binmhdf$Sig <- rep(' ', times = length(Items))
index <- order(as.numeric(binmhdf$P.value), decreasing = FALSE)
rank <- 1

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
sigcodes <- 'Significance Codes For Global Alphas: 0 \'***\' .001 \'**\' .01 \'*\' .05 \'.\' .1 \' \' 1\n'
	
#Finding dif items through significance tests
binmhdif <- binmhdf[as.numeric(binmhdf$P.value) < .05,]
binmhdifitems <- row.names(binmhdif)
binmhdifitemsout <- ifelse(length(row.names(binmhdif)) == 0,NA,paste(gsub('s','',row.names(binmhdif)),collapse=', '))

#Finding large effect size dif items
binmheffdif <- binmhdf[binmhdf$ETS.Code == 'C',]
binmheffdifitems <- row.names(binmheffdif)
binmheffdifitemsout <- ifelse(length(row.names(binmheffdif)) == 0,NA,paste(gsub('s','',row.names(binmheffdif)),collapse=', '))

#Outputting Binned MH DIF Table
print(binmhdf)
cat(sigcodes)
print_color(paste0('Items detected as DIF items: ',binmhdifitemsout,'\n'),'bold')

#Dissecting table information for help with plots
favorrefitems <- row.names(binmhdf[binmhdf$alpha > 1,])
largedifitems <- row.names(binmhdf[binmhdf$ETS.Code == 'C',])
print_color(paste0('Items with large DIF effect size: ',paste(largedifitems, collapse=', '),'\n'),'bred')
moddifitems <- row.names(binmhdf[binmhdf$ETS.Code == 'B',])
print_color(paste0('Items with moderate DIF effect size: ',paste(moddifitems, collapse=', '),'\n'),'byellow')


contab <- data.frame(Items = itemvec, Bins = binvec, Foc1 = n11vec, Ref1 = n10vec, Foc0 = n01vec, Ref0 = n00vec)
#print(head(contab,10))

pdf('MHTablePlots.pdf')
for (i in Items){
	pldf <- contab %>%
		filter(Items == i) %>%
		as_tibble()
	pldf <- melt(pldf, id = c('Items','Bins'))
	pldf <- pldf %>%
		rename(Cell = variable, Count = value) 

	pldf$Group <- ifelse(grepl('Ref',pldf$Cell),'Reference','Focal')
	pldf$Correct <- ifelse(grepl('1',pldf$Cell),'Correct','Incorrect')
		
	#Make treemap plots for every bin 	
	temp <- pldf %>%
		filter(Bins == 1) %>%
		group_by(Group) %>%
		mutate(Perc = round(100*Count/sum(Count),2)) 
	pl1 <- ggplot(temp, aes(area=Count, label=paste(Correct,'\n',Perc), subgroup=Group, fill=Group))+geom_treemap()+geom_treemap_subgroup_border(color='black')+geom_treemap_text(place='center')+labs(title=paste0('Bin 1'))+theme_bw()+theme(text=element_text(family='serif'))
	temp <- pldf %>%
		filter(Bins == 2) %>%
		group_by(Group) %>%
		mutate(Perc = round(100*Count/sum(Count),2)) 
	pl2 <- ggplot(temp, aes(area=Count, label=paste(Correct,'\n',Perc), subgroup=Group, fill=Group))+geom_treemap()+geom_treemap_subgroup_border(color='black')+geom_treemap_text(place='center')+labs(title=paste0('Bin 2'))+theme_bw()+theme(text=element_text(family='serif'))
	temp <- pldf %>%
		filter(Bins == 3) %>%
		group_by(Group) %>%
		mutate(Perc = round(100*Count/sum(Count),2)) 
	pl3 <- ggplot(temp, aes(area=Count, label=paste(Correct,'\n',Perc), subgroup=Group, fill=Group))+geom_treemap()+geom_treemap_subgroup_border(color='black')+geom_treemap_text(place='center')+labs(title=paste0('Bin 3'))+theme_bw()+theme(text=element_text(family='serif'))
	temp <- pldf %>%
		filter(Bins == 4) %>%
		group_by(Group) %>%
		mutate(Perc = round(100*Count/sum(Count),2)) 
	pl4 <- ggplot(temp, aes(area=Count, label=paste(Correct,'\n',Perc), subgroup=Group, fill=Group))+geom_treemap()+geom_treemap_subgroup_border(color='black')+geom_treemap_text(place='center')+labs(title=paste0('Bin 4'))+theme_bw()+theme(text=element_text(family='serif'))
	temp <- pldf %>%
		filter(Bins == 5) %>%
		group_by(Group) %>%
		mutate(Perc = round(100*Count/sum(Count),2)) 
	pl5 <- ggplot(temp, aes(area=Count, label=paste(Correct,'\n',Perc), subgroup=Group, fill=Group))+geom_treemap()+geom_treemap_subgroup_border(color='black')+geom_treemap_text(place='center')+labs(title=paste0('Bin 5'))+theme_bw()+theme(text=element_text(family='serif'))

	#Make barplot to compare counts for each bin in the contingency tables
	binplot <- ggplot(data=pldf, aes(x=Bins, y=Count, fill=Cell))+geom_bar(stat="identity")

	#Combine plots
	scoreplot <- plot_grid(binplot,pl1,pl2,pl3,pl4,pl5, ncol=2)
	if (i %in% largedifitems){
		if (i %in% favorrefitems){
			title <- ggdraw()+draw_label(paste0(i, ': Large DIF Effect Size Favoring Reference'),fontface='bold')+theme_bw()+theme(text=element_text(family='serif'), panel.border=element_blank(), plot.background=element_blank(), panel.background=element_blank())
		}else {
			title <- ggdraw()+draw_label(paste0(i, ': Large DIF Effect Size Favoring Focal'),fontface='bold')+theme_bw()+theme(text=element_text(family='serif'), panel.border=element_blank(), plot.background=element_blank(), panel.background=element_blank())
		}
	}else if (i %in% moddifitems){
		if (i %in% favorrefitems){
			title <- ggdraw()+draw_label(paste0(i, ': Moderate DIF Effect Size Favoring Reference'),fontface='bold')+theme_bw()+theme(text=element_text(family='serif'), panel.border=element_blank(), plot.background=element_blank(), panel.background=element_blank())
		}else {
			title <- ggdraw()+draw_label(paste0(i, ': Moderate DIF Effect Size Favoring Focal'),fontface='bold')+theme_bw()+theme(text=element_text(family='serif'), panel.border=element_blank(), plot.background=element_blank(), panel.background=element_blank())
		}
	}else {
		if (i %in% favorrefitems){
			title <- ggdraw()+draw_label(paste0(i, ': Negligible DIF Effect Size Favoring Reference'),fontface='bold')+theme_bw()+theme(text=element_text(family='serif'), panel.border=element_blank(), plot.background=element_blank(), panel.background=element_blank())
		}else {
			title <- ggdraw()+draw_label(paste0(i, ': Negligible DIF Effect Size Favoring Focal'),fontface='bold')+theme_bw()+theme(text=element_text(family='serif'), panel.border=element_blank(), plot.background=element_blank(), panel.background=element_blank())
		}
	}
	scoreplot <- plot_grid(title, scoreplot, ncol=1, rel_heights=c(.05,1))
	print(scoreplot)
}
dev.off()




#Curious about runtime
end <- Sys.time()
hrdiff <- as.numeric(difftime(end, start, units = 'hours'))
mindiff <- as.numeric(difftime(end, start, units = 'mins'))
secdiff <- as.numeric(difftime(end, start, units = 'secs'))
print_color(paste0('Runtime: ',floor(hrdiff),' hours ',floor(mindiff %% 60),' mins ',round(secdiff %% 60),' seconds\n'),'bgreen')


