###########################################################################
#
# Mark Cembrowski, Janelia Research Campus, Nov 5 2016
#
# This script plots the single cell molecule count values for a given
# gene or set of genes.
#
###########################################################################

singleCellBarplot <- function(theGenes,rMax=-1,quiet=F){
	subMat <- scCa1[theGenes,]

	if(!quiet){
		# Print maximum values for genes.
		print('Maximum vals across all cells:')
		maxes <- apply(subMat,1,max)
		names(maxes) <- theGenes
		print(maxes)
	}

	stats <- cbind(
		apply(subMat,1,mean),
		apply(subMat,1,sd),
		apply(subMat,1,mean)-apply(subMat,1,sd),
		apply(subMat,1,mean)+apply(subMat,1,sd) )
	colnames(stats) <- c('mean','sd','lo','hi')
	stats <- as.data.frame(stats)
	stats[,'gene'] <- rownames(stats)

	# Set x limits for plotting; this will prevent any genes that are
	# not present in the Ziesel data from being labeled NA.
	theLimits <- theGenes
	
	gg <- ggplot(stats,aes(x=gene,y=mean)) 
	gg <- gg + geom_bar(stat='identity',colour='black',fill='white')
	gg <- gg + geom_errorbar(aes(ymin=lo,ymax=hi),width=0.5)
	gg <- gg + theme_bw() 
	if(rMax>-1){
		gg <- gg + scale_y_continuous(lim=c(min(stats$lo),rMax))
	}
	gg <- gg + scale_x_discrete(limits=theLimits)
	gg <- gg + xlab('Gene') + ylab('Mean expression across CA1 PCs')
	
	print(gg)
}
