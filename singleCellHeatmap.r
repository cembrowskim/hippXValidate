####################################################################################
#
# Mark Cembrowski, Janelia Research Campus, Nov 11 2016
#
# This function renders a heatmap of genes for single cell data. 
#
####################################################################################

singleCellHeatmap <- function(gsnList,label=T,rMax=-1){

	theClasses <- c('pyramidal CA1','oligodendrocytes',
			'interneurons','astrocytes_ependymal',
			'microglia','endothelial-mural')

	subMeta <- subset(scMeta,tissue=='ca1hippocampus')

	geneInMat <- gsnList[gsnList%in%rownames(scHip)]
	subMat <- scHip[geneInMat,]
	
	allStats <- as.data.frame(matrix(nrow=0,ncol=length(geneInMat)))

	for (curClass in theClasses){
		toKeep <- subMeta$level1class==curClass

		curMat <- subMat[,toKeep]
		curStats <- apply(curMat,1,mean)

		allStats <- rbind(allStats,curStats)
		colnames(allStats) <- names(curStats)
		rownames(allStats)[nrow(allStats)] <- curClass
	}

	# Sort data for plotting.
	theOrder <- apply(allStats,2,max)
	allStats <- allStats[,order(theOrder)]

	# Order
	thres1 <- sum(apply(allStats,2,max)<1)
	thres01 <- sum(apply(allStats,2,max)<0.1)
	print(paste('From',ncol(allStats),'genes:'))
	print(paste(thres1,'less than 1 mol/avg'))
	print(paste(thres01,'less than 0.1mol/avg'))

	# Melt data frame for plotting.
	allStatsMelt <- as.data.frame(melt(as.matrix(allStats)))
	colnames(allStatsMelt) <- c('grp','gene','count')
	
	if(rMax>0){
		# Reset for plotting purposes in ggplot.
		allStatsMelt$count[allStatsMelt$count>rMax] <- rMax
	}

	# Plot.
	heatGr <- ggplot(allStatsMelt,aes(x=grp,y=gene)) 
	heatGr <- heatGr + geom_tile(aes(fill=count),colour='white')
	if(rMax>0){ 
		theRangeMax <- rMax + 0.0001 
	}else{
		theRangeMax <- max(allStatsMelt$count)+0.0001
	}


	# Define other stylistic features.
	heatGr <- heatGr + scale_fill_gradientn(colours=c('cyan','white','red'),
		values=c(0,0.5,1),guide='colorbar',limits=c(0,theRangeMax)) # Changed from green and purple.
	heatGr <- heatGr + theme(panel.background=element_blank())
	heatGr <- heatGr + theme(axis.text.x = element_text(angle=45,hjust=1))
	heatGr <- heatGr + labs(x='Cell type',y='Gene')
	heatGr <- heatGr + scale_x_discrete(limits=c('pyramidal CA1','interneurons',
			'oligodendrocytes','astrocytes_ependymal','endothelial-mural',
			'microglia'))
	print(heatGr)

	invisible(heatGr)
}
