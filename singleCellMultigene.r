###############################################################
#
# Mark Cembrowski, Janelia Research Campus, Nov 7 2016
#
# This script plots the average expression, +/- SD, across
# a list of genes for each cell type in the dataset.
#
###############################################################

singleCellMultigene <- function(theGenes,rMax=-1,quiet=T){

	theClasses <- c('pyramidal CA1','oligodendrocytes',
			'interneurons','astrocytes_ependymal',
			'microglia','endothelial-mural')

	subMeta <- subset(scMeta,tissue=='ca1hippocampus')

	allStats <- as.data.frame(matrix(nrow=0,ncol=4))
	allMeans <- as.data.frame(matrix(nrow=0,ncol=nrow(subMeta)))

	geneInMat <- theGenes[theGenes%in%rownames(scHip)]
	subMat <- scHip[geneInMat,]
	
	for (curClass in theClasses){
		toKeep <- subMeta$level1class==curClass

		# Print out total number of entries.
		if(!quiet){
			print(paste('For class',curClass,'the number of cells is',sum(toKeep)))	
		}

		curMat <- subMat[,toKeep]
		curStats <- c( mean(apply(curMat,1,mean)),
				sd(apply(curMat,1,mean)),
				mean(apply(curMat,1,mean))-sd(apply(curMat,1,mean)),
				mean(apply(curMat,1,mean))+sd(apply(curMat,1,mean)))
	
		allStats <- rbind(allStats,curStats)

		curMeans <- apply(curMat,1,mean)
		allMeans <- rbind(allMeans,curMeans)
		colnames(allMeans) <- names(curMeans)
		rownames(allMeans)[nrow(allMeans)] <- curClass
	}
	rownames(allStats) <- theClasses
	colnames(allStats) <- c('mean','sd','lo','hi')
	allStats[,'grp'] <- rownames(allStats)

	allMeansMelt <- melt(as.matrix(allMeans))
	colnames(allMeansMelt) <- c('grp','gene','mean')


	colourRule <- c('astrocytes_ependymal'='blue',
			'interneurons'='purple','microglia'='green',
			'oligodendrocytes'='brown','endothelial-mural'='red','pyramidal CA1'='grey')

	gg <- ggplot(allStats,aes(x=grp,y=mean)) + geom_bar(stat='identity',aes(colour=grp),fill='white')
	gg <- gg + geom_point(data=allMeansMelt,alpha=0.3,aes(colour=grp))
	gg <- gg + geom_errorbar(aes(ymin=mean,ymax=hi),width=0.5)
	gg <- gg + theme_bw() + scale_color_manual(values=colourRule)
	gg <- gg + scale_x_discrete(limits=c('pyramidal CA1','interneurons','oligodendrocytes',
				'astrocytes_ependymal','endothelial-mural','microglia'))
	if(rMax>0){
		gg <- gg + scale_y_continuous(lim=c(0,rMax))
	}
	print(gg)
	
}
