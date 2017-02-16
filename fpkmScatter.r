####################################################################################################
# 
# Mark Cembrowski, Janelia Research Campus, April 15 2015
#
# This script plots a scatterplot, coming the gene expression across two different populations.
# Additional options allow labeling of specific genes.
#
# INPUT:
# sampleX: character. String corresponding to the sample_name, plotted on x axis.
# sampleY: character. String corresponding to the sample_name, plotted on y axis.
#
# OPTIONAL INPUT: 
# gsnList=vector(): character vector.  A list of specific genes to highlight.
# fdr=-1: numeric.  FDR value to use if using DE as a means of labeling data points.
# foldThres=-1: numeric.  Minimum fold change between populations in order to label points.
# fpkmThres=-1: numeric.  Minimum FPKM value for a gene to be considered enriched.
# avgPass=T: logical.  If T, allows genes to be enriched whenever the average value between the
#	two samples exceeds provided threshold.  If F, requires enrichment to be present across
#	all replicates, instead of the average, a stricter comparison.
# maxX,maxY,minX,minY=-1: numerics.  Override axis limits on x and y axes.
# label=T: logical.  Switch to control whether highlighted genes are labeled.
#
# OUTPUT:
# A scatterplot is rendered and the corresponding matrix is returned.
#
#####################################################################################################

fpkmScatter <- function(sampleX,sampleY,gsnList=vector(),fdr=-1,foldThres=-1,fpkmThres=-1,avgPass=T,
				maxX=-1,maxY=-1,minX=-1,minY=-1,label=T){

	# Set parameters corresponding to the type of analysis.
	apriori <- F
	explore <- F
	if(length(gsnList)>0.1){ apriori <- T }
	if(fdr>0||foldThres>0){ explore <- T }
	if(apriori&&explore){
		# Check to make sure that the user also has not supplied a list of genes to plot.
		if(length(gsnList)>0.1){
			stop('Cannot simultaneously do exploratory and a priori labeling')
		}
	}

	# Produce an empty output, if not labeling any genes.
	if( (!apriori) && (!explore) ){
		enrX <- c()
		enrY <- c()
	}

	# Obtain enrichment for provided genes, if doing a priori analysis.
	if(apriori){
		expX <- fpkmPoolMat[unlist(symToId(gsnList)),sampleX]
		expY <- fpkmPoolMat[unlist(symToId(gsnList)),sampleY]
		priorRat <- expX/expY
		enrX <- unlist(symToId(gsnList[priorRat>1]))
		enrY <- unlist(symToId(gsnList[priorRat<1]))
	}
	
	# Recover list of genes to plot, if doing an exploratory analysis.
	if(explore){
		enrX <- getEnrGenes(sampleX,sampleY,fdr=fdr,foldThres=foldThres,fpkmThres=fpkmThres,
			avgPass=avgPass)
		enrY <- getEnrGenes(sampleY,sampleX,fdr=fdr,foldThres=foldThres,fpkmThres=fpkmThres,
			avgPass=avgPass)
	}
	
	# Generate an informative table for output.
	enrichedRegion <- c(rep(sampleX,length(enrX)),rep(sampleY,length(enrY)))

	geneId <- c(enrX,enrY)
	geneShortName <- idToSym(geneId)
	foldDiff <- fpkmPoolMat[geneId,sampleX]/fpkmPoolMat[geneId,sampleY]
	foldDiff <- pmax(foldDiff,1/foldDiff) # Take convention that is >1.
	sampleXFpkm <- fpkmPoolMat[geneId,sampleX]
	sampleYFpkm <- fpkmPoolMat[geneId,sampleY]

	enrOut <- as.data.frame(cbind(enrichedRegion,geneId,geneShortName,foldDiff,
		sampleXFpkm,sampleYFpkm))
	enrOut$foldDiff <- as.numeric(as.character(enrOut$foldDiff))
	enrOut <- enrOut[with(enrOut,order(enrichedRegion,-foldDiff)),]	
	colnames(enrOut)[5] <- paste(sampleX,'Fpkm',sep='')
	colnames(enrOut)[6] <- paste(sampleY,'Fpkm',sep='')

	if(apriori){
		toHighlight <- unlist(symToId(gsnList))
	}else{
		if(explore){
			toHighlight <- enrOut$geneId
		}else{
			toHighlight <- vector()
		}
	}
	highlight <- rownames(fpkmPoolMat)%in%toHighlight

	# Prep matrix for plotting.
	df <- fpkmPoolMat[,c(sampleX,sampleY)]
	df <- cbind(df,highlight)	


	# Do additive smoothing for logarithmic plot.
	df[,sampleX] <- df[,sampleX] + 1
	df[,sampleY] <- df[,sampleY] + 1

	# Sort so that highlighted genes show up on top layers of plot.
	df <- df[order(df$highlight),]

	# Attach gene short names to data frame.
	df[,'gene_short_name'] <- idToSym(rownames(df))
	
	# Adjust column names to properly interact with aes_string.
	colnames(df)[c(1,2)] <- c('fpkmX','fpkmY')

	# Determine correlation coefficient.
	xCorr <- round(cor(x=df[,'fpkmX'],y=df[,'fpkmY']),digits=2)

	# Plot.  NOTE: need to use aes_string, as we're plotting inside a function!
	scatPlot <- ggplot(df,aes_string(x='fpkmX',y='fpkmY'))
	scatPlot <- scatPlot + scale_colour_manual(values=c("black","red"))
	scatPlot <- scatPlot + geom_abline(slope=1,linetype='dashed',weight=0.1)

	# Plot most points as transparent, highlighted genes as nontransparent
	scatPlot <- scatPlot + geom_point(data=subset(df,highlight==0),aes(colour=factor(highlight)),size=1,alpha=0.1)
	scatPlot <- scatPlot + geom_point(data=subset(df,highlight==1),aes(colour=factor(highlight)),size=1)
	

	theMax <- max(df[,c('fpkmX','fpkmY')])
	theMin <- 0.9
	theMaxX <- theMax
	theMaxY <- theMax
	theMinX <- theMin
	theMinY <- theMin
	if (maxX>0){theMaxX <- maxX}
	if (maxY>0){theMaxY <- maxY}
	if (minX>0){if (minX<maxX){theMinX <- minX}}
	if (minY>0){if (minY<maxY){theMinY <- minY}}

	scatPlot <- scatPlot + scale_x_log10() + scale_y_log10()
	scatPlot <- scatPlot + coord_cartesian(xlim=c(theMinX,theMaxX),ylim=c(theMinY,theMaxY))
	scatPlot <- scatPlot+geom_hline(aes(yintercept=11))
	scatPlot <- scatPlot+geom_vline(aes(xintercept=11))	

	
	if(label){
		scatPlot <- scatPlot + geom_text(data=subset(df,highlight==1),
			aes_string(label="gene_short_name"),hjust=0,vjust=0,size=2,colour="black")
	}
	scatPlot <- scatPlot + labs(x=paste('FPKM+1, ',sampleX,sep=""),y=paste('FPKM+1, ',sampleY,sep=""),
		title=paste('FPKM scatter plot, selected genes in blue; CC=',xCorr,sep=""))
	scatPlot <- scatPlot + theme_bw()
	scatPlot <- scatPlot + theme(legend.position='none')

	print(scatPlot)
	
	
	invisible(enrOut)
}
