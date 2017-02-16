################################################################################################
#
# Mark Cembrowski, Janelia Research Campus, April 16 2015
#
# This function renders a heatmap of genes.  Note that one of gsnList OR (sampleX,sampleY) must
# be provided, routing the flow to either visualise a predetermined list of genes or genes that
# are enriched/depleted across samples.
#
# 111115: added an option that contains whether both (enriched+depleted) genes are returned,
# or just enriched genes.  This is toggled in the enrOnly=F option.
#
# OPTIONAL INPUT:
# gsnList=vector(): character.  One or more gene names to plot; can be ENS IDs or gene short names.
# sampleX=vector(): character.  One or more samples, to explore enrichment/depletion as a cohort
#	relative to sampleY.
# sampleY=vector(): character.  Analogous to sampleX.
# foldThres=-1: numeric.  Minimum value for fold differences, if screening by fold change.
# fpkmThres=-1: numeric.  Minimum FPKM value for a gene in the expressed population in order
#	to be accepted for enrichment.
# fdr=-1: numeric.  FDR threshold if screening by DE.
# avgPass=T: logical.  Controls whether analysis is done at a pooled level (T) or a replicate 
#	level (F).
# doNorm=F: logical.  Normalize expression in heatmap on a gene-by-gene basis.
# mask=c(): character.  List of samples to screen from output visualisation.
# rMax=-1: numeric.  Maximum range value.
# enrOnly=F: logical. If T, return genes enriched only in sample X; if F, return genes 
#	enriched in either sample X or sample Y.
#
# OUTPUT:
# A heatmap is rendered for the input genes and returned.
#
################################################################################################

fpkmHeatmap <- function(gsnList=vector(),sampleX=vector(),sampleY=vector(),foldThres=-1,
		fpkmThres=-1,fdr=-1,avgPass=T,doNorm=F,replicates=F,mask=c(),rMax=-1,unmask=c(),
		enrOnly=F,label=T){

	# Set parameters corresponding to the type of analysis.
	apriori <- F
	explore <- F
	if(length(gsnList)>0.1){ apriori <- T }
	if(fdr>0||foldThres>0){ explore <- T }
	if(apriori&&explore){
		stop('Cannot simultaneously do exploratory and a priori labeling')
	}
	if(explore){
		if( (length(sampleX)<0.1) || (length(sampleY)<0.1) ){
			stop('Need at least two samples for comparative, exploratory analysis.')
		}
	}

	# If presupplying a gene list:
	# Convert gene list to unique Ensembl IDs, if supplied with gene short names.
	if(apriori){
		if(!isGeneId(gsnList)){
			gsnList <- symToId(gsnList)[[1]]
		}

		# Obtain FPKM values for gsnList and associated it with gene_ids (outId=T).
		if(replicates){
			merged <- fpkmRepMat[gsnList,]
		}else{
			merged <- fpkmPoolMat[gsnList,]
		}
	}
	# If doing exploration:
	# Identify genes to plot.
	if(explore){
		enrX <- getEnrGenes(sampleX,sampleY,fdr=fdr,foldThres=foldThres,fpkmThres=fpkmThres,
			avgPass=avgPass)
		enrY <- getEnrGenes(sampleY,sampleX,fdr=fdr,foldThres=foldThres,fpkmThres=fpkmThres,
			avgPass=avgPass)
		if(enrOnly){
			if(replicates){
				merged <- fpkmRepMat[enrX,]
			}else{
				merged <- fpkmPoolMat[enrX,]
			}
		}else{
			if(enrOnly){
				merged <- fpkmRepMat[c(enrX,enrY),]
			}else{
				merged <- fpkmPoolMat[c(enrX,enrY),]
			}
		}
	}

	# If imposing a mask of samples, use here.
	if(length(mask)>0.1){merged <- .maskSamples(mask,replicates=replicates,fpkmMat=merged)}

	# If unmasking specific samples, use here.
	if(length(unmask)>0.1){
		allSamps <- colnames(fpkmPoolMat)
		toMask <- allSamps%in%unmask
		toMask <- allSamps[!toMask]
		merged <- .maskSamples(toMask,replicates=replicates,fpkmMat=merged)
	}

	# Populate matrix with gene_id, gene_short_name, and FPKM values.
	merged <- cbind(rownames(merged),idToSym(rownames(merged)),merged)
	colnames(merged)[c(1,2)] <- c('gene_id','gene_short_name')

	# Generate a normalised version of this matrix and extracted melted FPKM values for later.
	normed <- sweep(merged[,-c(1,2)],1,apply(merged[,-c(1,2)],1,max),'/') # Remove 1 and 2 columns (gene info); normalise
	normed <- cbind(merged[,c(1,2)],normed)
	fpkmNorm <- melt(normed,value.var='fpkm',id.vars=c('gene_id','gene_short_name'))$value

	# Melt merged and attached normalised value.
	merged <- melt(merged,value.var='fpkm',id.vars=c('gene_id','gene_short_name'))
	colnames(merged)[c(3,4)] <- c('sample_name','fpkm')
	merged <- cbind(merged,fpkmNorm)

	# Clean up data frame for plotting; change gene_short_name vector to factor in
	# order to manipulate levels for plotting.  The initial rev(.) call plots
	# everything in alphabetical order on the y axis.
	merged$gene_short_name <- as.factor(merged$gene_short_name)

	# Reorder levels for plotting of y axis, in the order received from gsnList.
	if(apriori){
		if(!isGeneId(gsnList)){
			merged$gene_short_name <- factor(merged$gene_short_name,levels=gsnList)
		}else{
			merged$gene_short_name <- factor(merged$gene_short_name,levels=idToSym(gsnList))
		}
	}else{
		# Will get a warning for reorganising levels.  This is expected, because the melt
		# operation produces redundancy in the gene_short_name column, which R will complain
		# about but still render things just fine
		merged$gene_short_name <- idToSym(merged$gene_id)
		suppressWarnings(merged$gene_short_name <- factor(merged$gene_short_name,levels=merged$gene_short_name))
	}
	
	
	preTrunc <- as.data.frame(round(merged$fpkm))
	colnames(preTrunc) <- 'fpkmPre'
	merged <- cbind(merged,preTrunc)

	if (rMax>0){
		# Track elements that are truncated.
		isTrunc <- merged$fpkm>rMax
		merged$fpkm <- pmin(merged$fpkm,rMax-0.0001)
	}


	# Define ggplot object.
	heatGr <- ggplot(merged,aes_string(y='gene_short_name',x='sample_name'))

	# Define tile colouring.
	if (doNorm){
		heatGr <- heatGr + geom_tile(aes(fill=fpkmNorm),colour='white')
		theRangeMax <- 1
	}else{
		heatGr <- heatGr + geom_tile(aes(fill=fpkm),colour='white')
		theRangeMax <- max(merged$fpkm)
	}
	if (rMax>0){
		theRangeMax <- rMax+0.0001
	}

	# Define other stylistic features.
	heatGr <- heatGr + scale_fill_gradientn(colours=c('cyan','white','red'),
		values=c(0,0.5,1),guide='colorbar',limits=c(0,theRangeMax)) # Changed from green and purple.
	heatGr <- heatGr + theme(panel.background=element_blank())
	heatGr <- heatGr + theme(axis.text.x = element_text(angle=45,hjust=1))
	heatGr <- heatGr + labs(x='Cell type',y='Gene')

	if(rMax>0){
		if(label){
			heatGr <- heatGr + geom_text(data=subset(merged,merged$fpkmPre>rMax),aes_string(label='fpkmPre'),size=2)
		}
	}
	print(heatGr)

	invisible(heatGr)
}
