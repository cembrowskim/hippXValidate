#########################################################################################
#
# Mark Cembrowski, Janelia Research Campus, April 15 2015
#
# This script plots, in detail, the expression of a given gene.
#
# INPUT:
# theGene: character.  Symbol or ID of gene to plot.
# replicates=F: logical.  Visualise individual replicate FPKM values as data points.
# mask=c(): character.  Sample names to hide from final visualisation. 
# rMax=-1: numeric.  Maximum y value.
#
#########################################################################################

fpkmBarplot <- function(theGene,replicates=F,mask=c(),unmask=c(),rMax=-1){
	if(length(mask)*length(unmask)>0.1){
		stop('Cannot use mask and unmask options at same time.')
	}

	if(!isGeneId(theGene)){
		# Transform to Ensembl ID.
		theGene <- symToId(theGene)[[1]] # Successful hits.
	}
	
	# Get FPKM value for gene.
	df <- subset(fpkmPool,gene_id==theGene)
	
	# If requested replicates, get replicates.
	if(replicates){
		dfRep <- subset(fpkmRep,gene_id==theGene)
	}

	# Mask samples, if selected.
	if(length(mask)>0.1){
		toMask <- grep(paste(mask,collapse='|'),df$sample_name)
		if(length(toMask)<0.1){
			stop('All requested masked sample(s) not found.')
		}
		
		df <- df[-toMask,]

		if(replicates){
			toMask <- grep(paste(mask,collapse='|'),dfRep$sample_name)
			dfRep <- dfRep[-toMask,]
		}
	}	
	
	# Show only unmasked samples, if selected.  This operation is the opposite
	# of mask; ie, only the unmasked samples will be shown.
	if(length(unmask)>0.1){
		toUnmask <- grep(paste(unmask,collapse='|'),df$sample_name)
		if(length(toUnmask)<0.1){
			stop('All requested unmasked sample(s) not found.')
		}

		df <- df[toUnmask,]

		if(replicates){
			toUnmask <- grep(paste(unmask,collapse='|'),dfRep$sample_name)
			dfRep <- dfRep[toUnmask,]
		}
	}	
	# Assign ordering to specific groups.
	theLimits <- c('dg_d','dg_v','ca3_d','ca3_v','ca2','ca1_d','ca1_i','ca1_v','pv','sst')

	# Account for masked samples.
	if(length(mask)>0.1){
		for (ii in 1:length(mask)){
			curMask <- mask[ii]
			toMask <- theLimits==curMask
			theLimits <- theLimits[!toMask]
		}
	}
	
	# Account for unmasked samples.
	if(length(unmask)>0.1){
		theLimits <- unmask
	}

	# Plot.
	p <- ggplot(df,aes(x=sample_name,y=fpkm))
	p <- p + geom_bar(stat='identity',aes(fill=sample_name,colour=sample_name))
	if(replicates){
		p <- p + geom_point(aes(x=sample_name,y=fpkm),size=3,
			shape=18,colour='black',data=dfRep)
	}

	p <- p + geom_errorbar(aes(ymin=conf_lo,ymax=conf_hi,group=1),colour='black',width=0.5)

	if(rMax>0){
		p <- p + coord_cartesian(ylim=c(-rMax/100,rMax)) # -rMax/100 to avoid clipping.
	}

	p <- p + xlab('Cell type') + ylab('FPKM')
	p <- p + ggtitle(idToSym(theGene))
	p <- p + theme_bw() + theme(legend.position='none')
	p <- p + scale_x_discrete(limits=theLimits)

	print(p)
}

