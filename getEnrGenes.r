#####################################################################################
#
# Mark Cembrowski, Janelia Research Campus, April 15 2015
#
# This script return sgenes that are enriched in one (or more) populations, relative
# to another set of one or more populations.  sampleX will be the enriched population,
# sampleY will be the depleted population.
#
#####################################################################################

getEnrGenes <- function(sampleX,sampleY,fdr=-1,foldThres=-1,fpkmThres=-1,avgPass=T){
	# Populate a matrix with the genes that fit supplied arguments.
	# Begin with fold change.
	if(foldThres>1){
		if(avgPass){
			enrXFold <- .apply(fpkmPoolMat[,sampleX],1,min)>
				(foldThres* (.apply(fpkmPoolMat[,sampleY],1,max)) )
		}else{
			minX <- .apply(fpkmRepMat[,grepl(paste(sampleX,collapse='|'),colnames(fpkmRepMat))],1,min)
			maxY <- .apply(fpkmRepMat[,grepl(paste(sampleY,collapse='|'),colnames(fpkmRepMat))],1,max)
			enrXFold <- minX>(foldThres*maxY)
	}
		enrXFold <- rownames(fpkmPoolMat)[enrXFold]
	}else{
		# Take an effective fold threshold of 1, preventing underexpressed
		# genes in sampleX from showing up if not screening for fold changes
		# (e.g., when using FDR to screen).  Otherwise, the FDR results
		# give bidirectional (ie, significantly up AND down regulated) genes
		# which is undesirable.
		enrXFold <- .apply(fpkmPoolMat[,sampleX],1,min)>
				( .apply(fpkmPoolMat[,sampleY],1,max) )
		enrXFold <- rownames(fpkmPoolMat)[enrXFold]
	}

	# Next, annotate all genes that cross threshold.
	if(fpkmThres>0){
		if(avgPass){
			enrXFpkm <- .apply(fpkmPoolMat[,sampleX],1,max)>fpkmThres
		}else{
			enrXFpkm <- fpkmRepMat[,grepl(paste(sampleX,collapse='|'),colnames(fpkmRepMat))]
			enrXFpkm <- .apply(enrXFpkm,1,min)>fpkmThres
			
		}
		enrXFpkm <- rownames(fpkmPoolMat)[enrXFpkm]
	}else{
		enrXFpkm <- rownames(fpkmPoolMat)
	}
	
	# Next, annotate all genes that are differentially expressed.
	enrXThres <- intersect(enrXFold,enrXFpkm)
	firstTest <- T
	if(fdr>0){
		for (ii in sampleX){
			for (jj in sampleY){
				# Recover pairwise DE tests.  May have a column name of sample1-sample2
				# OR sample2-sample1 ; cover both cases.
				theTest <- paste(ii,'.',jj,sep='')
				theCol <- which(colnames(qVal)%in%theTest)
		
				if(length(theCol)<0.1){
					theTest <- paste(jj,'.',ii,sep='')
					theCol <- which(colnames(qVal)%in%theTest)
				}
				enrXDeTemp <- qVal$gene_id[qVal[,theCol]<fdr]
				if(firstTest){
					enrXDe <- enrXDeTemp
					firstTest <- F
				}else{
					enrXDe <- intersect(enrXDe,enrXDeTemp)
				}
			}
		}
		enrX <- intersect(enrXThres,enrXDe)
		
	}else{
		enrX <- enrXThres
	}

	# Arrange things in a nice order, by fold change, for output.
	subEnr <- fpkmPoolMat[enrX,sampleX]
	subDep <- fpkmPoolMat[enrX,sampleY]
	foldEnrich <- .apply(subEnr,1,mean)/.apply(subDep,1,mean)

	enrX <- enrX[order(-foldEnrich)]

	invisible(enrX)
}

# This is a helper function that acts like apply when the argument X is a matrix/data
# frame, and just returns the supplied X when it is a vector.  This is to circumvent 
# issues when X can be either a vector or matrix/data frame when selecting an arbitrary
# number of columns from a data frame.
.apply <- function(X,MARGIN,FUN){
	if(length(dim(X))<0.1){
		return(X)
	}else{
		return(apply(X,MARGIN,FUN))
	}
}
