###################################################################
#
# This script looks at the number of cells that have expression
# of a given molecule below a given threshold.
#
###################################################################

singleCellBelowThres <- function(geneList,theThres=1){

	# Load signel cell data.
	source('~/Dropbox (HHMI)/research/reviews/shah2016/loadSingleCellAllHipp.r')
	
	loadSingleCellHipp()

	toKeep <- geneList[geneList%in%rownames(scHipp)]
	tempCount <- scHipp[toKeep,]
print(head(tempCount))
	tempThres <- apply(tempCount,1,max)<theThres

	print(paste('Of',nrow(tempCount),'cells,',
		sum(tempThres),'have a count less than',theThres))
}
