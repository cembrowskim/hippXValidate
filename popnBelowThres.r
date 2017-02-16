#####################################################################
#
# This script looks at a given list of genes (usually, the barcoded
# component of the 125 or 249 gene sets from Shah et al) and the
# number of genes that are across threshold.
#
#####################################################################

popnBelowThres <- function(geneList,theThres=10,thePops=c('dorsal','ventral',
			'intermediate','dg_d','dg_v','ca3_d','ca3_v','pv',
			'sst')){

	# Retrieve FPKM data.
	tempMat <- subFpkmMatrix(geneList)
	tempMat <- tempMat[,thePops]

	tempMax <- apply(tempMat,1,max)
	belowThres <- sum(tempMax < theThres)

	print(paste('For',length(tempMax),'genes, the total number below a threshold of',theThres,'is',belowThres))
}
