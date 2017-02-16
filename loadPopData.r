##########################################################################################
#
# Mark Cembrowski, Janelia Research Campus, Nov 11 2016
#
# This script loads the FPKM data from population-level RNA-seq.
#
##########################################################################################

# Load FPKM values and differential expression results.
loadPopData <- function(){
	print('Loading population RNA-seq data...')
	
	# Load the proper datasets.
	fpkmPool <<- read.table('./data/fpkmPool.txt',sep='\t',header=T)
	fpkmRep <<- read.table('./data/fpkmRep.txt',sep='\t',header=T)
	fpkmPoolMat <<- read.table('./data/fpkmPoolMat.txt',sep='\t',header=T,row.names=1)
	fpkmRepMat <<- read.table('./data/fpkmRepMat.txt',sep='\t',header=T,row.names=1)

	print('Done.')
}

