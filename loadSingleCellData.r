#####################################################################################
#
# Mark Cembrowski, Janelia Research Campus, Nov 11 2016
#
# Load single cell data from Zeisel et al., Science, 2015.
#
#####################################################################################

loadSingleCellData <- function(pop1='CA1Pyr1',pop2='CA1Pyr2') {
	print('Loading single-cell RNA-seq data...')

	# Load metadata.
	datLoc <- './data/expression_mRNA_17-Aug-2014.txt'
	scMeta <- read.table(datLoc,sep='\t',nrow=10,comment.char='') 
	rownames(scMeta) <- scMeta[,2]
	scMeta <- scMeta[,-c(1,2)]
	scMeta <- as.data.frame(t(scMeta))
	scMeta <<- scMeta

	# Load single cell data.
	scDataFull <- read.table(datLoc,sep='\t',skip=11,comment.char='')
	scNames <- scDataFull[,1]
	scDataFull <- scDataFull[,-c(1,2)]
	print('done loading data')

	# Remove duplicate entries.  These comprise just two entries, which won't affect the
	# overall analysis, nor are they critical loci.
	notDup <- !duplicated(scNames)
	scNames <- scNames[notDup]
	scDataFull <- scDataFull[notDup,]

	rownames(scDataFull) <- scNames

	# Filter according to CA1 PCs.
	keepCols1 <- scMeta$level2class == pop1
	keepCols2 <- scMeta$level2class == pop2
	keepCols <- keepCols1+keepCols2 > 0.1
	scData <- scDataFull[,keepCols]
	
	scCa1 <<- scData

	# Build another matrix corresponding to all hippocampus data.
	keepCols <- scMeta$tissue == 'ca1hippocampus'
	scData <- scDataFull[,keepCols]

	scHip <<- scData

	print('Done.')
}
