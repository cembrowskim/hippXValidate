##################################################################################
#
# Mark Cembrowski, Janelia Research Campus, Nov 10 2016
#
# This script loads the barcoded probe names from the gene sets of Shah et al.
#
# The 125 gene set is taken from mmc5 on the Neuron website, and the 249 gene set
# is taken from mmc6 on the Neuron website.
#
##################################################################################

loadGenes <- function(){
	# Load 125 gene set.	
	temp125 <- read.table('./data/barcoded125.txt',sep='\t',header=T)

	genesAll125 <- sort(as.character(temp125$gene))
	genesBar125 <- as.character(temp125$gene[temp125$barcode>0.1])

	# Capitalise first letter in genes.
	theFirst <- toupper(substr(genesAll125,1,1))
	theLast <- substr(genesAll125,2,nchar(genesAll125))
	genesAll125 <- paste(theFirst,theLast,sep='')
	genesAll125 <- sort(genesAll125)

	genesAll125 <<- genesAll125
	genesBar125 <<- genesBar125

	# Load 249 gene set.
	temp249 <- read.table('./data/barcoded249.txt',sep='\t',header=T)

	genesAll249 <- sort(as.character(temp249$gene))
	genesBar249 <- as.character(temp249$gene[temp249$barcode>0.1])

	# Capitalise first letter in genes.
	theFirst <- toupper(substr(genesAll249,1,1))
	theLast <- substr(genesAll249,2,nchar(genesAll249))
	genesAll249 <- paste(theFirst,theLast,sep='')
	genesAll249 <- sort(genesAll249)

	genesAll249 <<- genesAll249
	genesBar249 <<- genesBar249
}
	
