###########################################################################################
#
# Mark Cembrowski, Janelia Farm, October 26 2013
#
# This script takes a list of a genes (either in XLOC form or gene_short_name form) and 
# returns the corresponding fpkm matrix.
#
# This script is analogous to the fpkmMatrix call implemented within cummeRbund, but offers
# extraction of FPKM values only for the genes of interest.  Replicates can be obtained
# if the replicates=T option is used; this is analogous to the repFpkmMatrix call from
# cummeRbund.
#
# INPUT: a list of genes to generate FPKM values for (either gene_id or gene_short_name;
# 	the script will recognise the form provided based upon the first entry).
# OPTIONAL INPUT: 
#	'replicates': obtain FPKM values for replicates? [default=F]
#	'outId': use gene_id as row names in output? [default=F; the default row values
#		are those supplied in geneList; ie if gene_ids are used, then gene_ids are
#		returned, and if gene_short_names are used then gene_short_names are
#		returned.]  Turning this on forces row names to be gene_id.
#	'outShortName': use gene_short_name as row names in output [default=F; see note
#		for outId optional argument above].
# OUTPUT: an fpkm matrix is returned.
#
############################################################################################

subFpkmMatrix <- function(geneList,replicates=F,outId=F,outShortName=F){
	# Implement a simple check to make sure that both outId and outShortName are not
	# on at the same time.
	if( (outId+outShortName)>1.1 ){
		print('Cannot use outId=T and outShortName=T at the same time.  This is akin')
		print('to trying to force row names on output matrix to be both gene_ids and')
		print('gene_short_name at the same time.  Set at least one of these values F.')
		return()
	}

	# Generate whole FPKM matrix.
	if(replicates){
		theMat <- fpkmRepMat
	}else{
		theMat <- fpkmPoolMat
	}

	# Detect if geneList is gene_id or gene_short_name
	if(isGeneId(geneList)){
		theIds <- geneList
		idList <- T
	}else{
		# Find XLOCs corresponding to each gene_short_name
		symTrans <- symToId(geneList)
		theIds <- symTrans[[1]]
		missed <- symTrans[[2]]
		idList <- F
	}
	
	# Extract only the desired rows.
	theMat <- theMat[theIds,]

	# Reassign row names if supplied a list of gene_short_names
	if(!idList){
		if (!outId){
			# No override if outId=F
			rownames(theMat) <- geneList[!geneList%in%missed] # Remove misses.
		}
	}

	# Override row names if outShortName=T.
	if( outShortName && idList){
		rownames(theMat) <- idToSym(geneList)
	}

	invisible(theMat)
}
