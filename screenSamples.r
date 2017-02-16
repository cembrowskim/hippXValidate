#############################################################################
#
# Mark Cembrowski, Janelia Farm, January 29 2014
#
# This is a helper function that checks to make sure a list of sample names
# provided are all contained within the samples in the dataset.
#
# Search, as currently implement, is EXACT.  e.g., if the samples are
# c('ca3_mpp3','ca3_grp'), then providedNames='ca3' will be stopped, but
# 'providedNames='ca3_mpp3' will not be stopped.
#
# fpkmMat can be supplied if the full fpkmMatrix (or repFpkmMatrix) should
# not be used, or should not be regenerated (saves time for large
# matrices, especially when dealing with replicates).
#
# It halts execution if at least one sample is not found.
#
# INPUTS:
# providedNames: character.  List of samples to check for.
# 
# OPTIONAL ARGUMENTS:
# replicates=F: logical.  Use replicates in analysis.
# fpkmMat=matrix(): matrix.  Can supply a specific matrix for analysis.
#############################################################################

.screenSamples <- function(providedNames,replicates=F,fpkmMat=matrix()){

	if(prod(dim(fpkmMat))>1.1){
		theMat <- fpkmMat
		if(!replicates){
			namesInDb <- colnames(theMat)
		}else{
			namesInDb <- substr(colnames(theMat),1,nchar(colnames(theMat))-2)
		}
	}else{
		if(!replicates){
			theMat <- fpkmMatrix(genes(cuff_data))
			namesInDb <- colnames(theMat)
		}else{
			theMat <- repFpkmMatrix(genes(cuff_data))
			namesInDb <- colnames(theMat)
			# Remove replicate identifiers from the end of 
			# column names; e.g., turn 'ca3_mpp3_0' into
			# 'ca3_mpp3'
			namesInDb <- substr(namesInDb,1,nchar(namesInDb)-2)
			namesInDb <- unique(namesInDb)
		}
	}


	allInDb <- prod(providedNames%in%namesInDb)
	if(allInDb<0.1){
		stop('At least one provided name not in database.')
	}

	# Invisibly return full matrix, if going to be used by a downstream
	# script (e.g., masking).
	invisible(theMat)
}
