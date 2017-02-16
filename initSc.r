##############################################################################################
#
# Mark Cembrowski, Janelia Research Campus, 
#
# This script loads up the relevant tables and scripts for single-cell RNA-seq. 
#
##############################################################################################

# Load scripts.
	source('./code/loadSingleCellData.r') # Load single cell data.
	source('./code/singleCellBarplot.r')	# single cell barplot
	source('./code/singleCellPca.r') # pca
	source('./code/singleCellDispersion.r')	# single cell dispersion
	source('./code/singleCellMultigene.r') # single cell multigene plot
