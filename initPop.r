##############################################################################################
#
# Mark Cembrowski, Janelia Research Campus, 
#
# This script loads up the relevant tables and scripts for population RNA-seq. 
#
##############################################################################################

# Load scripts.
	source('./code/loadPopData.r')		# load data
	source('./code/dendro.r')		# dendrograms
	source('./code/fpkmBarplot.r') 		# single gene barplots
	source('./code/fpkmHeatmap.r') 		# multi gene heatmaps
	source('./code/fpkmScatter.r') 		# two transcriptomes scatter plot
	source('./code/geneBar125Heatmap.r')	# heatmap of genes across hipp
	source('./code/idToSym.r')		# convert gene IDs to symbols
	source('./code/isGeneId.r')		# check to see if a supplied ID is a legit gene ID
	source('./code/maskSamples.r')		# mask samples from analysis
	source('./code/screenSamples.r')	# make sure supplied sample names are in dataset.
	source('./code/subFpkmMatrix.r')	# allow easy obtaining of gene expression subsets
	source('./code/symToId.r')		# convert gene symbols to IDs

# Load libraries.
	library(ggplot2)
	library(reshape2)

