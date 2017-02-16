###################################################################
#
# Mark Cembrowski, Janelia Research Campus, Nov 11 2016
#
# This single-serving script plots the results of looking for
# barcoded genes across all hippocampal populations.
#
###################################################################

geneBar125Heatmap <- function(label=F){
	# Identify probes that are part of pop'n RNA-seq dataset.
	probesInPop <- genesBar125%in%unlist(idToSym(rownames(fpkmPoolMat)))

	# Extract these genes.
	probesInPop <- genesBar125[probesInPop]

	# Get maximum FPKM values for these genes.
	maxFpkm <- apply(subFpkmMatrix(probesInPop),1,max)

	# Sort by FPKM.
	maxSort <- probesInPop[order(maxFpkm)]

	# Plot.
	fpkmHeatmap(maxSort,rMax=10,label=label)
}
