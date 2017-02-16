#############################################################################
#
# Mark Cembrowski, Janelia Research Campus, Nov 6 2016
#
# This script plots does PCA for all hippocampal cells.
#
# If gsnList is supplied, then the PCA is performed on this list; otherwise,
# high CV genes are chosen across a range of expression values.
#
#############################################################################

singleCellPca <- function(gsnList=c(),x='PC1',y='PC2',quiet=T,
				xLim=c(),yLim=c(),alpha=0.2){
	subMat <- scHip

	if(length(gsnList)<0.1){

		# Identify genes to use for analysis via overdispersion.
		theAvg <- apply(subMat,1,mean)
		binMin <- 1 # molecules/cell to consider use for analysis.
		subMat <- subMat[theAvg>binMin,]
		theAvg <- apply(subMat,1,mean)
		theSd <- apply(subMat,1,sd)
		theCv <- theSd/theAvg
		binSize <- 118
		totGenesPerBin <- 3

		subMat <- subMat[order(theAvg),]
		theAvg <- theAvg[order(theAvg)]
		theSd <- theSd[order(theAvg)]
		theCv <- theSd[order(theAvg)]

		genesToKeep <- c()
		for (ii in 1:(floor(length(theAvg))/binSize)){
			firstInd <- 1+(ii-1)*binSize
			lastInd <- ii*binSize
				
			tempMat <- subMat[firstInd:lastInd,]
			tempCv <- apply(tempMat,1,sd)/apply(tempMat,1,mean)
			tempCv <- tempCv[order(-1*tempCv)]
		
			genesToKeep <- c(genesToKeep,names(tempCv)[1:totGenesPerBin])
		}	
	}else{
		genesToKeep <- gsnList
		# Ensure all genes provided are in single cell data.
		genesToKeep <- genesToKeep[genesToKeep%in%rownames(scHip)]
	}
	print(paste('Using',length(genesToKeep),'for PCA'))	
	
	# Do PCA.
	tempDat <- scHip[genesToKeep,]

	# Need to exclude genes if not expressed anywhere.
	isZero <- apply(tempDat,1,sum) < 0.1
	tempDat <- tempDat[!isZero,]
#	tempDat <- log10(tempDat+1)

	colourRule <- c('astrocytes_ependymal'='blue',
			'interneurons'='purple','microglia'='green',
			'oligodendrocytes'='brown','endothelial-mural'='red','pyramidal CA1'='grey')
	
	PC <- prcomp(t(tempDat),scale=T)
	if(!quiet){print(summary(PC))}
	dat <- data.frame(obsnames=row.names(PC$x),PC$x)	
	theClasses <- scMeta$level1class[colnames(scHip)]
	dat[,'grp'] <- theClasses[rownames(dat)]

	# Sort data frame for nice plotting.
	dat <- rbind(
		subset(dat,grp=='pyramidal CA1'),
		subset(dat,grp=='oligodendrocytes'),
		subset(dat,grp=='interneurons'),
		subset(dat,grp=='astrocytes_ependymal'),
		subset(dat,grp=='microglia'),
		subset(dat,grp=='endothelial-mural'))
#	dat <- dat[order(-dat$grp),]

#
#	dat$grp <- factor(dat$grp,levels=c('astrocytes_ependymal',
#		'interneurons','oligodendrocytes','endothelial-mural','microglia','pyramidal CA1'))

	plot <- ggplot(dat,aes_string(x=x,y=y))
	plot <- plot + geom_point(alpha=alpha, size = 2, aes(label = obsnames,colour=grp))
	plot <- plot + theme_bw()
	plot <- plot + scale_colour_manual(values=colourRule)
#	plot <- plot + scale_x_log10() + scale_y_log10()

	
	if(abs(length(xLim)-2)<0.1){plot <- plot + scale_x_continuous(lim=xLim)}
	if(abs(length(yLim)-2)<0.1){plot <- plot + scale_y_continuous(lim=yLim)}
#	plot <- plot + scale_x_continuous(lim=c(xLo,xHi)) 
#	plot <- plot + scale_y_continuous(lim=c(yLo,yHi))

#	if(xLabLo<0){
#		plot <- plot + geom_text(data=subset(datLab,xVal<xLabLo),aes(x=xVal,y=yVal,label=obsnames),hjust=0,vjust=0,size=2)
#	}
#	if(xLabHi>0){
#		plot <- plot + geom_text(data=subset(datLab,xVal>xLabHi),aes(x=xVal,y=yVal,label=obsnames),hjust=0,vjust=0,size=2)
#	}
#	if(yLabLo<0){
#		plot <- plot + geom_text(data=subset(datLab,yVal<yLabLo),aes(x=xVal,y=yVal,label=obsnames),hjust=0,vjust=0,size=2)
#	}
#	if(yLabHi>0){
#		plot <- plot + geom_text(data=subset(datLab,yVal>yLabHi),aes(x=xVal,y=yVal,label=obsnames),hjust=0,vjust=0,size=2)
#	}

	print(plot)
}
