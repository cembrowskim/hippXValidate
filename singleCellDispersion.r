##########################################################################
#
# Mark Cembrowski, Janelia Research Campus, Nov 11 2016
#
# Graph the dispersion of single-cell RNA-seq data.
#
##########################################################################

singleCellDispersion <- function(xMin=0.1,yMin=0.1,yMax=-1,gsnList=c(),
		doLabel=F,sdEnr=2,doHl=T,doBin=T,doFit=F,quiet=T){
	# Generate dataset contaning means and CVs.
	df <- scCa1
	theMean <- apply(df,1,mean)
	theCv <- apply(df,1,sd)/theMean

	df <- cbind(df,theMean,theCv)
	colnames(df)[length(colnames(df))-1] <- 'mean'
	colnames(df)[length(colnames(df))] <- 'cv'

	# Remove transcripts that have zero expression.
	df <- df[!is.na(df$cv),]

	# Bin results with a variable bin size.
	numElem <- 100 # Number of transcripts for each bin.
	df <- df[order(df$mean),] # Sort according to mean.
	binMean <- c()
	binCv <- c()
	binErr <- c()
	binAnnot <- c()

	# Generate gsnList.
	gsnList <- rownames(scCa1)[rownames(scCa1)%in%genesBar125]

	# Print out a quick summary.
	if(!quiet){
		zeiselSub <- scCa1[rownames(scCa1)%in%genesBar125,]
		meanSummary <- apply(zeiselSub,1,mean)
		meanSummary <- meanSummary[order(meanSummary)]
		print(meanSummary)
		print('Number of genes undetected in Zeisel:')
		print(sum(meanSummary<0.00000001))
		sumSummary <- apply(zeiselSub,1,sum)
		sumSummary <- sumSummary[order(sumSummary)]
		print(sumSummary)	
	}

	# Track results of gene short list, if supplied.
	gsnData <- as.data.frame(matrix(nrow=length(gsnList),ncol=7))
	colnames(gsnData) <- c('binMean','binSd','binMed','binMad',
				'geneCv','geneSd','geneMad')
	geneInd <- 1

	if(doBin){
		leftInd <- 1
		repeat{
			rightInd <- leftInd+numElem-1
			if(rightInd>nrow(df)){
				rightInd <- nrow(df)
			}
	
			curNames <- rownames(df)[leftInd:rightInd]
			curX <- df$mean[leftInd:rightInd]
			curY <- df$cv[leftInd:rightInd]
	
			binMean[length(binMean)+1] <- mean(curX)
			binCv[length(binMean)+1] <- median(curY)#mean(curY)
			binErr[length(binMean)+1] <- mad(curY)#sd(curY)
	
			toAdd <- intersect(curNames,gsnList)
			if(length(toAdd)>0.1){
				for (ii in 1:length(toAdd)){
					rownames(gsnData)[geneInd] <- toAdd[ii]
					gsnData[geneInd,'binMean'] <- mean(curY)
					gsnData[geneInd,'binSd'] <- sd(curY)
					gsnData[geneInd,'binMed'] <- median(curY)
					gsnData[geneInd,'binMad'] <- mad(curY)
					
					gsnData[geneInd,'geneCv'] <- df[toAdd[ii],]$cv
					gsnData[geneInd,'geneSd'] <- (df[toAdd[ii],]$cv-mean(curY))/sd(curY)
					gsnData[geneInd,'geneMad'] <- (df[toAdd[ii],]$cv-median(curY))/mad(curY)
			
					geneInd <- geneInd + 1
					
				}
			}
			if(sum(curY>(mean(curY)+sdEnr*sd(curY)))>0.1){
				binAnnot <- c(binAnnot,curNames[curY>(median(curY)+sdEnr*mad(curY))])
			}

			if(abs(rightInd-nrow(df))<0.1){
				break
			}
	
			leftInd <- leftInd + numElem
		}
		binLow <- binCv-2*binErr
		binHigh <- binCv+2*binErr
		binVals <- as.data.frame(cbind(binMean,binErr,binCv,binLow,binHigh))
	}

	if(length(gsnList)<0.1){
		annot <- rownames(df)%in%binAnnot
	}else{
		geneId <- gsnList
		#geneId <- symToId(gsnList) # assume gsn
		annot <- rownames(df)%in%geneId
	}	
	annotFlag <- annot
	annot[annotFlag] <- rownames(df)[annotFlag]

	annot[!annotFlag] <- ''
	
	df <- cbind(df,annot)

	# Scatter the remainder.
	gg <- ggplot(df,aes(x=mean,y=cv)) 

	gg <- gg + theme_bw()

	if(yMax<0){ yMax <- max(df$cv)*1.1 }

	gg <- gg + scale_x_log10(limits=c(xMin,max(df$mean)*1.1)) + scale_y_log10(limits=c(yMin,yMax))

	if(doFit){
		theLineX <- seq(0.01,1000,0.01)
		theLineP <- 1/sqrt(theLineX) # Poisson.
		theLineF <- theLineX^(-0.55)+0.64 # Fit line from Zeisel et al
		dfLine <- data.frame(x=theLineX,p=theLineP,f=theLineF)
		gg <- gg + geom_line(data=dfLine,aes(x=x,y=p),colour='red')
		gg <- gg + geom_line(data=dfLine,aes(x=x,y=f),colour='cyan')
	}
	gg <- gg + geom_point(alpha=0.1) 	

	
	if(doBin){
		gg <- gg + geom_line(data=binVals,aes(x=binMean,y=binCv),colour='green',alpha=0.5)
		gg <- gg + geom_line(data=binVals,aes(x=binMean,y=binLow),colour='green',alpha=0.5)
		gg <- gg + geom_line(data=binVals,aes(x=binMean,y=binHigh),colour='green',alpha=0.5)
	}
	
	

	if(doLabel){	
		gg <- gg + geom_text(aes(label=annot),hjust=0,vjust=0,size=2,colour='red')
	}
	if(doHl){
		gg <- gg + geom_point(data=subset(df,annot!=''),colour='red',alpha=1.0)
	}

	print(gg)

	invisible(gsnData)
}

