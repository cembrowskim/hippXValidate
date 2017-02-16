##############################################################################
#
# Mark Cembrowski, Janelia Research Campus, Apr 8 2015
#
# This makes a dendrogram, and allows for masking of samples; it is initially
# based off of the cummeRbund implementation.
#
##############################################################################

dendro <- function (logMode=T,pseudocount=1,replicates=F,theMask=c(),theMethod='complete'){
	if(replicates){
		fpkmMat <- fpkmRepMat
	}
	else{
		fpkmMat <- fpkmPoolMat
	}
	if(logMode){
		fpkmMat <- log10(fpkmMat + pseudocount)
	}

	if(length(theMask)>0.1){
		toMask <- grepl(paste(theMask,collapse='|'),colnames(fpkmMat))
		fpkmMat <- fpkmMat[,!toMask]
	}

	res <- .JSdist(.makeprobs(fpkmMat))
	res <- as.dendrogram(hclust(res,method=theMethod))
	plot(res, main = paste("All", deparse(substitute(genes(cuff_data))),sep = " "))
}

.JSdist <- function(mat){
    res <- matrix(0, ncol = dim(mat)[2], nrow = dim(mat)[2])
    col_js <- apply(mat, MARGIN = 2, .shannon.entropy)
    colnames(res) <- colnames(mat)
    rownames(res) <- colnames(mat)
    for (i in 1:dim(mat)[2]) {
        for (j in i:dim(mat)[2]) {
            a <- mat[, i]
            b <- mat[, j]
            JSdiv <- .shannon.entropy((a + b)/2) - (col_js[i] + 
                col_js[j]) * 0.5
            res[i, j] = sqrt(JSdiv)
            res[j, i] = sqrt(JSdiv)
        }
    }
    res <- as.dist(res)
    attr(res, "method") <- "JSdist"
    res
}

.makeprobs <- function(a){
    colSums <- apply(a, 2, sum)
    b <- t(t(a)/colSums)
    b[is.na(b)] = 0
    b
}

.shannon.entropy <- function(p){
    if (min(p) < 0 || sum(p) <= 0) 
        return(Inf)
    p.norm <- p[p > 0]/sum(p)
    -sum(log2(p.norm) * p.norm)
}
