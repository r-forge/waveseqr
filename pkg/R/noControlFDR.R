# Uses randomized algorithm to calculate ecdf of random background
# calculate p-values and corrects for multiple-testing

noControlFDR <- function(peaks, chip, outfile, winsize=200,
					minreads=8,samplesize=1e6, adj.method="fdr"){	
	cat("\nEstimating statistical significance in the absence of control ...\n")

	starttime <- proc.time()[3]
	
	# obtain ecdf from randomizedBackground
	background.mat <- randomizedBackground(peaks,chip,samplesize,winsize)

	progress <- paste("Total time elapsed = ",as.character(proc.time()[3] - starttime)," seconds\n",sep="")
    cat(progress)
	
	# remove peaks with reads < minreads
	read.ind <- as.numeric(background.mat[,4]) >= minreads
	background.mat.filt <- background.mat[read.ind,]
	p.val.filt <- background.mat.filt[,5]

	# adjust p-values for multiple-testing	
	cat("\nAdjusting p-values for multiple testing ...\n")	
	p.adj.filt <- p.adjust(as.numeric(p.val.filt),adj.method)
	p.val.mat <- cbind(as.character(p.val.filt),as.character(p.adj.filt))
	
	cat("\nPrinting peaks to file ...\n")
	write.table(cbind(as.character(background.mat.filt[,1]),background.mat.filt[,2:4],p.val.mat),
		outfile,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)	
}