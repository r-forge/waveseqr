# Uses two-sided binomial test to compare test sample to control
# and corrects p-values for multiple-testing

controlFDR <- function(peaks, chip, control, outfile, winsize=200,
					minreads=8, binom.sided="two.sided", adj.method="fdr"){		
	cat("\nEstimating statistical significance in presence of control ...\n")		
	
	starttime <- proc.time()[3]
	# compare chip with control	
	control.mat <- compareWithControl(peaks, chip, control, winsize, minreads, binom.sided)
	
	progress <- paste("Total time elapsed = ",as.character(proc.time()[3] - starttime)," seconds\n",sep="")
    cat(progress)
	
	# adjust p-values for multiple testing	
	cat("\nAdjusting p-values for multiple testing ...\n")
	p.val <- control.mat[,7]	
	p.adj.mat <- p.adjust(as.numeric(p.val),adj.method)
	p.val.mat <- cbind(as.character(p.val),as.character(p.adj.mat))
		
	cat("\nPrinting peaks to file ...\n")
	write.table(cbind(as.character(control.mat[,1]),control.mat[,2:6],p.val.mat),
			outfile,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
	
}