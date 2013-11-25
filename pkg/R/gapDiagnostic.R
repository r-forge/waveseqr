gapDiagnostic <- function(peaks, gap.max=15, topN=10000, winsize=200, plot.flag=TRUE){
	# if file
	if(!is.data.frame(peaks)){
		peaks <- read.table(peaks, sep="\t", header=FALSE, stringsAsFactors=FALSE)
	}
	
	peaks <- peaks[order(peaks[,1],peaks[,2]),]
	total.reads <- sum(peaks[,4])	
	if(dim(peaks)[1] < topN){
		topN <- dim(peaks)[1]
		cat("Adjusting peak number to ",topN,"\n")
	}
	
	peaks.ordered <- peaks[order(peaks[,4], decreasing=TRUE),]
	peak.frac <- sum(peaks.ordered[1:topN,4])/total.reads
	#cat(total.reads,"\t",peak.frac,"\n")
	#cat(unique(peaks[,1]),"\n")
	for(g in 1:gap.max){
		capture.output( 
			temp.peaks <- gapPeaks(peaks, outfile=NA, gap=g, winsize=winsize) 
		)
		#cat(dim(temp.peaks),"\n")
		if(dim(temp.peaks)[1] < topN){
			peak.num <- dim(temp.peaks)[1]
		} else {
			peak.num <- topN
		}
		peaks.ordered <- temp.peaks[order(temp.peaks[,4], decreasing=TRUE),]
		peak.frac <- c(peak.frac, sum(peaks.ordered[1:peak.num,4])/total.reads)
	}
	
	if(plot.flag){
		#dev.new()
		plot(0:gap.max, peak.frac, type="o",pch=2, col="red",
			xlim=c(0,gap.max), ylim=c(0,1),
			xlab="Gap Size", ylab="Fraction of total reads")
	}
	return(peak.frac)
}