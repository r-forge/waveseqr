getPeakReads_batch <- function(peaks, indir, design, outfile){	
	# get file list
	targets <- read.table(design, sep="\t", header=FALSE, stringsAsFactors=FALSE)
	
	# get peak list
	peak.mat <- read.table(peaks, sep="\t", header=FALSE)

	# file list
	file.list <- targets[,1]

	for(i in 1:length(file.list)){	
		reads <- paste(indir,file.list[i],sep="/")
		
		sample.name <- targets[i,2]
		cat("Processing ", file.list[i], "\n")
		
		reads.mat <- read.table(reads, sep="\t", header=FALSE)
		peak.reads <- getPeakReads(peak.mat, reads.mat)
		
		if(i == 1){
			peak.reads.mat <- cbind(peak.mat[,1:3], peak.reads)
			colnames.all <- c("chr","start","end", sample.name)
		} else {
			peak.reads.mat <- cbind(peak.reads.mat, peak.reads)
			colnames.all <- c(colnames.all, sample.name)
		}
	}

	colnames(peak.reads.mat) <- colnames.all
	
	write.table(peak.reads.mat, file=outfile, quote=FALSE, row.names=FALSE, sep="\t")
}