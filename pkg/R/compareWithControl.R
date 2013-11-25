compareWithControl <- function(peaks, chip, control, winsize=200, minreads=8, binom.sided="two.sided"){
	# This function reads peak, chip and control files, applies 
	# two-sided binomial test to compare chip with control
	# and outputs peaks with p-values
	
	cat("\nComparing chip with control ...\n")
	#input file with expected format: chr start end reads	
	cat("\tReading peaks file ...\n")
	peaks.file <- read.table(file=peaks,sep="\t",header=FALSE,stringsAsFactors=FALSE)

	peaks.read.ind <- peaks.file$V4 >= minreads
	peaks.file <- peaks.file[peaks.read.ind,]

	# file expected to be in bedGraph format: chr start end reads
	peaklen <- peaks.file$V3 - peaks.file$V2
	peaklen <- peaklen + 1
	peak.win <- peaklen/winsize

	# peak chip reads
	chip.peak.reads <- peaks.file$V4

	# peak locations
	peak.loc <- peaks.file$V2

	# read padded chip file	
	cat("\tReading chip file ...\n")
	chip.file <- read.table(file=chip,sep="\t",header=FALSE,stringsAsFactors=FALSE)
	chip.reads <- chip.file$V4
	total.chip.reads <- sum(chip.reads)

	# get chip reads in peak	
	cat("\tGetting reads in peaks ...\n")
	control.peak.reads <- getPeakReads(peaks.file, chip.file, winsize)

	# input padded control file	
	cat("\tReading control file ...\n")
	control.file <- read.table(file=control,sep="\t",header=FALSE,stringsAsFactors=FALSE)
	control.reads <- control.file$V4
	total.control.reads <- sum(control.reads)

	# get control reads in peak	
	cat("\tGetting control reads in peaks ...\n")
	control.peak.reads <- getPeakReads(peaks.file, control.file, winsize)
	
	cat("\tNormalizing reads to per million mapped reads ...\n")
	# normalize chip and chip peak reads
	chip.reads <- chip.reads*1e6/total.chip.reads
	chip.peak.reads <- chip.peak.reads*1e6/total.chip.reads
	peaks.file$V4 <- chip.peak.reads*1e6/total.chip.reads

	# normalize control reads
	control.reads <- control.reads*1e6/total.control.reads
	control.peak.reads <- control.peak.reads*1e6/total.control.reads
	control.file$V4 <- control.reads

	# matrix to contain chip and control reads
	reads.chip.control <- cbind(chip.peak.reads,control.peak.reads)

	# variable to contain fold-changes
	peak.ratio.all <- (chip.peak.reads + 1)/(control.peak.reads + 1)

	# use binomial test to compare chip and control data	
	cat("\tComparing chip and control reads ...\n")
	peak.num <- length(chip.peak.reads)
	peak.pval <- mat.or.vec(peak.num,1)
	for(j in 1:peak.num){
		test.result <- binom.test(floor(reads.chip.control[j,]),alternative=binom.sided)
		peak.pval[j] <- test.result$p.value
	}

	cbind(as.character(peaks.file$V1),peaks.file$V2,peaks.file$V3,
		chip.peak.reads, control.peak.reads, peak.ratio.all, peak.pval)
}
