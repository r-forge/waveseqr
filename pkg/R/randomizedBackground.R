randomizedBackground <- function(peaks, chip, samplesize=1e6, winsize=200){
# This function reads peak and read files, applies 
# randomized algorithm to simulate background distribution 
# obtains empirical cdf and outputs peaks with p-values

cat("\nRunning randomized background estimation ...\n")
#input file with expected format: chr start end reads
cat("\tReading peaks file ...\n")
infile <- read.table(file=peaks,sep="\t",header=FALSE, stringsAsFactors=FALSE)

# file in BedGraph format: chr start end reads
peaklen <- infile$V3 - infile$V2
peak.reads <- infile$V4

# read padded bedGraph file
cat("\tReading padded graph file ...\n")
chip.file <- read.table(file=chip,sep="\t",header=FALSE, stringsAsFactors=FALSE)

# variable to contain reads for randomized peaks
sample.reads.all <- 0

# list of chromosomes
chr.all <- levels(infile$V1)

# chromosome sizes
chr.all.sizes <- table(infile$V1)
chr.size.max <- max(chr.all.sizes)

cat("\tRandomizing peak locations ...\n")
# do for each chromosome
for(j in 1:length(chr.all)){
	chr <- chr.all[j]
	
	# print progress	
	cat("\t\t",chr,"\n")

	# get locations of reads on chromosome
	chr.index <- infile$V1 == chr
	chr.size <- sum(chr.index)
	peak.num <- floor(samplesize*chr.size/chr.size.max)

	# peak lengths for chr
	peak.chr <- peaklen[chr.index]

	# read locations for chr
	chip.chr.index <- chip.file$V1 == chr
	chip.chr <- chip.file$V4[chip.chr.index]
	chip.loc.chr <- chip.file$V2[chip.chr.index]
	chip.loc.start <- chip.loc.chr[1]
	
	#sample from peak length distribution
	peak.sample <- sample(peak.chr,peak.num,replace=TRUE)

	# correct for 0-based indexing and get number of windows in each peak
	peak.sample <- peak.sample+1
	peak.window <- peak.sample/winsize
	
	# randomize island locations
	peak.loc.sample <- sample(chip.loc.chr,peak.num,replace=TRUE)
	peak.loc.sample.index <- peak.loc.sample/winsize - chip.loc.start/winsize + 1;
	
	#create zero vector for reads
	sample.reads <- mat.or.vec(peak.num,1)
	
	for(i in 1:peak.num){
		temp.loc <- peak.loc.sample.index[i]
		temp.win <- peak.window[i]
		sample.reads[i] <- sum(chip.chr[temp.loc:temp.loc+temp.win])
	}
	
	if(j == 1){
		sample.reads.all <- sample.reads[sample.reads != 0]
	}
	else{
		sample.reads.all <- c(sample.reads.all,sample.reads[sample.reads != 0])
	}
}

x <- ecdf(sample.reads.all[sample.reads.all != 0])

pval <- x(peak.reads)
pval <- 1 - pval
cbind(as.character(infile$V1),infile$V2,infile$V3,peak.reads,pval)
}