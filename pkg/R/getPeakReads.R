# function to count reads in islands
getPeakReads <- function(peak.mat, reads.mat, winsize=200){
	# get peak lengths and locations
	peak.len <- peak.mat[,3] - peak.mat[,2]
	peak.len <- peak.len + 1
	peak.win <- peak.len/winsize

	# get peak locations & correct for 0-based indexing
	peak.loc <- peak.mat[,2]
	peak.loc.index <- peak.loc/winsize
	peak.loc.index <- peak.loc.index + 1

	# get reads
	reads <- reads.mat[,4]
	
	chr.all <- levels(peak.mat[,1])
	
	# vector to contain read counts
	peaknum <- length(peak.loc)
	read.counts.all <- mat.or.vec(peaknum,1)

	for(i in 1:length(chr.all)){
		chr <- chr.all[i]

		# do for each chromosome
		chr.index <- peak.mat[,1] == chr
		
		if(sum(chr.index) > 0){
			# get necessary info for chromsome
			peak.loc.index.chr <- peak.loc.index[chr.index]
			peak.win.chr <- peak.win[chr.index]	
			peak.num.chr <- length(peak.loc.index.chr)		
	
			# reads for chromosome
			reads.chr.index <- reads.mat[,1] == chr
			reads.chr <- reads[reads.chr.index]
	
			# zero vector for counts
			read.counts.chr <- mat.or.vec(peak.num.chr,1)
	
			for(j in 1:peak.num.chr){
				temp.loc <- peak.loc.index.chr[j]
				temp.win <- peak.win.chr[j]
				if((temp.loc + (temp.win-1)) <= length(reads.chr)){
					read.counts.chr[j] <- sum(reads.chr[temp.loc:(temp.loc + (temp.win-1))])			
				} else if(temp.loc <= length(reads.chr) && (temp.loc + (temp.win-1)) > length(reads.chr)){
					read.counts.chr[j] <- sum(reads.chr[temp.loc:length(reads.chr)])
				} else {
					read.counts.chr[j] <- 0
				}
				
			}			
			read.counts.all[chr.index] <- read.counts.chr
		}		
	}	
	read.counts.all
}