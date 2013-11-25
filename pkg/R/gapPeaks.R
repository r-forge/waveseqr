# input ungapped peaks file
# gap size, window size and output gapped peaks

gapPeaks <- function(peaks,outfile=NA,gap,winsize=200){

w <- options("warn")
options(warn=-1)
on.exit(options(w))
cat("\nConcatenating peaks separated by at most",as.character(gap),"non-significant windows ...\n")
starttime <- proc.time()[3]

#input file with expected format: chr start end reads p-score
if(!is.data.frame(peaks)){
	infile <- read.table(file=peaks,sep="\t",header=FALSE,stringsAsFactors=FALSE)
} else {
	infile <- as.data.frame(peaks)
}

# file expected to be in BedGraph format: chr start end reads
# get peak starts and ends
start <- infile[,2]
end <- infile[,3]
reads <- infile[,4]

if(is.factor(infile[,1])){
	infile[,1] <- as.character(infile[,1])
}
chr.all <- unique(infile[,1])
gap.size <- gap*winsize

chr.mat.all <- 0
cat("\tProcessing chromosome-by-chromosome ...\n")
#cat(chr.all)
#cat(head(infile),"\n")
#print(head(infile))
for(j in 1:length(chr.all)){
	chr <- chr.all[j]

	# print progress to screen
	cat("\t\t",chr,"\n")
	
	#analyze chromosome by chromosome
	chr.index <- infile[,1] == chr
	start.chr <- start[chr.index]
	end.chr <- end[chr.index]
	reads.chr <- reads[chr.index]
	
	if(sum(chr.index) > 1){
	chr.mat <- -1
	chr.temp <- 0
	i <- 1
	while(i <= length(start.chr)){
		# initialize
		if(i == 1){
			chr.temp <- c(start.chr[i],end.chr[i],reads.chr[i])
			i <- i + 1
		} else if(chr.temp[2] >= start.chr[i] - gap.size - 1){
			chr.temp[2] <- end.chr[i]
			chr.temp[3] <- chr.temp[3] + reads.chr[i]			
			i <- i + 1
		} else {
			if(chr.mat[1] == -1){				
				chr.mat <- chr.temp				
			} else {
				chr.mat <- rbind(chr.mat,chr.temp)
			}
			chr.temp <- cbind(start.chr[i],end.chr[i],reads.chr[i])
			i <- i + 1
		}		
	}
	# last peak
	if(!is.matrix(chr.mat)){
		chr.mat <- chr.temp				
	} else {
		chr.mat <- rbind(chr.mat,chr.temp)
	}
	
	chr.col <- infile[chr.index,1]
	#number of peaks
	if(is.matrix(chr.mat)){
		d <- dim(chr.mat)
		chr.col <- chr.col[1:d[1]]
	} else {
		chr.col <- chr
	}
	
	if(!is.data.frame(chr.mat.all)){
		if(is.matrix(chr.mat)){
			chr.mat.all <- cbind(as.character(chr.col),as.data.frame(chr.mat,stringsAsFactors=FALSE))		
		} else {
			chr.mat.all <- c(as.character(chr.col),chr.mat)
		}		
	} else{
		if(is.matrix(chr.mat)){
			chr.mat <- cbind(as.character(chr.col),as.data.frame(chr.mat,stringsAsFactors=FALSE))
			chr.mat.all <- rbind(chr.mat.all,chr.mat)
		}
		else{
			chr.mat <- cbind(chr.col,as.data.frame(t(chr.mat)))
			colnames(chr.mat) <- colnames(chr.mat.all)
			chr.mat.all <- rbind(chr.mat.all,chr.mat)
		}
	}
	} else {
		chr.mat <- infile[chr.index,1:4]
		if(!is.data.frame(chr.mat.all)){			
			chr.mat.all <- chr.mat
		} else{
			colnames(chr.mat) <- colnames(chr.mat.all)
			chr.mat.all <- rbind(chr.mat.all,chr.mat)
		}
	}	
}

if(!is.na(outfile)){
	write.table(chr.mat.all,file=outfile,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
} else {
	return(chr.mat.all)
}
progress <- paste("Total time elapsed = ",as.character(proc.time()[3] - starttime)," seconds\n",sep="")
cat(progress)
}