getThresholdDistribution <- function(infile, outdir, exptname, mother="morlet", N=5000, maxscale=12, 
							p.min=0.001, p.max=0.3){
# 
# get_thresdist(infile,outdir,exptname,N,scale,p.min,p.max,mother) 
#
# 	Obtains wavelet coefficient distributions for the specified
# 	input file using the continuous wavelet transform (CWT). 
# 	Monte Carlo sampling is performed with N samples of size 2^scale 
# 	and coefficients are obtained at a range of quantiles defined by 
# 	the range [1-p.min,1-p.max]. Sampling is performed chromosome-by-chromosome.
# 	
#	infile      input padded graph file (required)
#	outdir      output directory (required)
#   exptname    unique string to identify name of experiment (required)
#	N           number of samples (default=5000)
#	maxscale	maximum scale of wavelet decomposition. Also defines size of sample 
#               (2^sample) (default=12)
#	p.min        defines minimum quantile (1-p.min) to be output for wavelet coefficient 
#               distribution (default=0.3)
#	p.max        defines maximum quantile (1-p.max) as above (default=0.001)
#	mother      wavelet mother function used for the CWT. Tested choices are
#               "morlet","haar","gaussian1", "gaussian2" (default="morlet")
#

#   Apratim Mitra 2011-2012	

	# check if required parameters have been input
    if(missing(infile)) stop("Input file not specified!")
	else if(missing(outdir)) stop ("Output directory not specified!")
	else if(missing(exptname)) stop("Experiment name not specified!")

	# check folder hierarchy
    checkFolderHierarchy(outdir)
	
	# check that input file is a graph file
    if(regexpr("\\.graph",infile)[1] == -1) stop("Input file must be of bedGraph format!")
    
    # start time counter
    starttime <- proc.time()[3]
	    
    # sample length
    lensample <- 2^maxscale
    midpt <- lensample/2
    
    # p-value range preallocation
    prange <- seq((1-p.max),(1-p.min),p.min)
    pindex <- prange*N	
    
	filename <- sub(".+/","",infile)
    progress <- paste("Analyzing ",filename,"\n",sep="")
    cat(progress)    
    
    # file with chromosome information
    chrfile <- file.path(outdir,"padgraph",paste("chrlist_",exptname,".txt",sep=""))
    
    progress <- "\tGetting chromosome information ... "
    cat(progress) 
    
    # check to see if it exists already
    if(!file.exists(chrfile)){
		cmd <- paste("perl ",file.path(WS.perlpath,"get_chr_info.pl")," -i ",infile," -o ",chrfile,sep="")
		#cat("\n",cmd,"\n")
		system(cmd)		
    
		elapsedtime <- proc.time()[3] - starttime
        progress <- paste(as.character(elapsedtime)," seconds\n",sep="")
        cat(progress)
	} else {
        progress <- "skipped!\n"
        cat(progress)
    }    
	
    # read chrfile
    chrinfo <- read.table(chrfile, sep="\t", header=FALSE, stringsAsFactors=FALSE)
    
    # get chromosome info
    chrvec <- chrinfo[,1]
    numvec <- chrinfo[,2]
    chrnum <- length(numvec)
    totalreads <- numvec[chrnum]
    chrnum <- chrnum - 1
    
    #print summary
    progress <- paste("\t\tNumber of chromosomes = ",as.character(chrnum),"\n",sep="")
    cat(progress)

    progress <- paste("\t\tTotal reads in input file = ",as.character(totalreads),"\n",sep="")
    cat(progress)

    # open reads file
    fid <- file(infile,"r")
    
    #pre-allocation
    sample.reads <- mat.or.vec(maxscale,N)
    
    #conduct simulations
    skip <- 0
    currenttime <- proc.time()[3]

	for(k in 1:chrnum){
		progress <- paste("\t",chrvec[k]," ... ",sep="")
		cat(progress)

		#number of lines to be read
		numlines <- numvec[k]

		#read data from chromosome k
		chrdata <- scan(fid,what=list("","","",""),nlines=numlines,sep="\t",quiet=TRUE)
		reads <- as.numeric(chrdata[[4]])
		rm(chrdata)

		#reads normalized to per million mapped reads
		reads <- reads*1e6/totalreads

		# name of output file
		outfile <- file.path(outdir,"thresdist",paste(chrvec[k],"_",mother,"_",exptname,".txt",sep=""))

		currenttime <- proc.time()[3]
		N.chr <- floor(numlines*N/max(numvec[1:chrnum]))
		for(j in 1:N.chr){
			#print(N.chr)
			#if((j %% 50) == 0){
			#	test.time <- proc.time()[3] - currenttime
			#	cat(as.character(test.time),"\n")
			#}
			
			#obtain random sample of length=lensample from data y
			if(lensample < floor(numlines/2)){
				y1 <- sample(reads,lensample)
			} else {
				# skip chromosome if sample length < 0.5*numlines
				skip <- 1
				break
			}
			#sample.time <- proc.time()[3] - test.time
			#test.time <- proc.time()[3]
			
			#cwt with mother wavelet
			y1.cwt <- t(wavCWT(y1,scale.range=c(1,maxscale),wavelet=mother,n.scale=50))					
			
			#take Real part
			y1.cwt <- Re(as.matrix(y1.cwt))
			#print(dim(y1.cwt))
			#stop()
			
			# take slice from spectrum			
			sample.reads[,j] <- y1.cwt[,midpt]^2
		}

		# get distribution
		if(skip == 0){
			# output progress
			elapsedtime <- proc.time()[3] - currenttime
			currenttime <- proc.time()[3]
			progress <- paste(as.character(elapsedtime)," seconds\n",sep="")
			cat(progress)
			
			#sample.reads <- apply(sample.reads,1:2, function(x) x^2)
			sorted <- t(apply(sample.reads,1,function(x) x[order(x)]))
			thresdist <- sorted[,pindex]
			
			progress <- paste("\t\tPrinting to ",outfile,"\n",sep="")
			cat(progress)
			
			write.table(thresdist, outfile, sep="\t", row.names=FALSE,
				col.names=FALSE, quote=FALSE)			
		} else {
			progress <- "Skipped\n"
			cat(progress)
			skip <- 0
		}		
	}
    progress <- paste("Total time elapsed = ",as.character(proc.time()[3] - starttime)," seconds\n",sep="")
    cat(progress)
	close(fid)
}