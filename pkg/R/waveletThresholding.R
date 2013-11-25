waveletThresholding <- function(infile, outdir, exptname, mother="morlet", winsize=200,
               minsig=2, minsigscale=3, maxsigscale=-1, p.thres=0.2,p.min=0.001, 
               p.max=0.3){
#
# waveseq_thres(infile, outdir, exptname, 'mother', 'morl', 'minsigscale', 3,
#               'minsig', 2, 'pthres', 0.2, 'winsize', 200, 'p.min', 0.001, 
#               'p.max', 0.3) 
#   Apply precomputed wavelet coefficient threshold to continuous
#   wavelet transform (CWT) of specified input data. Thresholding 
#   performed chromosome-by-chromosome. Only scales > minscale considered
#   for peak detection. Number of significant scales has to be > minsig 
#   for the window to be considered significant.
#
#   infile      input padded graph file (required)
#   outdir      output directory (required)
#   exptname    unique string (matching thresdist files) for a
#               specific experiment (required)
#   mother      wavelet mother function (default='morl')
#   minscale    minimum scale considered for peak detection (default=3)
#   minsig      minimum number of significant scales for a window to be
#               considered significant (default=2)
#   pthres      p-value for calling significant windows (default=0.2)
#   winsize     window size for padded graph files (default=200)
#   p.min        defines minimum quantile (1-p.min) calculated for the wavelet 
#               coefficient distribution (default=0.001)
#	p.max        defines maximum quantile (1-p.max) as above (default=0.3)
#

#   Apratim Mitra 2011-2012	

	#check if required parameters have been input
    if(missing(infile)) stop("Input file not specified!")
	else if(missing(outdir)) stop ("Output directory not specified!")
	else if(missing(exptname)) stop("Experiment name not specified!")
    
	# start time counter
    starttime <- proc.time()[3]     
    
    # check that input file is a graph file
    if(regexpr("\\.graph",infile)[1] == -1) stop("Input file must be of bedGraph format!")
	    
    # check folder hierarchy
    checkFolderHierarchy(outdir)
    
	# write progress into log file
    filename <- sub(".+/","",infile)    
	    
    # print name of file being analyzed
	#progress <- paste("Analyzing ",filename,"\n",sep="")
    #cat(progress)
    
    # Get total reads, chromosome names and counts
	chrfile <- file.path(outdir,"padgraph",paste("chrlist_",exptname,".txt",sep=""))
    
    progress <- "\tGetting chromosome information ... "
    cat(progress)
    
    # check to see if it exists already
    if(!file.exists(chrfile)){
		cmd <- paste("perl ",file.path(WS.perlpath,"get_chr_info.pl"),
					" -i ",infile," > ",chrfile,sep="")
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
	    
    #read threshold distribution file    	
    cat("\tReading wavelet coefficient thresholds ... ")
    
    inmean <- file.path(outdir,"thresdist",paste("thres_avg_",mother,"_",exptname,".txt",sep=""))
    instd <- file.path(outdir,"thresdist",paste("thres_std_",mother,"_",exptname,".txt",sep=""))
    
    cat("done\n")
    
    # get column corresponding to pthres
	pindex <- seq((1-p.max),(1-p.min),by=p.min)
	pthresind <- length(pindex) - (p.max - p.thres)/p.min + 1
    thresmean <- read.table(inmean,sep="\t",header=FALSE)
	thresstd <- read.table(instd,sep="\t",header=FALSE)
	thresmat <- thresmean + thresstd
	    
    # vector containing coefficient thresholds
    pvec <- thresmat[,pthresind]	

    #maximum scale of wavelet decomposition
	maxscale <- dim(thresmean)[1]
    if(maxsigscale < 0){
        maxsigscale <- maxscale
    }

    # open data to be analyzed
    fin <- file(infile,"r")

    #open output file
    #outfile <- file.path(outdir,"peaks",paste(exptname,"_",mother,"_W",as.character(winsize),"_p",
	#				as.character(p.thres),"_nogap_min",as.character(minsigscale),"_max",
	#				as.character(maxsigscale),".txt",sep=""))
	outfile <- file.path(outdir,"peaks",paste(exptname,"_",mother,"_W",as.character(winsize),"_p",
					as.character(p.thres),"_nogap.txt",sep=""))
    fout <- file(outfile,"w+")

	currenttime <- proc.time()[3]
	
    # peak detection chromosome-by-chromosome
    for(j in 1:chrnum){
        progress <- paste("\t",chrvec[j]," ... ",sep="")
        cat(progress)
            
        #number of lines to be read
        numlines <- numvec[j]

        #read data from chromosome j
		chrdata <- scan(fin,what=list("","","",""),nlines=numlines,sep="\t",quiet=TRUE)
        reads.raw <- as.numeric(chrdata[[4]])
		rm(chrdata)
		cat("\n\t\tRead data ... ",as.character(proc.time()[3] - currenttime)," seconds\n")
		currenttime <- proc.time()[3]
		
        #reads normalized to per million mapped reads
        reads <- reads.raw*1e6/totalreads
		
		padlen <- 2^(floor(log2(length(reads))) + 1) - length(reads)
		reads <- c(reads,mat.or.vec(padlen,1))
		#cat(as.character(padlen),"\t",as.character(length(reads)),"\n")

        #perform cwt
		coefs <- t(wavCWT(reads,scale.range=c(1,maxscale),wavelet=mother,n.scale=50))
		#rm(reads)
		cat("\t\tCalculated CWT ... ",as.character(proc.time()[3] - currenttime)," seconds\n")
		currenttime <- proc.time()[3]
		
		coefs <- Re(as.matrix(coefs))
		#print(coefs[,1]^2)
		
        #power spectrum
        pwrspec <- t(apply(coefs,1,function(x) x^2))
		#print(pwrspec[1,])
		#stop()
		
        # preallocate for speed
        pwrspec1 <- mat.or.vec(dim(pwrspec)[1],dim(pwrspec)[2])
        rm(coefs)
		
        for(i in minsigscale:maxscale){
            zeroind <- pwrspec[i,] > pvec[i]
            pwrspec1[i,zeroind] <- pwrspec[i,zeroind]
        }

        #converting thresholded power spectrum to binary matrix
        pwrspec2 <- matrix(as.integer( pwrspec1 > 0 ),nrow=dim(pwrspec1)[1],ncol=dim(pwrspec1)[2])

		cat("\t\tCalculated power spectrum ... ",as.character(proc.time()[3] - currenttime)," seconds\n")
		currenttime <- proc.time()[3]
		
        #only consider scales >= minsigscale || scales <= maxsigscale
        pwrspec2[1:minsigscale,] <- 0
        pwrspec2[maxsigscale:maxscale,] <- 0
        pwrspec2 <- apply(pwrspec2,2,sum)

        #important tuning parameter - minimum of scales that have to be
        #significant to be marked as a potential peak
		#pwrspec2[ pwrspec2 <= minsig ] <- 0

        #convert binary vector to locations
        matsize <- length(pwrspec2)
        loc <- (1:matsize)*winsize
        loc[ pwrspec2 <= minsig ] <- 0

        #aggregate contiguous regions for putative peaks
		cat("\t\tAggregating contiguous regions for putative peaks ...")
        #locfinal = [];    
        peakloc <- c(0,0)    
        peakreads <- 0
        for(i in 1:numlines){
            # a new peak
            if(peakloc[1] == 0 && loc[i] > 0 && reads.raw[i] > 0){
                peakloc <- c(loc[i]-winsize, loc[i]-1)
                peakreads <- reads.raw[i]
            } else if(peakloc[1] > 0 && loc[i] > 0){
                # extend existing peak only if window contains reads
                if( reads.raw[i] > 0 ){
                    peakloc[2] <- loc[i] - 1
                    peakreads <- peakreads + reads.raw[i]
                #output peak and start afresh
                } else {
                    if(peakreads > 0){
						cat(paste(chrvec[j],as.character(format(peakloc[1],scientific=FALSE)),
							as.character(format(peakloc[2],scientific=FALSE)),paste(as.character(peakreads),"\n",sep=""),sep="\t"),
							file=fout,sep="",append=TRUE)
                        #locfinal = [locfinal; peakloc peakreads];
                    }
                    peakloc <- c(0,0)
					peakreads <- 0
                }
            } else if(peakloc[1] > 0 && loc[i] == 0){
                if( peakreads > 0 ){
					cat(paste(chrvec[j],as.character(format(peakloc[1],scientific=FALSE)),
						as.character(format(peakloc[2],scientific=FALSE)),paste(as.character(peakreads),"\n",sep=""),sep="\t"),
						file=fout,sep="",append=TRUE)
                    #locfinal = [locfinal; peakloc peakreads];
                }
                peakloc <- c(0,0)
				peakreads <- 0
            }
        }
		
		elapsedtime <- proc.time()[3] - currenttime
		currenttime <- proc.time()[3]
		progress <- paste(as.character(elapsedtime)," seconds\n",sep="")
		cat(progress)
    }
	progress <- paste("Total time elapsed = ",as.character(proc.time()[3] - starttime)," seconds\n",sep="")
    cat(progress)
	close(fout)
	close(fin)
}