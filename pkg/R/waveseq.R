waveseq <- function(chip, control=NA, outdir, exptname, preprocess=TRUE, redundancy=TRUE, thresdist=TRUE,
						peak.calling=TRUE, mother="morlet", winsize=200, fragsize=150, gap=0, minreads=8,
						samplesize=1e6, adj.method="fdr", p.thres=0.2, maxscale=12, N=5000, minsig=2, 
						minsigscale=3, maxsigscale = -1, p.min=0.001, p.max=0.3, binom.sided="two.sided"){
# 
#   Main function to run waveseq on specified data. Pre-process the input
#   BED file, run Monte Carlo sampling to estimate wavelet coefficient
#   thresholds, output putative peaks.
#
#   infile          Input BED file (required)
#   outdir          Ouput directory (required)
#   exptname        Unique string to identify name of experiment (required)
#   prep            Pre-process the input data (default = 1)
#   thresdist       Monte Carlo estimation of wavelet coefficient threshold
#                   (default = 1)   
#   waveseq         Peak detection with WaveSeq (default = 1)
#   winsize         Window size for summary bedGraph files (default=200)
#   fragsize        Average fragment size from sequencing experiment
#                   (default = 150)
#   redundancy      Remove redundancy (default = 1)
#   mother          Mother wavelet used for WaveSeq algorithm (default =
#                   'morl')
#   N               Number of samples in Monte Carlo step (default = 5000)
#   maxscale		Maximum scale for CWT computation (default = 12)
#   minscale        Only scales > minscale considered for peak detection 
#                   (default = 3)
#   minsig          Minimum number of significant scales for a window to be
#                   considered significant (default = 2)
#   pthres          P-value for calling significant windows (default=0.2)
#   pmin            Minimum quantile (1-pmin) calculated for the wavelet 
#                   coefficient distribution (default = 0.001)
#	pmax            Maximum quantile (1-pmax) as above (default = 0.3)
#

#   Apratim Mitra 2011-2012

	# check if required parameters have been input
    if(missing(chip)) stop("Input file not specified!")
	else if(missing(outdir)) stop ("Output directory not specified!")
	else if(missing(exptname)) stop("Experiment name not specified!")
	
	# check folder hierarchy
    checkFolderHierarchy(outdir) 
	
	logfile <- file.path(outdir,"log",paste(exptname,"_",mother,"_gap",as.character(gap),".log",sep=""))
    flog <- file(logfile,"w+")
	
	sink(file=flog, append=FALSE, split=TRUE)
	# print program version and parameter summary	
	cat("\nWaveSeq v 1.0\n\n")
	cat("Parameter list\n")
	cat("--------------\n")	
	cat("ChIP file :",chip,"\n")
	if(!is.na(control)){
		cat("Control file : ",control,"\n")	
	} else {
		cat("Control file : NA\n")	
	}	
	cat("Output directory : ",outdir,"\n")
	cat("Experiment name : ",exptname,"\n")
	cat("Mother wavelet : ",mother,"\n")
	cat("Window size : ",winsize,"\n")
	cat("Fragment size : ",fragsize,"\n")
	cat("Gaps : ",gap,"\n")	
	if(thresdist) cat("Max. iterations for wavelet threshold determination : ",N,"\n")
	cat("Max. scale for CWT computation : ",maxscale,"\n")
	if(maxsigscale == -1){
		cat("Scales to be considered for peak-calling : ",minsigscale,"-",maxscale,"\n")
		maxsigscale <- maxscale
	} else {
		cat("Scales to be considered for peak-calling : ",minsigscale,"-",maxsigscale,"\n")
	}
	cat("P-value for calling significant windows : ",p.thres,"\n")
	cat("Minimum reads in a peak : ",minreads,"\n")
	if(is.na(control))	cat("# of peaks sampled for FDR estimation : ",samplesize,"\n")	
	cat("Multiple-testing adjustment : ",adj.method,"\n\n")        
    
    # updating threshold calculation limits
    if(p.thres > p.max){
    	p.max <- p.thres + p.min
    } else if(p.thres < p.min){
    	p.min <- p.thres
    }
    chipfilename <- sub(".+/","",chip)	
    if(!is.na(control)) controlfilename <- sub(".+/","",control)
	
	# Preprocess the data   
    if(preprocess){        
        if(!is.na(control)){
	    progress <- paste("\nPreprocessing ",controlfilename," ...\n",sep="")
    	    cat(progress)		
            preprocess(control,outdir,exptname,winsize,fragsize,redundancy)
        }

        progress <- paste("\nPreprocessing ",chipfilename," ...\n",sep="")
        cat(progress)
        preprocess(chip,outdir,exptname,winsize,fragsize,redundancy)

    }    
	
	if(preprocess){
		# get name of padded graph file
		if(redundancy){
			chipfilename <- sub(".bed","-nodup-pad.graph",chipfilename,fixed=TRUE)
		} else {
			chipfilename <- sub(".bed","-pad.graph",chipfilename,fixed=TRUE)        
		}
		chippadfile <- file.path(outdir,"padgraph",chipfilename)
		
		if(!is.na(control)){
			if(redundancy){
				controlfilename <- sub(".bed","-nodup-pad.graph",controlfilename,fixed=TRUE)
			} else {
				controlfilename <- sub(".bed","-pad.graph",controlfilename,fixed=TRUE)        
			}
			controlpadfile <- file.path(outdir,"padgraph",controlfilename)		
		}
	} else {
		chippadfile <- chip
		controlpadfile <- control
	}
    
    # get wavelet coefficient distribution
    if(thresdist){
        progress <- "\nGetting wavelet coefficient distribution ...\n"
        cat(progress)
        getThresholdDistribution(chippadfile, outdir, exptname, mother=mother, N=N,
						maxscale=maxscale, p.min=p.min, p.max=p.max)
        
        # calculate mean and standard deviation of wavelet coefficients
        thresdir <- file.path(outdir,"thresdist")
        getConfidenceInterval(thresdir, exptname, mother=mother)
    }
        
    # perform waveseq peak calling
    if(peak.calling){
        progress <- "\nRunning WaveSeq peak detection algorithm ...\n"
        cat(progress)
        waveletThresholding(chippadfile, outdir, exptname, mother=mother, winsize=winsize,
						minsig=minsig, minsigscale=minsigscale, maxsigscale=maxsigscale,
						p.thres=p.thres, p.min=p.min, p.max=p.max)
    }
	#peaks <- file.path(outdir,"peaks",
	#				paste(exptname,"_",mother,"_W",as.character(winsize),"_p",
	#				as.character(p.thres),"_nogap_min",as.character(minsigscale),"_max",
	#				as.character(maxsigscale),".txt",sep=""))
	peaks <- file.path(outdir,"peaks",
					paste(exptname,"_",mother,"_W",as.character(winsize),"_p",
					as.character(p.thres),"_nogap.txt",sep=""))
	
	# gap peaks
	if(gap > 0){		
		if(regexpr("_nogap_",peaks,fixed=TRUE)[1] != -1){
			gapfile <- sub("_nogap_",paste("_gap",as.character(gap),"_",sep=""),peaks,fixed=TRUE)
		} else {	
			gapfile <- sub(".txt",paste("_gap",as.character(gap),".txt",sep=""),peaks,fixed=TRUE)
		}
		gapPeaks(peaks,gapfile,gap=gap,winsize=winsize)
	} else {
		gapfile <- peaks
	}	
	
	# estimate statistical significance of peaks
	if(is.na(control)){
		outfile <- sub(".txt","_nocontrol_fdr.txt",gapfile,fixed=TRUE)		
		noControlFDR(gapfile, chippadfile, outfile, winsize=winsize, minreads=minreads,
								samplesize=samplesize, adj.method=adj.method)
	} else {
		outfile <- sub(".txt","_control_fdr.txt",gapfile,fixed=TRUE)		
		controlFDR(gapfile, chippadfile, controlpadfile, outfile, winsize=winsize,
							minreads=minreads, adj.method=adj.method, binom.sided=binom.sided)
	}	
	cat("\nAll processes complete!\n\n")
	sink()
	closeAllConnections()
}