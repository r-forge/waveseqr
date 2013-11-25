preprocess <- function(infile,outdir,exptname,winsize=200,fragsize=150,redundancy=TRUE){
#
# preprocess(infile, outdir, exptname, winsize, fragsize, 
#               redundancy)
#
#   Preprocesses an input bed file to generate summary bedgraph 
#   and padded bedgraph files for use with WaveSeq. Also generates
#   chromosome information file.
#
#   infile      bed file to be preprocessed (required)
#   outdir      output directory (require)
#   exptname    unique string to identify name of experiment (required)
#   winsize     window size for calculating summary counts (default=200)
#   fragsize    Average fragment size from the sequencing experiment 
#               (default = 150)
#   redundancy  flag to remove redundancy before calculating summary counts
#               (default = 1)
#

#   Apratim Mitra 2011-2012

	# check if required parameters have been input
    if(missing(infile)) stop("Input file not specified!")
	else if(missing(outdir)) stop ("Output directory not specified!")
	else if(missing(exptname)) stop("Experiment name not specified!")
	
    # check folder hierarchy
    checkFolderHierarchy(outdir)
    
    # check that input file is a BED file
    if(regexpr('\\.bed',infile)[1] == -1) stop("Input file must be of BED format!")
    
    # start time counter
    starttime <- proc.time()[3]
    
    # input file name
    bedfile <- sub("\\s+","",infile)
    filename <- sub(".+/","",bedfile)  
    
    # print name of file being analyzed
    progress <- paste("Preparing ",filename," for analysis with WaveSeq\n",sep="")
    cat(progress)    

    # remove redundancy if required
    if(redundancy){
        # get name of nodup file
        filename <- sub(".bed","-nodup.bed",filename,fixed=TRUE)
        nodupfile <- file.path(outdir,"bed",filename)

        # print progress
        progress <- "\tRemoving redundancy ... "
		cat(progress)		

		cmd <- paste("perl ",file.path(get("WS.perlpath", envir=topenv(), inherits=FALSE),"rem_dup_reads.pl"),
					" -i ",bedfile," > ",nodupfile,sep="")
        system(cmd)

        elapsedtime <- proc.time()[3] - starttime        
		
        progress <- paste(as.character(elapsedtime)," seconds\n",sep="")
		cat(progress)		
        bedfile <- nodupfile
    }
    
    # convert file from BED to bedGraph format

    # get name of graph file
    filename <- sub(".bed",".graph",filename,fixed=TRUE)
    graphfile <- file.path(outdir,"graph",filename)

    progress <- "\tConverting from BED to bedGraph format ... "
	cat(progress)	

	cmd <- paste("perl ",file.path(get("WS.perlpath", envir=topenv(), inherits=FALSE),"bed2graph.pl")," -i ",
					bedfile," -w ",as.character(winsize)," -f ",
					as.character(fragsize)," > ",graphfile,sep="")
    system(cmd)

    elapsedtime <- proc.time()[3] - starttime

    progress <- paste(as.character(elapsedtime)," seconds\n",sep="")
	cat(progress)
    
    # Prepare graph files for use with WaveSeq
    # get name of padded graph file
    filename <- sub(".graph","-pad.graph",filename,fixed=TRUE)
    padfile <- file.path(outdir,"padgraph",filename)

    progress <- "\tPadding bedGraph files ... "
    cat(progress)
	
	cmd <- paste("perl ",file.path(get("WS.perlpath", envir=topenv(), inherits=FALSE),"pad_graph.pl")," -i ",
					graphfile," -w ",as.character(winsize)," > ",padfile,sep="")
    system(cmd)

	elapsedtime <- proc.time()[3] - starttime    

    progress <- paste(as.character(elapsedtime)," seconds\n",sep="")
	cat(progress)
    
    # get chromosome information
    chrfile <- paste("chrlist_",exptname,".txt",sep="")
	chrfile <- file.path(outdir,"padgraph",chrfile)
    progress <- "\tGetting chromosome information ... "
    cat(progress)

	cmd <- paste("perl ",file.path(get("WS.perlpath", envir=topenv(), inherits=FALSE),"get_chr_info.pl")," -i ",
					padfile," -o ", chrfile,sep="")
    system(cmd)

    elapsedtime <- proc.time()[3] - starttime    

    progress <- paste(as.character(elapsedtime)," seconds\n",sep="")
	cat(progress)
    
    # final progress message
	totaltime <- proc.time()[3] - starttime
    progress <- paste("Total time elapsed = ",as.character(totaltime)," seconds\n",sep="")
    cat(progress)
}