getConfidenceInterval <- function (indir, exptname, mother="morlet"){
#
# get_conf(indir,exptname,mother)
#
#   Obtains means and standard deviations from thresdist files
#   whose names match the specified string 'exptname'
#
#   indir       directory containing thresdist files (required)
#   exptname    unique string matching thresdist files for a
#               specific experiment (required)
#   mother      wavelet mother function used to calculate coefficient
#               distributions (default = 'morlet')
#

#   Apratim Mitra 2011-2012

	if(nargs() < 2){
		stop("get_conf() --> Not enough input arguments!",call.=FALSE)
	}
    # get list of thresdist files	
    files <- list.files(path=indir, pattern=paste(mother,exptname,sep="_"))
	files.ind <- grep("^chr",files)
	files <- files[files.ind]
    numfiles <- length(files)

    if(numfiles == 0) stop(paste("No files matching ",exptname," found in ",indir,sep=""))
     
	# read thresdist into list
	conf <- list()
    for(i in 1:numfiles){        
		infile <- file.path(indir,files[i])
        thresdist <- read.table(infile,sep="\t",header=FALSE)        
		conf[[i]] <- thresdist            
    }
	# convert list to 3D array
	conf.arr <- array(unlist(conf),dim=c(dim(thresdist),numfiles))
	
    # calculate mean and std deviation
    conf.avg <- apply(conf.arr, c(1,2), mean)
    conf.std <- apply(conf.arr, c(1,2), sd)
	conf.std[is.na(conf.std)] <- 0

    meanfile <- file.path(indir,paste("thres_avg_",mother,"_",exptname,".txt",sep=""))
	stdfile <- file.path(indir,paste("thres_std_",mother,"_",exptname,".txt",sep=""))    

	write.table(conf.avg, file=meanfile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
	write.table(conf.std, file=stdfile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)    
}