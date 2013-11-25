\name{waveseq}
\alias{waveseq}
\title{
	Run WaveSeq Analysis Pipeline
}
\description{	
	This is the main function of the WaveSeq algorithm modules implemented in R. There
	are six processing steps:\cr
	(1) Pre-process the input BED file.\cr
	(2) Perform Monte Carlo sampling to estimate wavelet coefficient thresholds using the 
		continuous wavelet transform (CWT).\cr
	(3) Apply thresholds to CWT of input data to obtain putative peaks.\cr
	(4) Concatenate peaks within a user-specified distance (optional)\cr
	(5) If control data is specified, compare reads within peaks to control with an exact binomial test.\cr
	(5) If control data is absent, use randomized algorithm to estimate background distribution and
		assign p-values to peaks.\cr
	(6) Correct p-values for multiple-testing. \cr	
}
\usage{
waveseq(chip, control=NA, outdir, exptname, preprocess = TRUE, redundancy = TRUE, thresdist = TRUE, 
	peak.calling = TRUE, mother = "morlet", winsize = 200, fragsize = 150, gap = 0, minreads = 8,
	samplesize = 1e6, adj.method = "fdr", p.thres = 0.2, maxscale = 12, N = 5000, minsig = 2, 
	minsigscale = 3, maxsigscale = -1, p.min = 0.001, p.max = 0.3, binom.sided = "two.sided")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chip}{
	Input file in BED format. If \code{preprocess = FALSE}, \code{chip} is assumed to be in 
	bedGraph format (required).
  }
  \item{control}{
	Padded bedGraph file for the control data being analyzed (default = NA).
  }  
  \item{outdir}{
	Output directory (required).
  }
  \item{exptname}{
	Unique string to identify name of experiment (required).
  }
  \item{preprocess}{
	Logical value to indicate whether input file should be preprocessed or not. If missing or \code{FALSE},
	input file is assumed to be in bedGraph format (default=TRUE).
  }
  \item{redundancy}{
	Logical value to indicate whether redundant reads should be removed 
	before calculating summary counts (default=TRUE).
  }
  \item{thresdist}{
	Logical value to indicate whether threshold estimation step is to be performed or not. If missing or
	\code{FALSE}, thresholds will be assumed to exist in \code{thresdist} folder of \code{outdir}
	(default=TRUE).
  }
  \item{peak.calling}{
	Logical value to indicate whether peak-calling step is to be performed or not (default=TRUE).
  }
  \item{mother}{
	Wavelet mother function used for the CWT. Tested choices are
	"morlet","haar","gaussian1", "gaussian2" (default="morlet").
  }
  \item{winsize}{
	Window size for padded graph files (default=200).
  }
  \item{fragsize}{
	Library size used for sequencing (default=150).
  }
  \item{gap}{
	Maximum number of non-significant windows separating peaks that are to be concatenated (default = 0).
  }
  \item{minreads}{
	Minimum reads within a peak (default = 8). This filtering is necessary for a more
	accurate estimation of FDR.
  }
  \item{samplesize}{
	Number of peaks to sample for building ecdf with the randomized algorithm (default = 1e6).
  }
  \item{adj.method}{
	Method for multiple testing adjustment. Can use any methods accepted by \code{\link{p.adjust}} 
	(default = "fdr").
  }
  \item{p.thres}{
	Threshold p-value for calling significant windows (default=0.2)
  }
  \item{maxscale}{
	Maximum scale of wavelet decomposition. Also defines size of sample (2^maxscale) (default=12).
  }
  \item{N}{
	Number of samples (default=5000)
  }
  \item{minsig}{
	Minimum number of significant scales for a window to be considered significant (default=2).
  }
  \item{minsigscale}{
	Minimum significant scale considered for peak detection (default=3).
  }
  \item{maxsigscale}{
	Minimum significant scale considered for peak detection (default=-1).
  }  
  \item{p.min}{
	Minimum quantile (1-p.min) to be output for wavelet coefficient distribution (default=0.3).
  }
  \item{p.max}{
	Maximum quantile (1-p.max) as above (default=0.001).
  }
  \item{binom.sided}{
	One-sided or two-sided binomial test, if control is not NA (default = "two.sided")
  }
}
\details{
	This is the main function to apply the WaveSeq algorithm to Next-generation sequencing
	data e.g. ChIP-Seq.
}
\value{
	This function concatenates peaks within a specified distance, estimates FDR 
	in the presence or absence of control data and writes peaks to file. This
	file has the following fields:
	\item{chromosome}{Chromosome name}
	\item{start}{Start position of peak (1-based indexing)}
	\item{end}{End position of peak}
	\item{ChIP.reads}{
		Number of ChIP reads in the peak (normalized to per 
		million mapped reads if \code{control} != NA)}	
	\item{Control.reads}{
		Number of control reads in the peak (normalized to per
		million mapped reads); Output only if \code{control} != NA}
	\item{Fold.Change}{
		Fold-change calculated as \code{ChIP.reads}/\code{Control.reads};
		Output only if \code{control} != NA}
	\item{p.value}{P-value for the peak}
	\item{adj.p.value}{P-values adjusted for multiple comparisons using \code{adj.method}}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
	Apratim Mitra
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{sample}}, \code{\link{ecdf}}, \code{\link{p.adjust}}, \code{\link{binom.test}}
}
\examples{
\dontrun{
# load package
library(WaveSeqR)

# get paths to data files
WS.path <- system.file(package="WaveSeqR")
chip.file <- file.path(WS.path,"inst","extdata","gabp_valouev2008_chr22-pad.graph")
control.file <- file.path(WS.path,"inst","extdata","rxnoip_valouev2008_chr22-pad.graph")
test.dir <- file.path(WS.path,"extdata")

# ChIP-Seq analysis with control
waveseq(chip=chip.file,
	outdir=test.dir,
	control=control.file,
	exptname="GABP_valouev_chr22_control",
	preprocess=FALSE,
	thresdist=FALSE,
	gap=0,
	winsize=200,
	minreads=8,			
	adj.method="fdr")
	
# ChIP-Seq analysis without control
waveseq(chip=chip.file,
	outdir=test.dir,
	exptname="GABP_valouev_chr22_nocontrol",
	preprocess=FALSE,
	thresdist=FALSE,
	gap=0,
	winsize=200,
	minreads=8,
	samplesize=1e6,
	adj.method="fdr")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ robust }
\keyword{ nonparametric }
\keyword{ distribution }
\keyword{ iteration }
\keyword{ htest }