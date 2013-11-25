\name{WaveSeq-package}
\alias{WaveSeq-package}
\alias{WaveSeq}
\docType{package}
\title{
WaveSeq: A Wavelet-based Algorithm for Detecting Histone Modification Enrichments
}
\description{
WaveSeq is a highly sensitive peak-caller for detecting punctate and diffuse enrichment regions 
from ChIP-Seq data. WaveSeq utilizes the continuous wavelet transform framework. The algorithm
involves the following steps: \cr
(1) Pre-processing input data.\cr
(2) Monte Carlo sampling to estimate wavelet coefficient thresholds.\cr
(3) Use computed thresholds to call putative enrichment regions.\cr
(4) Merging peaks within a user-specified distance \cr
(5) Estimating FDR in the presence or absence of control data \cr
(6) Correcting p-values for multiple testing.
}
\details{
\tabular{ll}{
Package: \tab WaveSeqR\cr
Type: \tab Package\cr
Version: \tab 1.0.2\cr
Date: \tab 2013-04-29\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
%% An overview of how to use the package, including the most important functions %%
}
\author{
Apratim Mitra, Jiuzhou Song

Maintainer: Apratim Mitra <amitra83@umd.edu>, Jiuzhou Song <songj88@umd.edu>
}
\references{
%% Literature or other references for background information %%
}
\keyword{ package }
\keyword{ robust }
\keyword{ nonparametric }
\keyword{ distribution }
\keyword{ iteration }
\keyword{ htest }

\examples{
# load package
library(WaveSeqR)

# get paths to data files
WS.path <- system.file(package="WaveSeqR")
chip.file <- file.path(WS.path,"extdata","gabp_valouev2008_chr22-15m-30m-pad.graph")
control.file <- file.path(WS.path,"extdata","rxnoip_valouev2008_chr22-15m-30m-pad.graph")
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
