\name{waveletThresholding}
\alias{waveletThresholding}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Perform Wavelet CWT Thresholding
}
\description{
	Apply precomputed wavelet coefficient threshold to continuous
	wavelet transform (CWT) of specified input data. Thresholding 
	performed chromosome-by-chromosome. Only scales > minscale considered
	for peak detection. Number of significant scales has to be > minsig 
	for the window to be considered significant.
}
\usage{
waveletThresholding(infile, outdir, exptname, mother = "morlet", winsize = 200, minsig = 2, minsigscale = 3, 
	maxsigscale = -1, p.thres = 0.2, p.min = 0.001, p.max = 0.3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infile}{
	Input padded graph file (required)
}
  \item{outdir}{
	Output directory (required)
}
  \item{exptname}{
	Unique string to identify name of experiment (required)
}
  \item{mother}{
	Wavelet mother function used for the CWT. Tested choices are
	"morlet","haar","gaussian1", "gaussian2" (default="morlet")
}
  \item{winsize}{
	Window size for padded graph files (default=200)
}
  \item{minsig}{
	Minimum number of significant scales for a window to be
	considered significant (default=2)
}
  \item{minsigscale}{
	Minimum significant scale considered for peak detection (default=3)
}
  \item{maxsigscale}{
	Minimum significant scale considered for peak detection (default=-1)
}
  \item{p.thres}{
	Threshold p-value for calling significant windows (default=0.2)
}
  \item{p.min}{
	Minimum quantile (1-p.min) to be output for wavelet coefficient 
	distribution (default=0.3)
}
  \item{p.max}{
	Maximum quantile (1-p.max) as above (default=0.001)
}
}
\details{
	This function writes the results from the wavelet thresholding into a file
	within the \code{peaks} folder of \code{outdir}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
	Apratim Mitra
}
\note{
	Function messages are piped to a log file which can be found in the \code{log}
	folder of \code{outdir}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{wavCWT}}
}
\examples{

}
