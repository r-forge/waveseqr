\name{getThresholdDistribution}
\alias{getThresholdDistribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wavelet Threshold Estimation
}
\description{	
	Obtains wavelet coefficient distributions for the specified	input file using the continuous 
	wavelet transform (CWT). Monte Carlo sampling is performed with N samples of size 2^maxscale 
 	and coefficients are obtained at a range of quantiles defined by the range [1-p.min,1-p.max]. 
	Statified sampling is performed chromosome-by-chromosome weighted by chromosome length.
}
\usage{
getThresholdDistribution(infile, outdir, exptname, mother = "morlet", N = 5000, maxscale = 12, p.min = 0.001, p.max = 0.3)
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
  \item{N}{
	Number of samples (default=5000)
}
  \item{maxscale}{
	Maximum scale of wavelet decomposition. Also defines size of sample 
	(2^maxscale) (default=12)
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
	This function outputs wavelet coefficient thresholds in increments of p.max
	for each chromosome to the \code{thresdist} folder of \code{outdir}. 
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
\section{Warning}{
	For large values of \code{N} and small values of \code{winsize}
	this step can be time-consuming.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{wavCWT}}
}
\examples{

}
\keyword{distribution}