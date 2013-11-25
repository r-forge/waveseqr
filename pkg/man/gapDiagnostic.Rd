\name{gapDiagnostic}
\alias{gapDiagnostic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to Select Suitable Gap Size
}
\description{
Calculates and plots fraction of total reads in \code{topN} peaks for a specified range 
of gap sizes. 
}
\usage{
gapDiagnostic(peaks, gap.max = 15, topN = 10000, winsize = 200, plot.flag = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peaks}{
	\code{\link{data.frame}} or file containing list of peaks output by WaveSeqR.
}
  \item{gap.max}{
	Maximum gap size for the gap diagnostic calculation (default = 15).
}
  \item{topN}{
	Number of peaks to perform above calculation (default=1000). Should be adjusted to equal
	the maximum number of expected peaks. 
}
  \item{winsize}{
	Window size used for peak-calling and bedGraph files.
}
  \item{plot.flag}{
	If \code{TRUE}, plots the gap diagnostic curve.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	Returns a vector of length \code{gap.max + 1} with fraction of total reads in \code{topN} peaks
	for gap sizes ranging from 0 to \code{gap.max}.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
	\code{\link{gapPeaks}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

