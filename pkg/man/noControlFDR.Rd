\name{noControlFDR}
\alias{noControlFDR}
\title{
	FDR Estimation without Control
}
\description{
	Estimates FDR for peaks in a ChIP sample using a randomized algorithm. 
	\code{noControlFDR} consists of three steps: \cr
	(1) Randomized background estimation using \code{\link{randomizedBackground}()}. \cr
	(2) Adjust p-values for multiple testing using \code{\link{p.adjust}}. \cr
	(3) Print peaks and p-values to file.\cr
}
\usage{
noControlFDR(peaks, chip, outfile, winsize = 200, 
			minreads = 8, samplesize = 1e+06, adj.method = "fdr")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peaks}{
	File containing peak locations. Must be in bedGraph format.
  }
  \item{chip}{
	Padded bedGraph file for the ChIP data being analyzed.
  }  
  \item{outfile}{
	Output file.
  }
  \item{winsize}{
	Window size for obtaining peaks (default = 200).
  }
  \item{minreads}{
	Minimum reads within a peak (default = 8). This filtering is necessary for a more
	accurate estimation of FDR.
  }
  \item{samplesize}{
	Number of peaks to sample for building ecdf with the randomized algorithm.
  }
  \item{adj.method}{
	Method for multiple testing adjustment. Can use any methods accepted by \code{\link{p.adjust}}
	(default = "fdr").
  }
}
\value{
	This function writes peaks to a file after FDR estimation in the absence of control. This
	file has the following fields:
	\item{chromosome}{Chromosome name}
	\item{start}{Start position of peak (1-based indexing)}
	\item{end}{End position of peak}
	\item{ChIP.reads}{
		Number of ChIP reads in the peak }	
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
	\code{\link{sample}}, \code{\link{ecdf}}, \code{\link{p.adjust}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ robust }
\keyword{ nonparametric }
\keyword{ distribution }
\keyword{ iteration }
