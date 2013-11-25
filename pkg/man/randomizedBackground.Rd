\name{randomizedBackground}
\alias{randomizedBackground}
\title{
	Randomized Background Estimation
}
\description{
	This function implements a randomized algorithm to generate the empirical cumulative 
	density function (ecdf) for the probability of a peak containing a certain number of reads. 
	It works by sampling a large number of peaks from the peak length distribution, randomizing 
	their locations on the genome and counting the number of reads within this 'random peak'. 
	The read counts within these random peaks are then used to generate the ecdf and assign 
	empirical p-values to each peak.
}
\usage{
randomizedBackground(peaks, chip, samplesize = 1e+06, winsize = 200)
}
\arguments{
  \item{peaks}{
	File containing peak locations. Must be in bedGraph format.
  }
  \item{chip}{
	Padded bedGraph file for the ChIP data being analyzed.
  }
  \item{samplesize}{
	Number of peaks to sample for building ecdf with the randomized algorithm.
  }
  \item{winsize}{
	Window size used for obtaining peaks(default = 200 bp).
  }
}
\details{
	For a detailed description of the bedGraph file format see 
	http://genome.ucsc.edu/goldenPath/help/bedgraph.html.
}
\value{
	Returns a matrix containing peak locations and their corresponding p-values 
	sorted by chromosome. Has the following columns: 
	\item{chromosome}{Chromosome name}
	\item{start}{Start position of peak (1-based indexing)}
	\item{end}{End position of peak}
	\item{reads}{Number of reads in the peak}
	\item{p-value}{P-value for the peak obtained from the estimated ecdf}
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

\section{Warning}{
	For small window sizes for the padded bedGraph files and consequently 
	large file sizes, reading the data may take a long time. Run time may 
	also be increased by larger values of \code{samplesize}.
}
\seealso{
	\code{\link{sample}}, \code{\link{ecdf}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ robust }
\keyword{ nonparametric }
\keyword{ distribution }
\keyword{ iteration }
