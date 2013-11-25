\name{compareWithControl}
\alias{compareWithControl}
\title{
	Compare ChIP to Control
}
\description{
	Part of the \code{\link{controlFDR}} function. This function reads in a list 
	of putative peaks, fetches reads in those regions from the control data and compares the 
	two using a two-sided exact Binomial test. Both sets of reads are normalized to per million 
	mapped reads before performing the test.
}
\usage{
compareWithControl(peaks, chip, control, winsize = 200, minreads = 8, binom.sided = "two.sided")
}
\arguments{
  \item{peaks}{
	File containing peak locations. Must be in bedGraph format.
  }
  \item{chip}{
	Padded bedGraph file for the ChIP data being analyzed.
  }
  \item{control}{
	Padded bedGraph file for the control data being analyzed.
  }
  \item{winsize}{
	Window size used for obtaining peaks(default = 200 bp).
  }
  \item{minreads}{
	Minimum reads in a peak (default = 8). This filtering is necessary for a more accurate 
	estimation of false discovery rate.
  }
  \item{binom.sided}{
	One-sided or two-sided binomial test, if control is not NA (default = "two.sided")
  }
}
\details{
	For a detailed description of the bedGraph file format see 
	http://genome.ucsc.edu/goldenPath/help/bedgraph.html.
}
\value{
	Returns a matrix containing peak locations, fold changes and their corresponding 
	p-values sorted by chromosome. Has the following columns:
		\item{chromosome}{Chromosome name}
		\item{start}{Start position of peak (1-based indexing)}
		\item{end}{End position of peak}
		\item{ChIP.reads}{Number of ChIP reads in the peak (normalized to per million mapped reads)}
		\item{Control.reads}{Number of control reads in the peak (normalized to per million mapped reads)}
		\item{Fold.change}{Fold change calculated as chip/control}
		\item{p-value}{P-value for the peak obtained from two-sided exact binomial test \code{\link{binom.test}}}
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
	large file sizes, reading the data may take a long time.
}

\seealso{
	\code{\link{binom.test}}, \code{\link{getPeakReads}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ distribution }
