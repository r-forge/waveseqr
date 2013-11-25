\name{gapPeaks}
\alias{gapPeaks}
\title{
	Concatenate Peaks within a User-specified Distance
}
\description{
	Part of the \code{\link{waveseq}} workflow. This function reads in a file containing
	peak locations and a user-defined distance (gap*winsize) and concatenates peaks separated
	by at most 'gap' windows. This is useful for data having broad enrichment patterns as a 
	single large enrichment region can be called as several disjoint peaks. Analyzes peaks
	chromosome-by-chromosome.
}
\usage{
gapPeaks(peaks, outfile, gap, winsize = 200)
}
\arguments{
  \item{peaks}{
	Data frame or file containing peak locations. Must be in bedGraph format.
  }
  \item{outfile}{
	Output file containing gapped peaks (default=NA).
  }
  \item{gap}{
	Maximum number of non-significant windows separating peaks that are to be concatenated.
  }
  \item{winsize}{
	Window size used for obtaining peaks(default = 200 bp).
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
If \code{outfile} is not 'NA', returns data frame with peak 
locations and reads in four columns: <chr> <start> <end> <reads>.
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
}
\examples{
library(WaveSeqR)

# load peak list into memory
data(gabp_peaks)

# get file names
WS.path <- system.file(package = "WaveSeqR")
gap.file <- file.path(WS.path,"data","gabp_valouev2008_chr22-15m-30m_peaks_gap2.txt")

# concatenate peaks separated by 2 non-significant windows
#gapPeaks(peaks=peaks.file, outfile=gap.file, gap=2,winsize=200)
gapPeaks(peaks=gabp_peaks, outfile=gap.file, gap=2,winsize=200)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}