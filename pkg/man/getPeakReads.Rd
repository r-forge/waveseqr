\name{getPeakReads}
\alias{getPeakReads}
\title{
	Get Reads Within a List of Specified Locations
}
\description{
	Part of the \code{\link{compareWithControl}} function. Extracts reads
	falling within specified locations from a data frame obtained from
	a padded bedGraph file.
}
\usage{
getPeakReads(peak.mat, reads.mat, winsize = 200)
}
\arguments{
  \item{peak.mat}{
	Data frame containing peak locations. Has four columns: chromosome, start, end, reads
  }
  \item{reads.mat}{
	Data frame containing read counts from a padded bedGraph file. Has the same columns as
	\code{peak.mat}.
  }
  \item{winsize}{
	Window size used for obtaining peaks(default = 200 bp).
  }
}
\details{
	%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A vector of same length as the number of rows of \code{peak.mat} containing read 
	counts extracted from \code{read.mat}.
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
	\code{read.mat} must be obtained from a padded bedGraph file as output by the
	WaveSeq preprocessing pipeline.
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}