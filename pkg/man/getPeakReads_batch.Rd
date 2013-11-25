\name{getPeakReads_batch}
\alias{getPeakReads_batch}
\title{
	Get Reads Within a List of Specified Locations in Batch
}
\description{
	Extracts reads falling within list of putative peaks from
	a list of padded bedGraph files and outputs to specified file.
	Uses the \code{getPeakReads} function.
}
\usage{
getPeakReads_batch(peaks, indir, design, outfile)
}
\arguments{
  \item{peaks}{
	File containing peak locations in first three columns: chromosome, start, end
  }
  \item{indir}{
	Input directory containing padded bedGraph files.
  }
  \item{design}{
	Tab-delimited file containing two columns: <file name>  <ID>. Files are assumed to be in
	\code{indir}.
  }
  \item{outfile}{
	Output file.
  }
}
\details{
	%%  ~~ If necessary, more details than the description above ~~
}
\value{
	%%
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
	\code{peaks} must be obtained as output by WaveSeqR.
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
	\code{\link{getPeakReads}}
}
\examples{
 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}