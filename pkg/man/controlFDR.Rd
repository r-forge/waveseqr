\name{controlFDR}
\alias{controlFDR}
\title{
	FDR Estimation with Control
}
\description{
	Estimates FDR for peaks in a ChIP sample by comparing it with a control data set. 
	\code{controlFDR} consists of three steps: \cr
	(1) Compare ChIP reads with control reads within peaks using \code{\link{compareWithControl}()}. \cr
	(2)	Adjust p-values for multiple testing using \code{\link{p.adjust}}. \cr
	(3)	Print peaks and p-values to file. \cr
}
\usage{
controlFDR(peaks, chip, control, outfile, winsize = 200,
			minreads = 8, binom.sided = "two.sided", adj.method = "fdr")
}
%- maybe also 'usage' for other objects documented here.
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
  \item{binom.sided}{
	One-sided or two-sided binomial test, if control is not NA (default = "two.sided")
  }
  \item{adj.method}{
	Method for multiple testing adjustment. Can use any methods accepted by \code{\link{p.adjust}}
	(default = "fdr").
  }
}
\value{
	This function writes peaks to a file after FDR estimation in the presence of control. This
	file has with the following fields:
	\item{chromosome}{Chromosome name}
	\item{start}{Start position of peak (1-based indexing)}
	\item{end}{End position of peak}
	\item{ChIP.reads}{
		Number of ChIP reads in the peak (normalized to per 
		million mapped reads)}	
	\item{Control.reads}{
		Number of control reads in the peak (normalized to per
		million mapped reads)}
	\item{Fold.Change}{
		Fold-change calculated as \code{ChIP.reads}/\code{Control.reads}}
	\item{p.value}{P-value for the peak}
	\item{adj.p.value}{P-values adjusted for multiple comparisons using \code{adj.method}}
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
	Apratim mitra
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{p.adjust}}, \code{\link{binom.test}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ distribution }
