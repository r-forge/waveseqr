\name{getConfidenceInterval}
\alias{getConfidenceInterval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get confidence intervals of Wavelet coefficients after Monte Carlo sampling
}
\description{
This function calculates the mean and standard deviation of wavelet coefficients
across chromosomes for each scale of the continuous wavelet transform after 
the Monte Carlo threshold estimation step.
}
\usage{
getConfidenceInterval(indir, exptname, mother = "morlet")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{indir}{
	Input directory containing threshold estimates
}
  \item{exptname}{
	Experiment name. This is a pattern used for getting a list of files from \code{indir}.
}
  \item{mother}{
	Mother wavelet used for the Monte Carlo threshold estimation. Default="morlet"
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
}
\keyword{utilities}