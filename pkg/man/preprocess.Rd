\name{preprocess}
\alias{preprocess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Data preprocessing
}
\description{
	Preprocesses an input bed file to generate summary bedgraph 
	and padded bedgraph files for use with WaveSeqR. Also generates
	chromosome information file.
}
\usage{
preprocess(infile, outdir, exptname, winsize = 200, fragsize = 150, redundancy = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infile}{
	BED file to be preprocessed (required)
}
  \item{outdir}{
	Output directory (required)
}
  \item{exptname}{
	Unique string to identify name of experiment (required)
}
  \item{winsize}{
	Window size for padded graph files (default=200)
}
  \item{fragsize}{
	Library size used for sequencing (default=150)
}
  \item{redundancy}{
	Logical value indicating whether redundant reads should be removed 
	before calculating summary counts (default=TRUE)
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
	Function messages are piped to a log file which can be found in the \code{log}
	folder of \code{outdir}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
