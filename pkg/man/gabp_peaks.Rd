\name{gabp_peaks}
\alias{gabp_peaks}
\docType{data}
\title{
Peak Calling results on GABP data
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Data frame containing the list of peaks detected by WaveSeqR
in a published GABP ChIP-Seq data set. The data corresponds
to the genomic region chr22:15000000-30000000 of the hg18 version 
of the human genome. Mother function was 'morlet' and gap size was 0.
}
\usage{data(gabp_peaks)}
\format{
  A data frame with 1590 observations on the following 4 variables.
  \describe{
    \item{\code{chr}}{Chromosome name}
    \item{\code{start}}{Start position}
    \item{\code{end}}{End position}
    \item{\code{reads}}{Number of reads}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
http://mendel.stanford.edu/sidowlab/downloads/quest/
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Valouev, A. et al (2008) \emph{Genome-wide analysis of transcription factor binding sites based on ChIP-Seq data.}.
Nat Methods. 2008 Sep;5(9):829-34.
}
\examples{
data(gabp_peaks)
## maybe str(gabp_peaks) ; plot(gabp_peaks) ...
}
\keyword{datasets}
