# set environment variables and check perl installation
WS.perlpath <- ""
WS.path <- ""

.onLoad <- function(libname, pkgname)
{
  WS.path <- system.file(package="WaveSeqR")
  local.perlpath <- file.path("exec")
  
  WS.perlpath <- file.path(WS.path, local.perlpath)
  
  assign("WS.perlpath", WS.perlpath, envir=topenv())
  assign("WS.path", WS.path, envir=topenv())

  perlOutput <- try(system("perl --version", intern=TRUE, ignore.stderr=TRUE), silent=TRUE)
  if (inherits(perlOutput, "try-error")) {
    stop("perl not found! Must be installed for WaveSeqR to be operative.")
  }
}
