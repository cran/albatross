\name{feems}
\alias{feems}
\alias{absorp}
\alias{cube}
\docType{data}
\title{
  Synthetic fluorescence excitation-emission matrices and absorbance spectra
}
\description{
  This dataset consists of twelve fluorescence and absorbance spectra
  simulated from three trilinear components, with scattering signal
  added and divided by a correction factor to simulate inner filter
  effect.
}
\usage{data("feems")}
\format{\describe{
  \item{feems}{
    A named list of 12 \code{\link{feem}} objects containing
    fluorescence data measured with excitation wavelengths between
    \eqn{230} nm and \eqn{350} nm (with a step of \eqn{2} nm) and
    emission wavelengths between \eqn{240} nm and \eqn{435} (with a step
    of \eqn{5} nm).
  }
  \item{cube}{
    A 12-sample \code{\link{feemcube}} object consisting of of 32 by 10
    FEEMs measured at the same wavelength range as above with inner
    filter effect corrected.
  }
  \item{absorp}{
    A 12-element named list containing absorbance spectra measured between
    \eqn{230} and \eqn{450} nm in \eqn{1} cm cells. Each element of the
    list is a two-column matrix. The first column contains the
    wavelengths and the second column contains the absorbance values.
  }
}}
\examples{
  data(feems)
  plot(cube)
  plot(feems$a)
  matplot(
    absorp[[1]][,1],
    sapply(absorp, function(x) x[,2]),
    type = 'l', lty = 1
  )
}
\keyword{datasets}
