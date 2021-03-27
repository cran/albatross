\name{feems}
\alias{feems}
\alias{absorp}
\docType{data}
\title{
  Fluorescence excitation-emission matrices and absorbance spectra
}
\description{
  This dataset contains: \itemize{
    \item Twelve FEEMs with anti-Stokes zone missing and scattering
    signal not handled, captured at different wavelength grids and with
    some points missing.
    \item Twelve absorbance spectra for purposes of IFE correction.
  }
}
\usage{data("feems")}
\format{\describe{
  \item{feems}{
    A named list of 12 \code{\link{feem}} objects containing fluorescence
    data measured at different wavelength grids (excitation wavelengths
    between \eqn{230} nm and \eqn{500} or \eqn{550} nm; emission
    wavelengths between \eqn{240} nm and \eqn{600} or \eqn{650} nm).
    Intensity at \eqn3{
      \lambda_\mathrm{em} < \lambda_\mathrm{ex} + 10\;\mathrm{nm}
    }{\lambda.em < \lambda.ex + 10 nm}{
      &lambda;<sub>em</sub> &lt; &lambda;<sub>ex</sub> + 10
      <span style="font-style: normal;">nm</span>
    } is not measured.
  }
  \item{absorp}{
    A named list of 12 2-column matrices containing absorbance spectra
    measured between \eqn{230} and \eqn{650} nm in \eqn{1} cm cells.
  }
}}
\examples{
  data(feems)
  plot(feems$a)
  plot(absorp$a)
}
\keyword{datasets}
