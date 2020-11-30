\name{feemscatter}
\alias{feemscatter}
\alias{feemscatter.feem}
\alias{feemscatter.list}
\alias{feemscatter.feemcube}
\title{Handle scattering signal in FEEMs}
\description{
  Remove or interpolate scattering signal in individual FEEM objects,
  FEEM cube objects, or lists of them.
}
\usage{
  feemscatter(x, \dots)
  \method{feemscatter}{list}(x, \dots, cl, progress = TRUE)
  \method{feemscatter}{feemcube}(x, \dots, cl, progress = TRUE)
  \method{feemscatter}{feem}(
    x, widths, method = c("omit", "pchip", "loess", "kriging"),
    add.zeroes = 30, Raman.shift = 3400, \dots
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An individual FEEM object, FEEM cube object, or a list of them, to
    handle the scattering signal in.
  }
  \item{widths}{
    A numeric vector of length 4 containing the widths (in nm) of the
    scattering signal, in the following order:
    \enumerate{
      \item Rayleigh scattering
      \item Raman scattering
      \item Rayleigh scattering, \eqn{2 \lambda}
      \item Raman scattering, \eqn{2 \lambda}
    }
    Set a width to \eqn{0} if you don't want to handle this particular
    kind of scattering signal.
  }
  \item{method}{
    A string choosing \emph{how} to handle the scattering signal:
    \describe{
      \item{\code{"omit"}}{Replace it with \code{NA}.}
      \item{\code{"pchip"}}{
        Interpolate it line-by-line using piecewise cubic Hermitean
        polynomials (\code{\link{pchip}}).
      }
      \item{\code{"loess"}}{
        Interpolate it by fitting a locally weighted polynomial surface
        (\code{\link{loess}}). In this case the remaining parameters are
        passed verbatim to \code{\link{loess}}, which may be used to set
        parameters such as \code{span}.
      }
      \item{\code{"kriging"}}{
        Interpolate it by means of ordinary or simple Kriging, as
        implemented in \pkg{pracma} function \code{\link{kriging}}.
        Pass a \code{type} argument to choose between the two methods.
      }
    }
  }
  \item{add.zeroes}{
    Set intensities at
    \eqn{%
      \lambda_\mathrm{em} < \lambda_\mathrm{ex}
      - \mathtt{add.zeroes}\:\mathrm{nm}%
    }{%
      \lambda.em < \lambda.ex - add.zeroes (nm)
    } to \eqn{0} unless they have been measured. Set to \code{NA} to
    disable this behaviour.
  }
  \item{Raman.shift}{
    Raman shift of the scattering signal of water,
    \eqn{\textrm{cm}^{-1}}{cm^-1}.
  }
  \item{\dots}{
    Passed verbatim from \code{feemscatter} generics to
    \code{feemscatter.feem}. If \code{"loess"} method is selected,
    remaining options are passed to \code{\link{loess}} (the \code{span}
    parameter is of particular interest there).
  }
  \item{cl}{
    If not \code{\link{missing}}, a \pkg{parallel} cluster object to
    run the scattering correction code on or \code{NULL} for the default
    cluster object registered via \code{\link{setDefaultCluster}}.
  }
  \item{progress}{
    Set to \code{FALSE} to disable the progress bar.
  }
}
\details{
  The \code{"pchip"} method is implemented as described in the article
  by Bahram et al., 2006 (see the reference below): each emission
  spectrum at different excitation wavelengths is considered one by one.
  The last excitation spectrum (highest emission wavelength) is
  interpolated first to handle the corner case of extrapolating
  \eqn{2 \lambda} scattering near the high end of the emission
  wavelength range. Zero-valued points are inserted at the lowest
  emission wavelength (unless the data is already present) to prevent
  the extrapolation from blowing up.

  The \code{"loess"} method feeds the whole FEEM except the area to be
  interpolated to \code{\link{loess}}, then asks it to predict the
  remaining part of the spectum.

  The \code{"kriging"} method is much more computationally expensive than
  the previous two, but, on some spectra, provides best results, not
  affected by artifacts resulting from line-by-line one-dimensional
  interpolation (\code{pchip}) or varying degrees of smoothness in
  different areas of the spectrum (\code{loess}).
}
\value{
  An object of the same kind (FEEM object / FEEM cube / list of them)
  with scattering signal handled as requested.
}
\references{
  M. Bahram, R. Bro, C. Stedmon, A. Afkhami, Handling of Rayleigh and Raman
  scatter for PARAFAC modeling of fluorescence data using interpolation,
  Journal of Chemometrics. 20 (2006) 99-105.
  \doi{10.1002/cem.978}

  Press, W. H., A. A. Teukolsky, W. T. Vetterling, and B. P. Flannery (2007).
  Numerical recipes: The Art of Scientific Computing (3rd Ed.). Cambridge
  University Press, New York, Sect. 3.7.4, pp. 144-147.
}

\seealso{
  \code{\link{feem}}, \code{\link{feemcube}}
}
\examples{
  data(feems)
  plot(x <- feemscatter(
    feems[[1]], widths = c(25, 25, 20, 20), method = "loess", Raman.shift = 3500,
    span = .03, control = loess.control(trace.hat = 'approximate')
  ))
}
\keyword{method}
