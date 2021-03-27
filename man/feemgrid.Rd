\name{feemgrid}
\alias{feemgrid}
\alias{feemgrid.feem}
\alias{feemgrid.feemcube}
\alias{feemgrid.list}
\title{
  Interpolate FEEMs on a given wavelength grid
}
\description{
  Use interpolation to change the wavelength grid of a single FEEM or
  unify the grid of a collection of them.
}
\usage{
  feemgrid(x, ...)
  \method{feemgrid}{feem}(
    x, emission, excitation,
    method = c("whittaker", "loess", "kriging", "pchip"), \dots
  )
  \method{feemgrid}{feemcube}(
    x, emission, excitation, \dots, progress = TRUE
  )
  \method{feemgrid}{list}(
    x, emission, excitation, \dots, progress = TRUE
  )
}
\arguments{
  \item{x}{
    A \code{\link{feem}} object, a \code{\link{feemcube}}, or a list of
    \code{\link{feem}} objects.
  }
  \item{emission, excitation}{
    Desired wavelength grid, as numeric vectors. Must be specified for
    a single FEEM.  If not specified for a collection of FEEMs, all
    wavelengths falling in the range of the intersection all wavelengths
    intervals are chosen.
  }
  \item{method}{
    Interpolation method, see \code{\link{feemscatter}} for details.
  }
  \item{\dots}{
    Passed from generics to \code{feemgrid.feem}, then to interpolation
    methods.  See \code{\link{feemscatter}} for details.
  }
  \item{progress}{
    Set to \code{FALSE} to disable the progress bar.
  }
}
\details{
  The algorithm doesn't know how to distinguish between \code{NA}s
  that haven't been measured and \code{NA}s that resulted from
  combining different wavelength grids, so it tries to interpolate all
  of them.  As a result, leaving large areas of the spectrum undefined
  (e.g. anti-Stokes area) is not recommended, since it would result in
  extrapolation and introduce strong artefacts.
}
\value{
  An object of the same kind (FEEM object / FEEM cube / list of them)
  with emission and excitation wavelengths as requested.
}

\seealso{
  \code{\link{feemscatter}}
}
\examples{
  data(feems)
  x <- feemscatter(feems$a, rep(25, 4))
  y <- feemgrid(x, seq(240, 600, 5), seq(230, 550, 10))
  plot(plot(x, main = 'Original'    ), split = c(1, 1, 2, 1), more = TRUE)
  plot(plot(y, main = 'Interpolated'), split = c(2, 1, 2, 1))
}
