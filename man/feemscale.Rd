\name{feemscale}
\alias{feemscale}
\alias{feemscale.feem}
\alias{feemscale.feemcube}
\alias{feemscale.list}
\title{Rescale FEEM spectra to a given norm and remember the scale factor}
\description{
  Given a norm function (typically, standard deviation), scale the
  intensities in FEEM objects to it and optionally remember the scale
  factor.
}
\usage{
  feemscale(x, \dots)
  \method{feemscale}{feem}(x, norm = sd, remember = TRUE, \dots)
  \method{feemscale}{feemcube}(x, \dots)
  \method{feemscale}{list}(x, \dots)
}
\arguments{
  \item{x}{
    A FEEM object, a FEEM cube object, or a list of anything compatible
    with \code{feemscale} generic.
  }
  \item{norm}{
    A function taking a numeric matrix and returning its norm. Typically,
    \code{\link{sd}} or \code{\link{sumsq}}.
  }
  \item{remember}{
    Whether to remember the scale factor. If \code{FALSE}, the scale
    factor in the returned object is unchanged.
  }
  \item{\dots}{
    Passed as-is to \code{feemscale} and, eventually, to
    \code{feemscale.feem}. Also passed as-is to the norm function,
    which is useful for passing \code{na.rm = TRUE} to functions like
    \code{\link{sd}} or \code{\link{sumsq}}.
  }
}
\value{
  \code{feemscale.feem}: a FEEM object with intensities divided by scale
  factor (\code{norm(x)}) and its \code{scale} attribute multiplied by
  the scale factor.

  \code{feemscale.feemcube}: a FEEM cube built from FEEM objects scaled
  as described above.

  \code{feemscale.list}: a list consisting of results of
  \code{feemscale} generic applied to its elements.
}

\seealso{
  \code{\link{feem}}
}
\examples{
  feemscale(feem(matrix(1:42, 6), 1:6, 1:7))
}
\keyword{methods}
