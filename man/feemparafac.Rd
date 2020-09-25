\name{feemparafac}
\alias{feemparafac}
\title{
  Compute PARAFAC on a FEEM cube object
}
\description{
  This function forwards its arguments to \code{\link{parafac}} from the
  \pkg{multiway} package, optionally rescales the result and attaches
  a few attributes.
}
\usage{
  feemparafac(X, \dots, rescale = 3, retries = 10)
}
\arguments{
  \item{X}{
    A FEEM cube object. The per-sample factors will be multiplied by
    \code{attr(X, 'scales')} stored in it.
  }
  \item{\dots}{
    Passed as-is to \code{\link{parafac}}.
  }
  \item{rescale}{
    Rescale the resulting factors to leave all the variance in the given
    mode: emission, excitation, or sample (default). Set to \code{NA}
    to disable.
  }
  \item{retries}{
    Retry for given number of tries until \code{\link{parafac}} returns
    a successfully fitted model or stops due to the iteration number
    limit. Raise a fatal error if all tries were unsuccessful.
  }
}
\details{
  The function tries hard to guarantee the convergence flag to be
  \code{0} (normal convergence) or \code{1} (iteration number limit
  reached), but never \code{2} (a problem with the constraints). A fatal
  error is raised if repeated runs of \code{\link{parafac}} do not
  return a (semi-)successfully fitted model.

  The \code{output} option is fixed to \code{"best"} value. Obtaining
  a list of alternative solutions can therefore be achieved by running:

  \code{replicate(\var{n}, feemparafac(..., nstart = 1))}
}
\value{
  An object of classes \code{feemparafac} and \code{parafac} with the
  \code{cube} attribute added containing a reference to the original
  FEEM cube object.
}
\references{
  R. Bro, PARAFAC. Tutorial and applications, Chemometrics and Intelligent
  Laboratory Systems. 38 (1997) 149-171. \doi{10.1016/S0169-7439(97)00032-4}
}

\seealso{
  \code{\link{parafac}} for the \code{parafac} class structure;
  \code{\link{fitted.feemparafac}}, \code{\link{residuals.feemparafac}},
  \code{\link{plot.feemparafac}} for methods specific to values returned
  from this function.
}
\examples{
  data(feems)
  cube <- feemscale(
    feemscatter(
      feemcube(feems, FALSE)[(1:45)*4,(1:13)*4,],
      rep(24, 4)), na.rm = TRUE
  )
  plot(
    factors <- feemparafac(
      cube, nfac = 2, const = rep('nonneg', 3)
    )
  )
}
