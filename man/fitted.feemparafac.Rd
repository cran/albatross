\name{fitted.feemparafac}
\alias{fitted.feemparafac}
\alias{residuals.feemparafac}
\title{Extract fitted PARAFAC values or residuals}
\description{
  \code{fitted} calculates an approximation of a FEEM cube fitted by
  PARAFAC.

  \code{residuals} returns the difference between \code{fitted} and the
  original data as a FEEM cube.
}
\usage{
  \method{fitted}{feemparafac}(object, ...)
  \method{residuals}{feemparafac}(object, ...)
}
\arguments{
  \item{object}{
    An object returned by \code{\link{feemparafac}}.
  }
  \item{\dots}{
    No arguments besides those described above are allowed.
  }
}
\details{
  The output of \code{\link{fitted.parafac}} from \pkg{multiway} package
  is rescaled back according to saved scales (typically those from
  \code{\link{feemscale}}, e.g. \code{\link{sd}} norms of each spectrum)
  from the original cube in order to be comparable with it.

  The output of \code{residuals} is
  \eqn3{
    \mathbf{X} - \hat{\mathbf{X}}
  }{X - fitted(X)}{
    <b>X</b> - <b>X&#x0302;</b>
  }.
}
\value{
  A FEEM cube object.
}
\seealso{
  \code{\link{feemcube}}, \code{\link{fitted.parafac}},
  \code{\link{resid}}.
}
\examples{
  data(feems)
  cube <- feemscale(
    feemscatter(
      feemcube(feems, FALSE)[1:36*5, 1:11*5,],
      rep(24, 4)), na.rm = TRUE
  )
  factors <- feemparafac(cube, nfac = 2, const = rep('nonneg', 3))
  # calls plot.feemcube for estimated spectra
  plot(fitted(factors))
  plot(resid(factors))
}
\keyword{method}
