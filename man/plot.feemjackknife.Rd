\name{plot.feemjackknife}
\alias{plot.feemjackknife}
\title{Plot jack-knifing results}
\description{Plot the results of PARAFAC jack-knifing of a dataset.}
\usage{
  \method{plot}{feemjackknife}(
    x, kind = c('estimations', 'RIP', 'IMP'), \dots
  )
}
\arguments{
  \item{x}{
    An object returned by \code{\link{feemjackknife}}.
  }
  \item{kind}{\describe{
    \item{\code{"estimations"}}{
      Plot the loadings from every leave-one-out model.
    }
    \item{\code{"RIP"}}{
      Produce a Resample Influence Plot, i.e. plot mean squared
      difference between loadings in overall and leave-one-out models
      against mean squared residuals in leave-one-out models.
    }
    \item{\code{"IMP"}}{
      Produce an Identity Match Plot, i.e. plot scores in leave-one-out
      models against scores in overall model.
    }
  }}
  \item{\dots}{
    When \code{kind} is \code{"RIP"} or \code{"IMP"}, pass a \code{q}
    argument to specify the quantile of residual values (for RIP) or
    absolute score differences (IMP) above which sample names (or
    numbers) should be plotted.

    Passed as-is to \code{\link{xyplot}}.
  }
}
\details{
  The function provides sane defaults for \code{\link{xyplot}} parameters
  \code{xlab}, \code{ylab}, \code{scales}, \code{as.table}, but they can
  be overridden.
}
\value{
  A \pkg{lattice} plot object. Its \code{print} or \code{plot} method
  will draw the plot on an appropriate plotting device.
}
\references{
  J. Riu, R. Bro, Jack-knife technique for outlier detection and
  estimation of standard errors in PARAFAC models, Chemometrics and
  Intelligent Laboratory Systems. 65 (2003) 35-49.
  \doi{10.1016/S0169-7439(02)00090-4}
}
\seealso{
  \code{\link{feemjackknife}}
}
\examples{\donttest{
  data(feems)
  cube <- feemscale(
    feemscatter(feemcube(feems, FALSE), rep(24, 4)), na.rm = TRUE
  )
  plot( # jack-knifing takes a long time
    jk <- feemjackknife(cube, nfac = 2, const = rep('nonneg', 3))
  )
  plot(jk, 'IMP')
  plot(jk, 'RIP')
}}
\keyword{method}
\keyword{hplot}
