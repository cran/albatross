\name{plot.feemsplithalf}
\alias{plot.feemsplithalf}
\title{Plot results of split-half analysis}
\description{
  Present the results of repeatedly splitting the dataset into halves
  and fitting PARAFAC models to each half.
}
\usage{
  \method{plot}{feemsplithalf}(x, kind = c("tcc", "factors"), \dots)
}
\arguments{
  \item{x}{
    An object returned by \code{\link{feemsplithalf}}.
  }
  \item{kind}{
    \describe{
      \item{\code{"tcc"}}{
        Display statistics of between-half TCCs for different numbers
        of components. The smallest TCC is chosen between emission- and
        excitation-mode values, but otherwise they are not aggregated.
        Components with the same number have the same colour.
      }
      \item{\code{"factors"}}{
        Plot the resulting factors themselves on panels per each number
        of components and each mode (emission or excitation). Components
        with the same number have the same colour.
      }
    }
  }
  \item{\dots}{
    Passed as-is to \code{\link{xyplot}}.
  }
}
\details{
  The function provides sane defaults for \code{\link{xyplot}} parameters
  \code{xlab}, \code{ylab}, \code{as.table}, but they can be overridden.
}
\value{
  A \pkg{lattice} plot object. Its \code{print} or \code{plot} method
  will draw the plot on an appropriate plotting device.
}
\seealso{
  \code{\link{feemsplithalf}}, \code{\link{xyplot}}.
}
\examples{\donttest{
  data(feems)
  cube <- feemscale(
    feemscatter(feemcube(feems, FALSE), rep(24, 4)), na.rm = TRUE
  )
  # will take time
  sh <- feemsplithalf(cube, 2:4, random = 5, const = rep('nonneg', 3))
  plot(sh)
  plot(sh, 'factors')
}}
\keyword{hplot}
\keyword{method}
