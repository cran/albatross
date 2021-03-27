\name{plot.feem}
\alias{plot.feem}
\alias{plot.feemcube}
\title{Plot a FEEM object}
\description{
  Plot a 2D fluorescence intensity surface as a false-colour image.
}
\usage{
\method{plot}{feem}(
  x,
  xlab = quote(lambda[em] * ", nm"), ylab = quote(lambda[ex] * ", nm"),
  cuts = 128, col.regions = marine.colours(256), ...
)
\method{plot}{feemcube}(
  x,
  xlab = quote(lambda[em] * ", nm"), ylab = quote(lambda[ex] * ", nm"),
  cuts = 128, col.regions = marine.colours(256), as.table = TRUE, ...
)
}
\arguments{
  \item{x}{
    An FEEM object.
  }
  \item{xlab}{
    The x-axis label for the plot, with a sane default.
  }
  \item{ylab}{
    The y-axis label for the plot, with a sane default.
  }
  \item{cuts}{
    The number of distinct levels the intensity would be divided into,
    areas between them assigned different colours.
  }
  \item{col.regions}{
    The palette to take the colours from, a character vector of R colour
    specifications.
  }
  \item{as.table}{
    Whether to draw the panels left to right, top to bottom. (Otherwise
    they are drawn left to right, bottom to top.)
  }
  \item{\dots}{
    Passed as-is to \code{\link{levelplot}}.
  }
}
\value{
  A \pkg{lattice} plot object. Its \code{print} or \code{plot} method
  will draw the plot on an appropriate plotting device.
}

\seealso{
  \code{\link{levelplot}}
}
\examples{
  plot(feem(matrix(1:42/42, nrow = 7), 320 + 1:7, 300 + 1:6))
}
\keyword{hplot}
\keyword{method}
