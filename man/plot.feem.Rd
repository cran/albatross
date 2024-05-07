\name{plot.feem}
\alias{plot.feem}
\alias{plot.feemcube}
\alias{plot.feem.resid}
\title{Plot a FEEM object}
\description{
  Plot a 2D fluorescence intensity surface as a pseudo-colour image.
}
\usage{
\method{plot}{feem}(
  x, xlab, ylab, cuts = 128,
  col.regions = marine.colours(256), \dots
)
\method{plot}{feemcube}(
  x, xlab, ylab, cuts = 128,
  col.regions = marine.colours(256), as.table = TRUE, \dots
)
\method{plot}{feem.resid}(
  x, \dots, at, col.regions = diverging.colours(256)
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
    The palette to take the colours from, a character vector of \R
    colour specifications.
  }
  \item{at}{
    The breakpoints in the intensity values. For residual plots,
    automatically set in a symmetric manner, ignoring the \code{cuts}
    argument.
  }
  \item{as.table}{
    Whether to draw the panels left to right, top to bottom. (Otherwise
    they are drawn left to right, bottom to top.)
  }
  \item{\dots}{
    Passed as-is to \code{\link[lattice]{levelplot}}.
  }
}
\value{
  A \CRANpkg{lattice} plot object. Its \code{print} or \code{plot} method
  will draw the plot on an appropriate plotting device.
}

\seealso{
  \code{\link[lattice]{levelplot}}
}
\examples{
  plot(feem(matrix(1:42/42, nrow = 7), 320 + 1:7, 300 + 1:6))
}
\keyword{hplot}
\keyword{method}
