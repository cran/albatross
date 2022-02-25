\name{[.feem}
\alias{[.feem}
\alias{[<-.feem}
\title{Extract or replace parts of FEEM objects}
\description{
  Extract or replace parts of FEEM spectra. Returns FEEM objects unless
  dimensions should be dropped. When assigning from a FEEM object,
  requires wavelengths to match and warns if scale factors differ.
}
\usage{
  \method{[}{feem}(x, i, j, drop = TRUE)
  \method{[}{feem}(x, i, j) <- value
}
\arguments{
  \item{x}{
    A FEEM object.
  }
  \item{i, j}{
    Row and column indices, respectively. As in usual R subsetting
    (see \link{Extract}), may be integer, logical or character vectors,
    or missing.
  }
  \item{drop}{
    Coerce result to the lowest possible dimension (dropping the
    \code{feem} class if so).
  }
  \item{value}{
    An array-like object to assign values from. When assigning from
    FEEM objects, wavelengths are required to match and warnings are
    issued if scale factors don't match. Use vector subsetting (zero or
    one argument inside the brackets) to disable the check.
  }
}
\value{
  For \code{[}: If \code{drop} is \code{TRUE} and at least one of the
  index arguments chooses only one element along its axis, a named
  numeric vector. Otherwise, a FEEM object.

  For \code{[<-}: a FEEM object.
}

\seealso{
  \code{\link{feem}}, \code{\link{[.feemcube}}
}
\examples{
  (z <- feem(matrix(1:40, ncol = 8), 66 + 1:5, 99 + 1:8, 3))
  str(z[1:4, 1:2])
  str(z[1,, drop = TRUE])
  z[2:3, 4:5] <- feem(matrix(1:4, 2), 66 + 2:3, 99 + 4:5, 3)
  z
}
\keyword{array}
\keyword{method}
