\name{[.feemcube}
\alias{[.feemcube}
\alias{[<-.feemcube}
\title{Extract or replace parts of FEEM cubes}
\description{
  Extract or replace single intensities, vectors of them, whole FEEM
  spectra or even data cubes or their parts from a FEEM cube.
}
\usage{
  \method{[}{feemcube}(x, i, j, k, drop = TRUE)
  \method{[}{feemcube}(x, i, j, k) <- value
}
\arguments{
  \item{x}{
    A FEEM cube object.
  }
  \item{i, j, k}{
    Row, column and sample indices, respectively. As usual, may be
    integer, logical or character vectors. Omitting a parameter results
    in choosing the whole axis.
  }
  \item{drop}{
    Coerce result to the lowest possible dimension (dropping the
    \code{feemcube} class).
  }
  \item{value}{
    An array-like object to assign values from. When assigning from
    FEEM or FEEM cube objects, wavelengths are required to match and
    warnings are issued if scale factors don't match. Use vector
    subsetting (zero or one argument inside the brackets) to disable the
    check.
  }
}
\value{
  For \code{[}: If choosing multiple values along each axis or
  \code{drop} is \code{FALSE}, a FEEM cube object. If choosing only one
  sample but multiple wavelengths, a FEEM object. Otherwise, a named
  numeric matrix or vector, depending on the dimensions chosen.

  For \code{[<-}: a FEEM cube object.
}
\seealso{
  \code{\link{feemcube}}, \code{\link{[.feem}}
}
\examples{
  z <- feemcube(array(1:385, c(5, 7, 11)), 1:5, 1:7, 1:11)
  str(z[1:4, 1:2, 1:2])
  z[2:3, 4:5, 3] <- feem(matrix(1:4, 2), 2:3, 4:5, 3)
  z[,,3]
}
\keyword{array}
\keyword{method}
