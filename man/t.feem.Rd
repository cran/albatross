\name{t.feem}
\alias{t.feem}
\title{Transpose a FEEM}
\description{
  Transpose a FEEM object and return an unclassed matrix.
}
\usage{
  \method{t}{feem}(x)
}
\arguments{
  \item{x}{
    A \code{\link{feem}} object.
  }
}
\value{
  A \code{\link{matrix}} without a \code{\link{class}} attribute.
}
\details{
  Like in \code{\link{aperm}}, transposing makes sense for a matrix, but
  violates the \code{\link{feem}} class invariants. Therefore, the class
  attribute has to be removed.
}
\examples{
  t(feem(matrix(1:40, ncol = 8), 1:5, 1:8))
}
\keyword{methods}
\keyword{utilities}
