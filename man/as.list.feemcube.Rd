\name{as.list.feemcube}
\alias{as.list.feemcube}
\title{
  Transform a FEEM cube object into a list of FEEM objects
}
\description{
  Return a list of FEEM objects comprising it. Used internally in
  \code{.feemcube} methods of the package generics and in
  \code{\link{as.data.frame.feemcube}}, but may be useful elsewhere.
}
\usage{
  \method{as.list}{feemcube}(x, \dots)
}
\arguments{
  \item{x}{
    A FEEM cube object.
  }
  \item{\dots}{
    No arguments besides those specified above are allowed.
  }
}
\value{
  A named list of FEEM objects comprising \code{x}.
}

\seealso{
  \code{\link{as.list}}; \code{\link{feemcube}} and its \code{list}
  constructor.
}
\examples{
  str(as.list(feemcube(array(1:60, 3:5), 1:3, 1:4)))
}
\keyword{method}
