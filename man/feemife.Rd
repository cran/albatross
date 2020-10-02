\name{feemife}
\alias{feemife}
\alias{feemife.feem}
\alias{feemife.feemcube}
\alias{feemife.list}
\title{Absorbance-based inner filter effect correction}
\description{
  Use absorbance data to correct inner-filter effect in FEEM objects
  and collections of them.
}
\usage{
  feemife(x, \dots)
  \method{feemife}{feem}(x, absorbance, abs.path = 1, \dots)
  \method{feemife}{feemcube}(x, absorbance, abs.path, \dots)
  \method{feemife}{list}(x, absorbance, abs.path, \dots)
}
\arguments{
  \item{x}{
    A FEEM object, a FEEM data cube, or a list of them.
  }
  \item{absorbance}{
    If \code{x} is a FEEM object: a two-column matrix-like object suitable
    for \code{\link{xy.coords}} containing the absorbance spectrum of the
    sample: the wavelengths in the first column and the unitless
    absorbance values in the second column.

    Otherwise, this could be a list of such objects or a multi-column
    matrix-like objects. If \code{x} contains names of the samples
    (is a named list or had \code{names} specified when calling
    \code{\link{feemcube}}), \code{absorbance} is a named list or has
    named columns, and all samples from \code{x} can be looked up in
    \code{absorbance}, results of this lookup are used. If name lookup
    fails but (given \eqn{N}-sample \code{x}) \code{absorbance} has
    exactly \eqn{N+1} columns or is an \eqn{N}-element list, absorbances
    are supposed to be present in the same order as the samples in
    \code{x}.
  }
  \item{abs.path}{
    If \code{x} is a FEEM object, a number specifying the length of the
    optical path used when measuring the absorbance, cm.

    Otherwise, a named vector containing the names from \code{x}, or a
    vector of exactly same length as the number of FEEMs in \code{x}:
    same lookup rules apply as for \code{absorbance} argument.

    If not set, assumed to be \eqn{1}.
  }
  \item{\dots}{
    No parameters besides those described above are allowed.
  }
}
\details{
  The formula used is:

  \deqn{
    I_{corr}(\lambda_{em}, \lambda_{ex}) = %
    I_{orig}(\lambda_{em}, \lambda_{ex}) 10^{ %
      \frac{A(\lambda{em}) + A(\lambda{ex})}{2 L_{abs}} %
    }
  }{I.corrected = I.orig * 10^((A.em + A.ex) / (2 * abs.path))}
}
\value{
  An object of the same kind as \code{x}, with inner filter effect
  corrected.
}
\references{
  J.R. Lakowicz, Principles of Fluorescence Spectroscopy, 3rd ed.,
  Springer US, 2006.
  \url{https://www.springer.com/la/book/9780387312781}

  D.N. Kothawala, K.R. Murphy, C.A. Stedmon, G.A. Weyhenmeyer,
  L.J. Tranvik, Inner filter correction of dissolved organic matter
  fluorescence: Correction of inner filter effects, Limnology and
  Oceanography: Methods. 11 (2013) 616-630.
  \doi{10.4319/lom.2013.11.616}
}
\examples{
  data(feems)

  dataset <- feemcube(feems, FALSE)
  str(dataset)
  str(absorp)
  plot(feemife(dataset,absorp) / dataset)
}
\keyword{method}
