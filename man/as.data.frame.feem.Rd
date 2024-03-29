\name{as.data.frame.feem}
\alias{as.data.frame.feem}
\alias{as.data.frame.feemcube}
\title{Transform a FEEM object into a data.frame}
\description{
  Transform a FEEM object from its matrix form accompanied by vectors of
  wavelengths into a three-column form consisting of \eqn3{
    (\lambda_\mathrm{em}, \lambda_\mathrm{ex}, I)
  }{(\lambda.em, \lambda.ex, intensity)}{
    (&lambda;<sub>em</sub>, &lambda;<sub>ex</sub>, I)
  } tuples, which could be useful for
  export or plotting with \CRANpkg{lattice} or \CRANpkg{ggplot2}.
}
\usage{
  \method{as.data.frame}{feem}(x, row.names = NULL, optional = FALSE, ...)
  \method{as.data.frame}{feemcube}(x, ...)
}
\arguments{
  \item{x}{
    A FEEM object, or a FEEM cube object.
  }
  \item{row.names}{
    Passed to \code{\link[base]{data.frame}}. If default of \code{NULL}
    is used, \code{\link[base]{data.frame}} will generate sequential
    integer \code{row.names}.
  }
  \item{optional}{
    This option is required for compatibility with
    \code{\link[base]{as.data.frame}} generic, but is ignored, since the
    column names are already syntactic and the row names are generated
    by \code{\link[base]{data.frame}} automatically by default.
  }
  \item{\dots}{
    Passed as-is to \code{\link[base]{data.frame}}.
  }
}
\details{
  Rows where intensity is \code{NA} are omitted from the output.
}
\value{
  A \code{\link[base]{data.frame}} containing three numeric columns:
  \item{emission  }{Emission wavelength, nm.}
  \item{excitation}{Excitation wavelength, nm.}
  \item{intensity }{
    Fluorescence intensity at
    \eqn3{
      (\lambda_\mathrm{em}, \lambda_\mathrm{ex})
    }{(\lambda.em, \lambda.ex)}{
      (&lambda;<sub>em</sub>, &lambda;<sub>ex</sub>)
    }
  }
  \item{sample}{
    For FEEM cube objects, the unique name of the sample possessing this
    tuple of values, a factor. If the original object didn't have any
    names, sequential integers are used instead. If the original object
    had non-unique names, sequence numbers are appended to them using
    \code{\link[base]{make.unique}}.
  }
}

\seealso{
  \code{\link{feem.data.frame}}
}

\examples{
  z <- feem(matrix(1:42, nrow = 7), 1:7, 1:6)
  head(as.data.frame(z))
}
