\name{feem}
\alias{feem}
\alias{feem.matrix}
\alias{feem.data.frame}
\title{
  Create a fluorescence excitation-emission matrix object
}
\description{
  Functions to create fluorescence excitation-emission matrix objects from
  R matrices coupled with excitation and emission wavelengths or three-column
  \code{data.frame}s containing
  \eqn{(\lambda_\mathrm{em}, \lambda_\mathrm{ex}, I)}{%
  (\lambda.em, \lambda.ex, intensity)} tuples.
}
\usage{
  feem(x, \dots)
  \method{feem}{matrix}(x, emission, excitation, scale = 1, \dots)
  \method{feem}{data.frame}(
    x, scale = 1, emission = 'emission',
    excitation = 'excitation', intensity = 'intensity', \dots
  )
}
\arguments{
  \item{x}{
    The matrix or three-column \code{data.frame} to convert into a FEEM object.

    If converting a matrix, its rows should correspond to different
    fluorescence emission wavelengths specified in \code{emission} argument;
    conversely, its columns should correspond to excitation wavelengths
    specified in \code{excitation} argument.

    If converting a \code{data.frame}, it should have exactly three columns,
    one for emission wavelengths, excitation wavelength, and intensity. The
    names of the columns are expected to be \code{'emission'},
    \code{'excitation'}, and \code{'intensity'}, respectively, but can be
    overridden using namesake arguments.
  }
  \item{emission}{
    If conversing a matrix, this should be a vector of emission wavelengths,
    each wavelength corresponding to a row of the matrix.

    If converting a \code{data.frame}, this optional argument specifies the
    name of the column containing the emission wavelengths.
  }
  \item{excitation}{
    If conversing a matrix, this should be a vector of excitation wavelengths,
    each wavelength corresponding to a column of the matrix.

    If converting a \code{data.frame}, this optional argument specifies the
    name of the column containing the excitation wavelengths.
  }
  \item{intensity}{
    If converting a \code{data.frame}, this optional argument specifies the
    name of the column containing the fluorescence intensities.
  }
  \item{scale}{
    The scale value of a EEM is preserved through the analysis procedure to
    make it easy to undo the scaling after runnining PARAFAC. If the EEM
    has been pre-multiplied prior to creating the FEEM object, you can set
    the multiplier here.
  }
  \item{\dots}{
    Extra arguments besides those specified above are not allowed.
  }
}
\value{
  A FEEM object is a matrix with the following attributes added:
  \item{emission}{
    Fluorescence emission wavelengths corresponding to the rows of the
    matrix, nm.
  }
  \item{excitation}{
    Fluorescence excitation wavelengths corresponding to the columns of
    the matrix, nm.
  }
  \item{dimnames}{
    Dimension names, copies of information above. Used only for
    presentation purposes.
  }
  \item{scale}{
    Scale factor, preserved through the analysis, which may be used
    later to undo the scaling. Initially \eqn{1}.
  }
}

\seealso{
  FEEM methods: \code{\link{plot.feem}}, \code{\link{as.data.frame.feem}},
  \code{\link{[.feem}}, \code{\link{feemife}}, \code{\link{feemscale}},
  \code{\link{feemscatter}}.
}
\examples{
  feem(matrix(1:40, ncol = 8), 1:5, 1:8)
  feem(
    data.frame(x = 1:10, y = 21:30, z = 31:40),
    emission = 'x', excitation = 'y', intensity = 'z'
  )
}
\keyword{methods}
\keyword{utilities}
