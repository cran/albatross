\name{feemcube}
\alias{feemcube}
\alias{feemcube.array}
\alias{feemcube.list}
\title{Build a data cube of FEEMs}
\description{
  This function builds tagged 3-dimensional arrays of fluorescence
  excitation-emission spectra. Given a list of FEEM objects, it can
  determine the range of their wavelengths. Otherwise, the object is
  created from the supplied numeric array and vectors of wavelengths
  and sample names.
}
\usage{
  feemcube(x, \dots)
  \method{feemcube}{list}(x, all.wavelengths, \dots)
  \method{feemcube}{array}(x, emission, excitation, scales, names = NULL, \dots)
}
\arguments{
  \item{x}{
    A list of FEEM objects, possibly named.

    Alternatively, a numeric array.
  }
  \item{all.wavelengths}{
    Logical, a flag specifying whether to include wavelengths not present
    in \emph{all} of the samples. If \code{FALSE}, only those wavelength
    present in all of the samples are included.
  }
  \item{emission}{
    Numeric vector of emission wavelengths. Should correspond to the
    first dimension of the array \code{x}.
  }
  \item{excitation}{
    Numeric vector of excitation wavelengths. Should correspond to the
    second dimension of the array \code{x}.
  }
  \item{scales}{
    Numeric vector of scale factors corresponding to the spectra in the
    array. Should correspond to the third dimension of the array
    \code{x}. If missing, assumed to be all \eqn{1}.
  }
  \item{names}{
    Character vector of names of the samples. Should correspond to the
    third dimension of the array \code{x}.
  }
  \item{\dots}{
    Additional arguments besides those specified above are not allowed.
  }
}
\details{
  \code{feemcube.list} can be used to build FEEM data cubes from lists
  of FEEM objects even if their wavelength grids do not exactly match.
  The missing wavelengths may be set to \code{NA}
  (\code{all.wavelengths = TRUE}) or omitted from the cube
  (\code{all.wavelengths = FALSE}). See \code{\link{feemgrid}} if you
  need to adjust the wavelength grid of a list of EEMs before making it
  into a FEEM cube.
}
\value{
  A FEEM data cube is a numeric three-dimensional array with the
  following attributes:
  \item{emission}{
    Fluorescence emission wavelengths corresponding to the first
    dimension of the array, nm.
  }
  \item{excitation}{
    Fluorescence excitation wavelengths corresponding to the second
    dimension of the array, nm.
  }
  \item{dimnames}{
    Dimension names, copies of information above. Used only for
    presentation purposes.
  }
  \item{scales}{
    Scale factors of the samples, corresponding to the third dimension
    of the array. Assumed to be \eqn{1} missing if not specified by the
    user.
  }
}
\seealso{
  FEEM cube methods: \code{\link{[.feemcube}}, \code{\link{plot.feemcube}},
  \code{\link{as.data.frame.feemcube}}, \code{\link{as.list.feemcube}},
  \code{\link{feemife}}, \code{\link{feemscale}}, \code{\link{feemscatter}}.
}
\examples{
  # array form
  feemcube(array(1:24, 4:2), 1:4, 1:3)
  # list form
  feemcube(replicate(2, feem(matrix(1:6, 2), 1:2, 1:3), FALSE), TRUE)
}
\keyword{method}
