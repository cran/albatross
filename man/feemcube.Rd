\name{feemcube}
\alias{feemcube}
\alias{feemcube.array}
\alias{feemcube.list}
\alias{feemcube.feemparafac}
\alias{feemcube.feemsplithalf}
\alias{feemcube.feemjackknife}
\alias{feemcube.feemflame}
\alias{as.list.feemcube}
\title{Data cubes of fluorescence excitation-emission matrices}
\description{
  Given a list of \code{\link{feem}} objects or a 3-way array, build
  tagged 3-dimensional arrays of fluorescence excitation-emission
  spectra. Extract the data cube from the corresponding model objects.
  Transform the data cube into a list of \code{\link{feem}} objects.
}
\usage{
  feemcube(x, \dots)
  \method{feemcube}{list}(x, all.wavelengths, \dots)
  \method{feemcube}{array}(x, emission, excitation, scales, names = NULL, \dots)
  \method{feemcube}{feemparafac}(x, \dots)
  \method{feemcube}{feemsplithalf}(x, \dots)
  \method{feemcube}{feemjackknife}(x, \dots)
  \method{feemcube}{feemflame}(x, \dots)
  \method{as.list}{feemcube}(x, \dots)
}
\arguments{
  \item{x}{\describe{
    \item{feemcube}{
      A list of FEEM objects, possibly named, or a numeric array.

      Alternatively, a \code{\link{feemparafac}},
      \code{\link{feemjackknife}}, \code{\link{feemsplithalf}}, or a
      \code{\link{feemflame}} object.
    }
    \item{as.list.feemcube}{
      A \code{feemcube} object.
    }
  }}
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

  \code{feemcube.feemparafac}, \code{feemcube.jackknife}, and
  \code{feemcube.feemsplithalf} return the data cube originally passed
  to the corresponding functions.
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
    of the array. Assumed to be \eqn{1} if not specified by the user.
  }

  \code{as.list.feemcube}: A named list of FEEM objects comprising
  \code{x}.
}
\seealso{
  FEEM cube methods: \code{\link{[.feemcube}},
  \code{\link{plot.feemcube}}, \code{\link{as.data.frame.feemcube}},
  \code{\link{feemife}}, \code{\link{feemscale}},
  \code{\link{feemscatter}}.
}
\examples{
  # array form
  feemcube(
    array(1:24, c(4, 3, 2)), # 3-way array obtained elsewhere
    seq(340, 400, len = 4), seq(250, 300, len = 3) # wavelengths
  )
  # list form
  feemcube(
    replicate(2, feem( # list of feem objects
      matrix(1:6, 2), c(340, 400), c(250, 275, 300)
    ), FALSE),
    TRUE
  )
  str(as.list(feemcube(array(1:60, 3:5), 1:3, 1:4)))
}
\keyword{method}
