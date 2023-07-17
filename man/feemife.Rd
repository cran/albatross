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
  \method{feemife}{feemcube}(x, absorbance, abs.path, \dots, progress = FALSE)
  \method{feemife}{list}(x, absorbance, abs.path, \dots, progress = FALSE)
}
\arguments{
  \item{x}{
    A FEEM object, a FEEM data cube, or a list of them.
  }
  \item{absorbance}{
    If \code{x} is a FEEM object: a two-column matrix-like object
    containing the absorbance spectrum of the sample: the wavelengths
    in the first column and the unitless absorbance values in the
    second column.

    Otherwise, this could be a list of such objects or a multi-column
    matrix-like object. If \code{x} contains names of the samples
    (is a named list or had \code{names} specified when calling
    \code{\link{feemcube}}), \code{absorbance} is a named list or has
    named columns, and all samples from \code{x} can be looked up in
    \code{absorbance}, results of this lookup are used. If \code{x}
    or \code{absorbance} isn't named, but (given \eqn{N}-sample
    \code{x}) \code{absorbance} has exactly \eqn{N+1} columns or is an
    \eqn{N}-element list, absorbance spectra are assumed to be present
    in the same order as the samples in \code{x}. Otherwise, an error
    is raised.
  }
  \item{abs.path}{
    If \code{x} is a FEEM object, a number specifying the length of the
    optical path used when measuring the absorbance, cm.

    Otherwise, a named vector containing the names from \code{x}, or a
    vector of exactly same length as the number of FEEMs in \code{x}:
    same lookup rules apply as for \code{absorbance} argument.

    If not set, assumed to be \eqn{1}.
  }
  \item{progress}{
    Set to \code{TRUE} to enable a progress bar (implemented via
    \code{\link[utils]{txtProgressBar}}).
  }
  \item{\dots}{
    No parameters besides those described above are allowed.
  }
}
\details{
  If you receive errors alleging that some names don't match, but
  are absolutely sure that the absorbance spectra and path lengths are
  present in the same order as in \code{x}, remove the names from either
  of the objects, e.g. by passing \code{\link[base]{unname}(absorbance)}.

  The formula used to correct for inner filter effect is:

  \deqn3{
    I_\mathrm{corr}(\lambda_\mathrm{em}, \lambda_\mathrm{ex}) =
    I_\mathrm{orig}(\lambda_\mathrm{em}, \lambda_\mathrm{ex}) 10^{
      \frac{A(\lambda_\mathrm{em}) + A(\lambda_\mathrm{ex})}{2 L_\mathrm{abs}}
    }
  }{I.corrected = I.orig * 10^((A.em + A.ex) / (2 * abs.path))}{
    I<sub>corr</sub>(&lambda;<sub>em</sub>, &lambda;<sub>ex</sub>) =
    I<sub>orig</sub>(&lambda;<sub>em</sub>, &lambda;<sub>ex</sub>) &middot;
    10<sup>
      <sup>A(&lambda;<sub>em</sub>) +
        A(&lambda;<sub>ex</sub>)</sup>/<sub>2L<sub>abs</sub>
      </sub>
    </sup>
  }
}
\value{
  An object of the same kind as \code{x}, with inner filter effect
  corrected.
}
\references{
  \reference{Lakowicz2006}

  \reference{Kothawala2013}
}
\examples{
  data(feems)

  str(cube)
  str(absorp)
  plot(feemife(cube,absorp) / cube)
}
\keyword{method}
