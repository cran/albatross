\name{feemparafac}
\alias{feemparafac}
\alias{plot.feemparafac}
\title{
  Compute PARAFAC on a FEEM cube object
}
\description{
  This function forwards its arguments to \code{\link{parafac}} from the
  \pkg{multiway} package, optionally rescales the result and attaches
  a few attributes.
}
\usage{
  feemparafac(
    X, \dots, rescale = 3, retries = 10, subset = TRUE, envir = NULL
  )
  \method{plot}{feemparafac}(x, type = c("image", "lines"), \dots)
}
\arguments{
  \item{X}{
    A FEEM cube object. The per-sample factors will be multiplied by
    \code{attr(X, 'scales')} stored in it.

    If \code{envir} is \code{NULL} (by default), this should be just a
    value. If \code{envir} is given, this should be a name of the value
    to \code{\link{get}} from the environment.
  }
  \item{\dots}{\describe{
    \item{feemparafac}{Passed as-is to \code{\link{parafac}}.}
    \item{plot.feemparafac}{
      Passed as-is to \pkg{lattice} functions \code{\link{levelplot}} and
      \code{\link{xyplot}}.
    }
  }}
  \item{rescale}{
    Rescale the resulting factors to leave all the variance in the given
    mode: emission, excitation, or sample (default). Set to \code{NA}
    to disable.
  }
  \item{retries}{
    Retry for given number of tries until \code{\link{parafac}} returns
    a successfully fitted model or stops due to the iteration number
    limit. Raise a fatal error if all tries were unsuccessful.
  }
  \item{subset}{
    An integer or logical vector choosing the samples from \code{X}, as
    in \code{feemparafac(X[,,subset], \dots)}. Defaults to the whole cube.
  }
  \item{envir}{
    An environment to look up \code{X} in.
  }
  \item{x}{
    An object returned by \code{\link{feemparafac}}.
  }
  \item{type}{
    Given a fitted PARAFAC model:

    \deqn{%
      X_{ijk} = \sum_r A_{ir} B_{jr} C_{kr}}{%
      X[i,j,k] = sum(A[i,] * B[j,] * C[k,])%
    }

    With \eqn{\mathbf{A}}{A} corresponding to fluorescence emission
    spectra, \eqn{\mathbf{B}}{B} corresponding to fluorescence
    excitation spectra, and \eqn{\mathbf{C}}{C} corresponding to the
    scores of the components in different samples, the following plots
    can be produced:
    \describe{
      \item{\code{"image"}}{
        Plot the factors (\dQuote{loadings}) as a series of false-colour
        images of outer products
        \eqn{\mathbf{A}_r \otimes \mathbf{B}_r}{A[,r] \%o\% B[,r]}.
      }
      \item{\code{"lines"}}{
        Plot the factors \eqn{\mathbf A_r}{A[,r]} and \eqn{%
          \mathbf B_r}{B[,r]} as functions of wavelengths, with each
        pair of factors on a different panel.
      }
    }
    More plot kinds may be added in the future.
  }
}
\details{
  The function tries hard to guarantee the convergence flag to be
  \code{0} (normal convergence) or \code{1} (iteration number limit
  reached), but never \code{2} (a problem with the constraints). A fatal
  error is raised if repeated runs of \code{\link{parafac}} do not
  return a (semi-)successfully fitted model.

  The \code{output} option is fixed to \code{"best"} value. Obtaining
  a list of alternative solutions can therefore be achieved by running:

  \code{replicate(\var{n}, feemparafac(..., nstart = 1))}

  The \code{subset} and \code{envir} options are useful to repeatedly
  perform PARAFAC on different subsets of the same FEEM cube, e.g. in
  jack-knifing or split-half analysis. Since \code{feemparafac} keeps
  a reference to the its \code{X} and \code{envir} arguments, the use
  of \code{subset} should ensure that the same FEEM cube is referenced
  from multiple \code{feemparafac} objects instead of creating copies
  of its subsets. Additionally, environment objects are not duplicated
  on \code{\link{save}} or \code{\link{load}}, so storing \code{X} in an
  environment and passing it to multiple invocations of \code{feemparafac}
  will save a lot of memory when the results are serialised together.

  \code{plot.feemparafac} provides sane defaults for \pkg{lattice} options
  such as \code{xlab}, \code{ylab}, \code{as.table}, \code{auto.key},
  \code{type}, \code{cuts}, \code{col.regions}, but they can be
  overridden.
}
\value{
  \item{feemparafac}{
    An object of classes \code{feemparafac} and \code{parafac} with the
    following attributes added: \describe{
      \item{cube}{A copy of the \code{X} argument.}
      \item{subset}{A copy of the \code{subset} argument.}
      \item{envir}{A copy of the \code{envir} argument.}
    }
    To access the cube that was factored, use \code{\link{get}} on
    \code{cube} attribute if \code{envir} attribute is not \code{NULL},
    then extract the samples in question using the \code{subset}
    attribute. Normally, this is only needed to save memory when
    performing cross-validation; otherwise, \code{attr(pf, 'cube')}
    is enough.
  }
  \item{plot.feemparafac}{
    A \pkg{lattice} plot object. Its \code{print} or \code{plot} method
    will draw the plot on an appropriate plotting device.
  }
}
\references{
  R. Bro, PARAFAC. Tutorial and applications, Chemometrics and Intelligent
  Laboratory Systems. 38 (1997) 149-171. \doi{10.1016/S0169-7439(97)00032-4}
}

\seealso{
  \code{\link{parafac}} for the \code{parafac} class structure;
  \code{\link{fitted.feemparafac}}, \code{\link{residuals.feemparafac}},
  for methods specific to values returned from this function.

  \code{\link{levelplot}}, \code{\link{xyplot}}.
}
\examples{
  data(feems)
  cube <- feemscale(
    feemscatter(
      feemcube(feems, FALSE)[(1:45)*4,(1:13)*4,],
      rep(24, 4)), na.rm = TRUE
  )
  factors <- feemparafac(cube, nfac = 2, const = rep('nonneg', 3))
  plot(factors, 'image')
  plot(factors, 'line')
}
\keyword{hplot}
\keyword{method}
