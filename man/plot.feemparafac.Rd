\name{plot.feemparafac}
\alias{plot.feemparafac}
\title{Plot components of fitted PARAFAC model}
\description{
  Plot the factors of a fitted PARAFAC model as functions of emission
  and excitation wavelengths.
}
\usage{
  \method{plot}{feemparafac}(x, type = c("image", "lines"), ...)
}
\arguments{
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
  \item{\dots}{
    Passed as-is to \pkg{lattice} functions \code{\link{levelplot}} and
    \code{\link{xyplot}}.
  }
}
\details{
  The function provides sane defaults for \pkg{lattice} options such as
  \code{xlab}, \code{ylab}, \code{as.table}, \code{auto.key},
  \code{type}, \code{cuts}, \code{col.regions}, but they can be
  overridden.
}
\value{
  A \pkg{lattice} plot object. Its \code{print} or \code{plot} method
  will draw the plot on an appropriate plotting device.
}

\seealso{
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
