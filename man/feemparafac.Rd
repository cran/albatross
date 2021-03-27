\name{feemparafac}
\alias{feemparafac}
\alias{plot.feemparafac}
\alias{coef.feemparafac}
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
  \method{coef}{feemparafac}(
    object, type = c(
      "all", "scores", "loadings", "emission", "excitation", "samples"
    ), \dots
  )
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
    \item{coef.feemparafac}{
      No other parameters are allowed.
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
  \item{x, object}{
    An object returned by \code{\link{feemparafac}}.
  }
  \item{type}{
    Given a fitted PARAFAC model:

    \deqn3{
      X_{i,j,k} = \sum_r A_{i,r} B_{j,r} C_{k,r}
    }{X[i,j,k] = sum(A[i,] * B[j,] * C[k,])}{
      X<sub>i,j,k</sub> = &#x2211;<sub>r</sub>
        A<sub>i,r</sub> B<sub>j,r</sub> C<sub>k,r</sub>
    }

    With \eqn3{\mathbf{A}}{A}{<b>A</b>} corresponding to fluorescence
    emission spectra, \eqn3{\mathbf{B}}{B}{<b>B</b>} corresponding to
    fluorescence excitation spectra, and \eqn3{\mathbf{C}}{C}{<b>C</b>}
    corresponding to the scores of the components in different samples,
    the following plots can be produced:
    \describe{
      \item{image}{
        Plot the factors (\dQuote{loadings}) as a series of false-colour
        images of outer products
        \eqn3{
          \mathbf{a}_r \times \mathbf{b}_r^\top
        }{A[,r] \\\%o\\\% B[,r]}{
          <b>a</b><sub>r</sub> &times; <b>b</b><sub>r</sub><sup>
            <span style="font-style: normal;">T</span>
          </sup>
        }
      }
      \item{lines}{
        Plot the factors \eqn3{\mathbf a_r}{A[,r]}{<b>a</b><sub>r</sub>}
        and \eqn3{\mathbf b_r}{B[,r]}{<b>b</b><sub>r</sub>} as functions
        of wavelengths, with each pair of factors on a different panel.
      }
    }

    Fitted PARAFAC coefficients can be returned in the following forms:
    \describe{
      \item{emission, excitation, samples}{
        Return the contents of \eqn3{\mathbf{A}}{A}{<b>A</b>},
        \eqn3{\mathbf{B}}{B}{<b>B</b>} or \eqn3{\mathbf{C}}{C}{<b>C</b>},
        respectively, as a \code{\link{data.frame}} with three
        columns, the first one (named \code{wavelength} or \code{sample})
        containing the independent variable
        (\eqn3{\lambda_\mathrm{em}}{lambda.em}{&lambda;<sub>em</sub>} /
        \eqn3{\lambda_\mathrm{ex}}{lambda.ex}{&lambda;<sub>ex</sub>} /
        sample name or number),
        the second one (named \code{value}) containing the values and
        the third one (named \code{factor}) containing the factor numbers.
      }
      \item{scores}{Same as samples.}
      \item{loadings}{
        Same as \dQuote{emission} and \dQuote{excitation} combined using
        \code{\link{rbind}}, with a fourth column (\code{mode}) added,
        naming the kinds of loadings.
      }
      \item{all}{
        A list with names \dQuote{emission}, \dQuote{excitation},
        \dQuote{samples} containing results of
        \code{coef(object, "emission")}, \code{coef(object, "excitation")},
        \code{coef(object, "samples")}, respectively.
      }
    }
  }
}
\details{
  \code{feemparafac} tries hard to guarantee the convergence flag to be
  \code{0} (normal convergence) or \code{1} (iteration number limit
  reached), but never \code{2} (a problem with the constraints). A fatal
  error is raised if repeated runs of \code{\link{parafac}} do not
  return a (semi-)successfully fitted model.

  The \code{output} option is fixed to \code{"best"} value. Obtaining
  a list of alternative solutions can therefore be achieved by running:

  \code{replicate(\var{n}, feemparafac(..., nstart = 1), simplify = FALSE)}

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
    \code{\link{rownames}} are added from the original data cube to
    the \code{A}, \code{B}, \code{C} components of the list returned
    by \code{\link{parafac}}.

    If \code{envir} attribute is not \code{NULL}, use \code{\link{get}}
    together with the \code{cube} attribute to access the FEEM cube.
    Remember to use the \code{subset} attribute to access the samples that
    took part in the decomposition. Normally, this is only used to save
    memory when performing cross-validation; otherwise,
    \code{attr(pf, 'cube')} is enough.
  }
  \item{plot.feemparafac}{
    A \pkg{lattice} plot object. Its \code{print} or \code{plot} method
    will draw the plot on an appropriate plotting device.
  }
  \item{coef.feemparafac}{
    A \code{\link{data.frame}} or a list of them (only if \code{type}
    is \dQuote{all}). See the description of the \code{type} argument
    for more information.
  }
}
\references{
  \bibentry{'Article',
    author = person('Rasmus', 'Bro'),
    title = 'PARAFAC. Tutorial and applications',
    journal = 'Chemometrics and Intelligent Laboratory Systems',
    volume = 38,
    number = 2,
    pages = '149-171',
    year = 1997,
    doi = '10.1016/S0169-7439(97)00032-4'
  }
}

\seealso{
  \code{\link{parafac}} for the \code{parafac} class structure;
  \code{\link{fitted.feemparafac}}, \code{\link{residuals.feemparafac}},
  \code{\link{write.openfluor}} for methods specific to values returned
  from this function.
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
  head(coef(factors, 'loadings'))
  str(coef(factors, 'all'))
}
\keyword{hplot}
\keyword{method}
