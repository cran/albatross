\name{feemparafac}
\alias{feemparafac}
\alias{plot.feemparafac}
\alias{coef.feemparafac}
\alias{fitted.feemparafac}
\alias{residuals.feemparafac}
\alias{reorder.feemparafac}
\alias{rescale} % re-exporting is tricky business
\alias{rescale.feemparafac}
\alias{print.feemparafac}
\title{
  Compute PARAFAC on a FEEM cube object and access the results
}
\description{
  \code{feemparafac} forwards its arguments to
  \code{\link[multiway]{parafac}} from the \CRANpkg{multiway} package,
  optionally rescales the result and attaches a few attributes.
  Resulting objects of class \code{feemparafac} can be accessed using
  methods presented below.
}
\usage{
  feemparafac(
    X, \dots, const = rep('nonneg', 3),
    rescale = 3, retries = 10, subset = TRUE, envir = NULL
  )
  \method{plot}{feemparafac}(x, type = c("image", "lines"), \dots)
  \method{coef}{feemparafac}(
    object, type = c(
      "all", "scores", "loadings", "emission", "excitation", "samples"
    ), \dots
  )
  \method{fitted}{feemparafac}(object, \dots)
  \method{residuals}{feemparafac}(object, \dots)
  \method{reorder}{feemparafac}(x, neworder, like, \dots)
  \method{rescale}{feemparafac}(x, mode, newscale, absorb, like, \dots)
}
\arguments{
  \item{X}{
    A FEEM cube object. The per-sample factors will be multiplied by
    \code{attr(X, 'scales')} stored in it.

    If \code{envir} is \code{NULL} (by default), this should be just a
    value. If \code{envir} is given, this should be a name of the value
    to \code{\link[base]{get}} from the environment.
  }
  \item{\dots}{\describe{
    \item{feemparafac}{Passed as-is to \code{\link[multiway]{parafac}}.}
    \item{plot.feemparafac}{
      Passed as-is to \CRANpkg{lattice} functions
      \code{\link[lattice]{levelplot}} and
      \code{\link[lattice]{xyplot}}.
    }
    \item{reorder.feemparafac, rescale.feemparafac}{
      Forwarded to the respective \CRANpkg{multiway} functions.
    }
    \item{coef.feemparafac, fitted.feemparafac, residuals.feemparafac}{
      No other parameters are allowed.
    }
  }}
  \item{const}{
    A character vector of length 3 specifying the constraints for all
    modes of \code{X}, passed to \code{\link[multiway]{parafac}}.
    Defaults to non-negativity. See \code{\link[CMLS]{const}} for more
    information.
  }
  \item{rescale}{
    Rescale the resulting factors to leave all the variance in the given
    mode: emission, excitation, or sample (default). Set to \code{NA}
    to disable.
  }
  \item{retries}{
    Retry for given number of tries until
    \code{\link[multiway]{parafac}} returns a successfully fitted model
    or stops due to the iteration number limit. Raise a fatal error if
    all tries were unsuccessful.
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
    emission loadings, \eqn3{\mathbf{B}}{B}{<b>B</b>} corresponding to
    fluorescence excitation loadings, and \eqn3{\mathbf{C}}{C}{<b>C</b>}
    corresponding to the scores of the components in different samples,
    the following plots can be produced:
    \describe{
      \item{image}{
        Plot the factors (\dQuote{loadings}) as a series of pseudo-colour
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
        respectively, as a \code{\link[base]{data.frame}} with three
        columns, the first one (named \code{wavelength} or
        \code{sample}) containing the independent variable
        (\eqn3{\lambda_\mathrm{em}}{lambda.em}{&lambda;<sub>em</sub>} /
        \eqn3{\lambda_\mathrm{ex}}{lambda.ex}{&lambda;<sub>ex</sub>} /
        sample name or number),
        the second one (named \code{value}) containing the values and
        the third one (named \code{factor}) containing the factor numbers.
      }
      \item{scores}{Same as samples.}
      \item{loadings}{
        Same as \dQuote{emission} and \dQuote{excitation} combined using
        \code{\link[base]{rbind}}, with a fourth column (\code{mode})
        added, naming the kinds of loadings.
      }
      \item{all}{
        A list with names \dQuote{emission}, \dQuote{excitation},
        \dQuote{samples} containing results of
        \code{coef(object, "emission")}, \code{coef(object, "excitation")},
        \code{coef(object, "samples")}, respectively.
      }
    }
  }
  \item{neworder}{
    A permutation of integers between \eqn{1} and \code{ncol(x$A)} (the
    number of components) specifying the new order of factors.
    Forwarded to \code{\link[multiway]{reorder.parafac}}. Incompatible
    with the \code{like} argument.
  }
  \item{like}{
    A \code{feemparafac} object.

    In \code{reorder}, the factors in \code{x} will be reordered to
    match the factors in \code{like} according to the smallest of the
    cosine similarities (\code{\link[multiway]{congru}}) for the
    emission and excitation wavelengths.

    In \code{rescale}, every factor matrix
    \eqn3{\mathbf{A}}{A}{<b>A</b>} from \code{x} that was specified in
    \code{mode} will be multiplied by scaling factors
    \eqn3{\mathbf{c}}{c}{<b>c</b>} so in order to minimise
    \eqn3{
      \| \mathbf{A}_\mathtt{like}
         - \mathrm{diag}(\mathbf c) \mathbf{A}_\mathtt{x} \|^2
    }{||like[[A]] - diag(c) \\\%*\\\% x[[A]]||^2}{
      ||<b>A</b><sub><code>like</code></sub>
        - <b>c</b><b>A</b><sub><code>x</code></sub>||<sup>2</sup>}.
  }
  \item{mode}{
    The modes to rescale, with \code{"A"}, \code{"B"}, \code{"C"}
    corresponding to emission, excitation, and samples, respectively.
    When \code{like} is specified, defaults to \code{c("A", "B")}.
    Forwarded to \code{\link[multiway]{rescale.parafac}}.
  }
  \item{newscale}{
    The desired root-mean-square for each column of the modes being
    rescaled.
    Forwarded to \code{\link[multiway]{rescale.parafac}}. Incompatible
    with the \code{like} argument.
  }
  \item{absorb}{
    The mode that should absorb the inverse rescaling coefficients.
    When \code{like} is specified, defaults to \code{"C"}.
    Forwarded to \code{\link[multiway]{rescale.parafac}}.
  }
}
\details{
  \code{feemparafac} tries hard to guarantee the convergence flag to be
  \code{0} (normal convergence) or \code{1} (iteration number limit
  reached), but never \code{2} (a problem with the constraints). A fatal
  error is raised if repeated runs of \code{\link[multiway]{parafac}} do
  not return a (semi-)successfully fitted model.

  After the PARAFAC decomposition is calculated, the scores are
  multiplied by the \code{scales} attribute of the \code{X} object,
  making them represent the cube with scaling undone. Use
  \code{\link{feemscale}(remember = FALSE)} if you don't want to undo
  the scaling.

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
  on \code{\link[base]{save}} or \code{\link[base]{load}}, so storing
  \code{X} in an environment and passing it to multiple invocations of
  \code{feemparafac} will save a lot of memory when the results are
  serialised together.

  \code{plot.feemparafac} provides sane defaults for \CRANpkg{lattice}
  options such as \code{xlab}, \code{ylab}, \code{as.table},
  \code{auto.key}, \code{type}, \code{cuts}, \code{col.regions}, but
  they can be overridden.
}
\value{
  \item{feemparafac}{
    An object of classes \code{feemparafac} and \code{parafac} with the
    following attributes added: \describe{
      \item{cube}{A copy of the \code{X} argument.}
      \item{subset}{A copy of the \code{subset} argument.}
      \item{envir}{A copy of the \code{envir} argument.}
    }
    \code{\link[base]{rownames}} are added from the original data cube
    to the \code{A}, \code{B}, \code{C} components of the list returned
    by \code{\link[multiway]{parafac}}.

    Use \code{\link{feemcube}} on the return value to access the original
    data cube.
  }
  \item{plot.feemparafac}{
    A \CRANpkg{lattice} plot object. Its \code{print} or \code{plot}
    method will draw the plot on an appropriate plotting device.
  }
  \item{coef.feemparafac}{
    A \code{\link[base]{data.frame}} or a list of them (only if
    \code{type} is \dQuote{all}). See the description of the \code{type}
    argument for more information.
  }
  \item{fitted.feemparafac}{
    A \code{\link{feemcube}} object comparable to \code{X} as it was
    decomposed by \code{\link[multiway]{parafac}}, ignoring the scaling.
  }
  \item{resid.feemparafac}{
    A \code{\link{feemcube}} object equal to \eqn3{
      \mathbf{X} - \hat{\mathbf{X}}
    }{X - fitted(X)}{<b>X</b> - <b>X&#x0302;</b>}, with an additional
    class \code{feem.resid} set. Objects of this class are plotted with
    a different default palette, see \code{\link{plot.feem.resid}}.
  }
}
\references{
  \reference{Bro1997}
}

\seealso{
  \code{\link[multiway]{parafac}} class structure;
  \code{\link{write.openfluor}}, \code{\link{feemcube}} for methods
  specific to values returned from this function.
}
\examples{
  data(feems)
  cube <- feemscale(feemscatter(cube, rep(14, 4)), na.rm = TRUE)
  (factors <- feemparafac(cube, nfac = 3))
  plot(factors, 'image')
  plot(factors, 'line')
  head(coef(factors, 'loadings'))
  str(coef(factors, 'all'))
  str(feemcube(factors)) # original cube is retained
  plot(fitted(factors))
  plot(resid(factors))
}
\keyword{hplot}
\keyword{method}
