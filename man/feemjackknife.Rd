\name{feemjackknife}
\alias{feemjackknife}
\alias{plot.feemjackknife}
\alias{coef.feemjackknife}
\title{Jack-knife outlier detection in PARAFAC models}
\description{
  Perform leave-one-out fitting + validation of PARAFAC models on a
  given FEEM cube.
}
\usage{
  feemjackknife(cube, ..., progress = TRUE)
  \method{plot}{feemjackknife}(
    x, kind = c('estimations', 'RIP', 'IMP'), \dots
  )
  \method{coef}{feemjackknife}(
    object, kind = c('estimations', 'RIP', 'IMP'), \dots
  )
}
\arguments{
  \item{cube}{
    A \code{\link{feemcube}} object.
  }
  \item{progress}{
    Set to \code{FALSE} to disable the progress bar.
  }
  \item{x, object}{
    An object returned by \code{\link{feemjackknife}}.
  }
  \item{kind}{
    Chooses what to plot (when called as \code{plot(\dots)}) or return
    as a \code{\link[base]{data.frame}} (when called as
    \code{coef(\dots)}):
    \describe{
      \item{estimations}{
        Produce the loadings from every leave-one-out model.
      }
      \item{RIP}{
        Produce a Resample Influence Plot, i.e. mean squared difference
        between loadings in overall and leave-one-out models plotted
        against mean squared residuals in leave-one-out models.
      }
      \item{IMP}{
        Produce an Identity Match Plot, i.e. scores in leave-one-out
        models plotted against scores in the overall model.
      }
    }
  }
  \item{\dots}{\describe{
    \item{feemjackknife}{
      Passed as-is to \code{\link{feemparafac}} and, eventually, to
      \CRANpkg{multiway} function \code{\link[multiway]{parafac}}.
    }
    \item{plot.feemjackknife}{
      When \code{kind} is \dQuote{RIP} or \dQuote{IMP}, pass a \code{q}
      argument to specify the quantile of residual values (for RIP)
      or absolute score differences (IMP) above which sample names (or
      numbers) should be plotted. Default value for \code{q} is \eqn{0.9}.

      Remaining arguments are passed as-is to
      \code{\link[lattice]{xyplot}}.
    }
    \item{coef.feemjackknife}{
      No further parameters are allowed.
    }
  }}
}
\details{
  The function takes each sample out of the dataset, fits a PARAFAC
  model without it, then fits the outstanding sample to the model with
  emission and excitation factors fixed:

  \deqn3{
    \hat{\mathbf{c}} =
      (\mathbf{A} \ast \mathbf{B})^{+} \times \mathrm{vec}(\mathbf{X})
  }{c.hat = pseudoinverse(krprod(A, B)) \\\%*\\\% vec(X)}{
    <b>c&#x0302;</b> = (
      <b>A</b> &ast; <b>B</b>
    )<sup>+</sup> &times;
    <span style="font-style: normal;">vec</span>(<b>X</b>)
  }

  The individual leave-one-out models (fitted loadings
  \eqn3{\mathbf A}{A}{<b>A</b>}, \eqn3{\mathbf B}{B}{<b>B</b>} and scores
  \eqn3{\mathbf C}{C}{<b>C</b>}) are reordered according to best Tucker's
  congruence coefficient match and rescaled by minimising \eqn3{
    || \mathbf A \, \mathrm{diag}(\mathbf s_\mathrm A) -
      \mathbf A^\mathrm{orig} ||^2
  }{sum((A * s.A - A.orig)^2)}{
    ||<b>A</b>
      <span style="font-style: normal;">diag</span>(<b>s</b><sub>A</sub>)
      - <b>A</b><sup>orig</sup>||<sup>2</sup>
  } and \eqn3{
    || \mathbf{B} \, \mathrm{diag}(\mathbf s_\mathrm B) -
      \mathbf B^\mathrm{orig} ||^2
  }{sum((B * s.B - B.orig)^2)}{
    ||<b>B</b>
      <span style="font-style: normal;">diag</span>(<b>s</b><sub>B</sub>)
      - <b>B</b><sup>orig</sup>||<sup>2</sup>
  } over \eqn3{\mathbf s_\mathrm A}{s.A}{<b>s</b><sub>A</sub>} and
  \eqn3{\mathbf s_\mathrm B}{s.B}{<b>s</b><sub>B</sub>}, subject to \eqn3{
    \mathrm{diag}(\mathbf s_\mathrm A) \times
    \mathrm{diag}(\mathbf s_\mathrm B) \times
    \mathrm{diag}(\mathbf s_\mathrm C) = \mathbf I
  }{s.A * s.B * s.C = 1}{
    <span style="font-style: normal;">diag</span>(<b>s</b><sub>A</sub>)
    &times;
    <span style="font-style: normal;">diag</span>(<b>s</b><sub>B</sub>)
    &times;
    <span style="font-style: normal;">diag</span>(<b>s</b><sub>C</sub>)
    = <b>I</b>
  }, to make them comparable.

  Once the models are fitted, resample influence plots and identity
  match plots can be produced from resulting data to detect outliers.

  To conserve memory, \code{feemjackknife} puts the user-provided
  \code{cube} in an environment and passes it via \code{envir} and
  \code{subset} options of \code{\link{feemparafac}}. This means that,
  unlike in \code{\link{feemparafac}}, the \code{cube} argument has
  to be a \code{\link{feemcube}} object and passing \code{envir} and
  \code{subset} options to \code{feemjackknife} is not supported. It
  is recommended to fully name the parameters to be passed to
  \code{\link{feemparafac}} to avoid problems.

  \code{plot.feemjackknife} provides sane defaults for
  \code{\link[lattice]{xyplot}} parameters \code{xlab}, \code{ylab},
  \code{scales}, \code{as.table}, but they can be overridden.
}
\value{\describe{
  \item{feemjackknife}{
    A list of class \code{feemjackknife} containing the following
    entries:\describe{
    \item{overall}{
      Result of fitting the overall \code{cube} with
      \code{\link{feemparafac}}.
    }
    \item{leaveone}{
      A list of length \code{dim(cube)[3]} containing the reduced dataset
      components. Every \code{\link{feemparafac}} object in the list has
      an additional \code{Chat} attribute containing the result of
      fitting the excluded spectrum back to the loadings of the reduced
      model.
    }
  }}

  \item{plot.feemjackknife}{
    A \CRANpkg{lattice} plot object. Its \code{print} or \code{plot} method
    will draw the plot on an appropriate plotting device.
  }

  \item{coef.feemjackknife}{
    A \code{\link[base]{data.frame}} containing various columns,
    depending on the value of the \code{kind} argument:

    \describe{
      \item{estimations}{\describe{
        \item{loading}{Values of the loadings.}
        \item{mode}{
          The axis of the loadings, \dQuote{Emission} or \dQuote{Excitation}.
        }
        \item{wavelength}{
          Emission or excitation wavelength the loading values correspond to.
        }
        \item{factor}{The component number.}
        \item{omitted}{
          The sample (name if \code{cube} had names, integer if it didn't)
          that was omitted to get the resulting loading values.
        }
      }}

      \item{RIP}{\describe{
        \item{msq.resid}{
          Mean squared residual value for the model with a given sample
          omitted.
        }
        \item{Emission}{
          Mean squared difference in emission mode loadings between the
          overall model and the model with a given sample omitted.
        }
        \item{Excitation}{
          Mean squared difference in excitation mode loadings between
          the overall model and the model with a given sample omitted.
        }
        \item{omitted}{
          The sample (name if \code{cube} had names, integer if it didn't)
          that was omitted from a given model.
        }
      }}

      \item{IMP}{\describe{
        \item{score.overall}{Score values for the overall model.}
        \item{score.predicted}{
          Score values estimated from the loadings of the model missing
          a given sample.
        }
        \item{factor}{The component number.}
        \item{omitted}{
          The sample (name if \code{cube} had names, integer if it didn't)
          that was omitted from a given model.
        }
      }}
    }
  }
}}
\references{
  \reference{Riu2003}
}

\seealso{
  \code{\link{feemparafac}}, \code{\link{feemcube}}
}
\examples{\donttest{
  data(feems)
  cube <- feemscale(feemscatter(cube, rep(14, 4)), na.rm = TRUE)
  # takes a long time
  jk <- feemjackknife(cube, nfac = 3)
  # feemparafac methods should be able to use the environment and subset
  plot(jk$leaveone[[1]])
  plot(jk)
  plot(jk, 'IMP')
  plot(jk, 'RIP')
  head(coef(jk))
}}
\keyword{method}
\keyword{hplot}
