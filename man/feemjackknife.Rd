\name{feemjackknife}
\alias{feemjackknife}
\alias{plot.feemjackknife}
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
}
\arguments{
  \item{cube}{
    A \code{\link{feemcube}} object.
  }
  \item{progress}{
    Set to \code{FALSE} to disable the progress bar.
  }
  \item{x}{
    An object returned by \code{\link{feemjackknife}}.
  }
  \item{kind}{\describe{
    \item{\code{"estimations"}}{
      Plot the loadings from every leave-one-out model.
    }
    \item{\code{"RIP"}}{
      Produce a Resample Influence Plot, i.e. plot mean squared
      difference between loadings in overall and leave-one-out models
      against mean squared residuals in leave-one-out models.
    }
    \item{\code{"IMP"}}{
      Produce an Identity Match Plot, i.e. plot scores in leave-one-out
      models against scores in overall model.
    }
  }}
  \item{\dots}{\describe{
    \item{feemjackknife}{
      Passed as-is to \code{\link{feemparafac}} and, eventually, to
      \code{\link{parafac}}
    }
    \item{plot.feemjackknife}{
      When \code{kind} is \code{"RIP"} or \code{"IMP"}, pass a \code{q}
      argument to specify the quantile of residual values (for RIP) or
      absolute score differences (IMP) above which sample names (or
      numbers) should be plotted.

      Remaining arguments are passed as-is to \code{\link{xyplot}}.
    }
  }}
}
\details{
  The function takes each sample out of the dataset, fits a PARAFAC
  model without it, then fits the outstanding sample to the model with
  emission and excitation factors fixed.

  The individual leave-one-out models are reordered according to best
  Tucker's congruence coefficient match and rescaled by minimising \eqn{%
    || \mathbf A \, \mathrm{diag}(\mathbf s_\mathrm A) - %
      \mathbf A^\mathrm{orig} ||^2 %
  }{sum((A * s.A - A.orig)^2)} and \eqn{ %
    || \mathbf{B} \, \mathrm{diag}(\mathbf s_\mathrm B) - %
      \mathbf B^\mathrm{orig} ||^2 %
  }{sum((B * s.B - B.orig)^2)} over \eqn{\mathbf s_\mathrm A}{s.A} and
  \eqn{\mathbf s_\mathrm B}{s.B}, subject to \eqn{%
    \mathrm{diag}(\mathbf s_\mathrm A) \cdot %
    \mathrm{diag}(\mathbf s_\mathrm B) \cdot %
    \mathrm{diag}(\mathbf s_\mathrm C) = \mathbf I}{%
    s.C = 1 / (s.A * s.B) %
  }, to make them comparable.

  Once the models are fitted, resample influence plots and identity
  match plots can be produced from resulting data to detect outliers.

  It is recommended to fully name the parameters to be passed to
  \code{\link{feemparafac}} to avoid problems.

  \code{plot.feemjackknife} provides sane defaults for \code{\link{xyplot}}
  parameters \code{xlab}, \code{ylab}, \code{scales}, \code{as.table}, but
  they can be overridden.
}
\value{\describe{
  \item{feemjackknife}{
    An list of class \code{feemjackknife} containing the following
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
    A \pkg{lattice} plot object. Its \code{print} or \code{plot} method
    will draw the plot on an appropriate plotting device.
  }
}}
\references{
  J. Riu, R. Bro, Jack-knife technique for outlier detection and
  estimation of standard errors in PARAFAC models, Chemometrics and
  Intelligent Laboratory Systems. 65 (2003) 35-49.
  \doi{10.1016/S0169-7439(02)00090-4}
}

\seealso{
  \code{\link{feemparafac}}
}
\examples{\donttest{
  data(feems)
  cube <- feemscale(
    feemscatter(feemcube(feems, FALSE), rep(24, 4))[1:30*6, 1:9*6,],
    na.rm = TRUE
  )
  # takes a long time
  jk <- feemjackknife(cube, nfac = 3, const = rep('nonneg', 3))
  plot(jk)
  plot(jk, 'IMP')
  plot(jk, 'RIP')
  plot(jk$leaveone[[1]])
}}
\keyword{method}
\keyword{hplot}
