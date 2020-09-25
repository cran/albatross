\name{feemjackknife}
\alias{feemjackknife}
\title{Jack-knife outlier detection in PARAFAC models}
\description{
  Perform leave-one-out fitting + validation of PARAFAC models on a
  given FEEM cube.
}
\usage{
  feemjackknife(cube, ...)
}
\arguments{
  \item{cube}{
    A \code{\link{feemcube}} object.
  }
  \item{\dots}{
    Passed as-is to \code{\link{feemparafac}} and, eventually, to
    \code{\link{parafac}}.
  }
}
\details{
  The function takes each sample out of the dataset, fits a PARAFAC
  model without it, then fits the outstanding sample to the model with
  emission and excitation factors fixed.

  The individual leave-one-out models are reordered according to best
  Tucker's congruence coefficient match and rescaled so that:

  \deqn{\sum_r \left( %
    \sum_i (S_{1,r} A_{i,r} - A^\mathrm{orig}_{i,r})^2 + %
    \sum_j (S_{2,r} B_{j,r} - B^\mathrm{orig}_{j,r})^2   %
    \right) \rightarrow \min_\mathbf{S}
  }{%
    sum( %
      (A * S[1,] - A.orig)^2, %
      (B * S[2,] - B.orig)^2  %
    ) -> min over S %
  }

  subject to \eqn{%
    S_{3,r} = \frac{1}{S_{1,r} S_{2,r}} \; \forall r}{%
    S[3,] = 1 / apply(S[1:2,], 2, prod)%
  }, to make them more comparable.

  Once the models are fitted, resample influence plots and identity
  match plots can be produced from resulting data to detect outliers.

  It is recommended to fully name the parameters to be passed to
  \code{\link{feemparafac}} to avoid problems.
}
\value{
  An list of class \code{feemjackknife} containing the following
  entries:
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
}
\references{
  J. Riu, R. Bro, Jack-knife technique for outlier detection and
  estimation of standard errors in PARAFAC models, Chemometrics and
  Intelligent Laboratory Systems. 65 (2003) 35-49.
  \doi{10.1016/S0169-7439(02)00090-4}
}

\seealso{
  \code{\link{plot.feemjackknife}}, \code{\link{feemparafac}}
}
\examples{\donttest{
  data(feems)
  cube <- feemscale(
    feemscatter(feemcube(feems, FALSE), rep(24, 4)),
    na.rm = TRUE
  )
  # takes a long time
  jk <- feemjackknife(cube, nfac = 3, const = rep('nonneg', 3))
  plot(jk)
}}
