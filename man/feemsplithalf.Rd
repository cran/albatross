\name{feemsplithalf}
\alias{feemsplithalf}
\alias{print.feemsplithalf}
\title{Split-half analysis of PARAFAC models}
\description{
  This function validates PARAFAC with different numbers of components
  by means of splitting the data cube in halves, fitting PARAFAC to them
  and comparing the results [1].
}
\usage{
  feemsplithalf(cube, nfac, splits, random, \dots)
  \method{print}{feemsplithalf}(x, \dots)
}
\arguments{
  \item{cube}{
    A \code{\link{feemcube}} object.
  }
  \item{nfac}{
    An integer vector of numbers of factors to check.
  }
  \item{splits}{
    Number of parts to split the data cube into. Must be even. After
    splitting, all ways to recombine the parts into non-intersecting
    halves are enumerated [2], the halves are subjected to PARAFAC
    decomposition and compared against each other.

    The number of PARAFAC models fitted is \eqn{%
      {{\mathtt{splits}} \choose {\mathtt{splits}/2}} }{%
      choose(splits, splits/2)%
    }.

    Mutually incompatible with \code{random} parameter.
  }
  \item{random}{
    Number of times to shuffle the dataset, split into non-intersecting
    halves, fit a PARAFAC model to each of the halves and compare halves
    against each other.

    The number of PARAFAC models fitted is \eqn{2 \cdot \mathtt{random}}{%
      2 * random}.

    Mutually incompatible with \code{splits} parameter.
  }
  \item{x}{
    An object returned by \code{feemsplithalf}.
  }
  \item{\dots}{
    \code{feemsplithalf}: remaining options are passed to
    \code{\link{feemparafac}} and, eventually, to \code{\link{parafac}}.

    \code{print.feemsplithalf}: no options are allowed.
  }
}
\details{
  Pass either \code{splits} or \code{random} parameter, but not both, as
  they are mutually exclusive.

  As the models are fitted, they are compared to the first model of the
  same number of factors (Tucker's congruence coefficient is calculated
  using \code{\link{congru}} for emission and excitation mode factors,
  then the smallest value of the two is chosen for the purposes of
  matching). The models are first reordered according to the best match
  by TCC value, then rescaled [3] so that:

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
  }.
}
\value{
  An object of class \code{feemsplithalf}, containing named fields:
  \item{factors}{
    A \code{\link{list}} of \code{\link{feemparafac}} objects containing
    the factors of the halves. The list has dimensions, the first one
    corresponding to the halves (always 2), the second to different
    numbers of factors (as many in \code{nfac}) and the third to
    different groupings of the samples (depends on \code{splits} or
    \code{random}).
  }
  \item{tcc}{
    A named list containing arrays of Tucker's congruence coefficients
    between the halves. Each entry in the list corresponds to an element
    in the \code{nfac} argument. The dimensions of each array in the
    list correspond to, in order: the factors (1 to \var{nfac}), the
    modes (emission or excitation) and the groupings of the samples
    (depending on \code{splits} or \code{random}).
  }
  \item{nfac}{
    A copy of \code{nfac} argument.
  }
}
\references{\enumerate{
  \item W.S. DeSarbo, An Application of PARAFAC to a Small Sample
  Problem, Demonstrating Preprocessing, Orthogonality Constraints, and
  Split-Half Diagnostic Techniques (Appendix), Social Science Research
  Network, Rochester, NY, 1984.
  \url{https://papers.ssrn.com/abstract=2783446}

  \item K.R. Murphy, C.A. Stedmon, D. Graeber, R. Bro, Fluorescence
  spectroscopy and multi-way techniques. PARAFAC, Analytical Methods.
  5 (2013) 6557. \doi{10.1039/c3ay41160e}

  \item J. Riu, R. Bro, Jack-knife technique for outlier detection and
  estimation of standard errors in PARAFAC models, Chemometrics and
  Intelligent Laboratory Systems. 65 (2003) 35-49.
  \doi{10.1016/S0169-7439(02)00090-4}
}}

\seealso{
  \code{\link{plot.feemsplithalf}}, \code{\link{feemparafac}},
  \code{\link{parafac}}, \code{\link{congru}}, \code{\link{feemcube}}.
}
\examples{\donttest{
  data(feems)
  cube <- feemscale(
    feemscatter(feemcube(feems, FALSE), rep(24, 4)), na.rm = TRUE
  )
  (sh <- feemsplithalf( # takes a long time
    cube, 2:4, splits = 4, # 4 splits => S4C6T3
    # the rest is passed to multiway::parafac
    const = rep('nonneg', 3) # setting ctol and maxit is recommended
  ))
  plot(sh)
}}
\keyword{method}
