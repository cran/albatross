\name{feemsplithalf}
\alias{feemsplithalf}
\alias{plot.feemsplithalf}
\alias{print.feemsplithalf}
\title{Split-half analysis of PARAFAC models}
\description{
  This function validates PARAFAC with different numbers of components
  by means of splitting the data cube in halves, fitting PARAFAC to them
  and comparing the results [1].
}
\usage{
  feemsplithalf(cube, nfac, splits, random, \dots)
  \method{plot}{feemsplithalf}(x, kind = c("tcc", "factors"), \dots)
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
  \item{kind}{
    \describe{
      \item{\code{"tcc"}}{
        Display statistics of between-half TCCs for different numbers
        of components. The smallest TCC is chosen between emission- and
        excitation-mode values, but otherwise they are not aggregated.
        Components with the same number have the same colour.
      }
      \item{\code{"factors"}}{
        Plot the resulting factors themselves on panels per each number
        of components and each mode (emission or excitation). Components
        with the same number have the same colour.
      }
    }
  }
  \item{\dots}{\describe{
    \item{feemsplithalf}{
      Remaining options are passed to \code{\link{feemparafac}} and,
      eventually, to \code{\link{parafac}}.
    }
    \item{plot.feemsplithalf}{Passed as-is to \code{\link{xyplot}}.}
    \item{print.feemsplithalf}{No options are allowed.}
  }}
}
\details{
  Pass either \code{splits} or \code{random} parameter, but not both, as
  they are mutually exclusive.

  As the models are fitted, they are compared to the first model of the
  same number of factors (Tucker's congruence coefficient is calculated
  using \code{\link{congru}} for emission and excitation mode factors,
  then the smallest value of the two is chosen for the purposes of
  matching). The models are first reordered according to the best match
  by TCC value, then rescaled [3] by minimising \eqn{%
    || \mathbf A \, \mathrm{diag}(\mathbf s_\mathrm A) - %
      \mathbf A^\mathrm{orig} ||^2 %
  }{sum((A * s.A - A.orig)^2)} and \eqn{ %
    || \mathbf{B} \, \mathrm{diag}(\mathbf s_\mathrm B) - %
      \mathbf B^\mathrm{orig} ||^2 %
  }{sum((B * s.B - B.orig)^2)} over \eqn{\mathbf s_\mathrm A}{s.A} and
  \eqn{\mathbf s_\mathrm B}{s.B}, subject to \eqn{%
    \mathrm{diag}(\mathbf s_\mathrm A) \, %
    \mathrm{diag}(\mathbf s_\mathrm B) \, %
    \mathrm{diag}(\mathbf s_\mathrm C) = \mathbf I}{%
    s.C = 1 / (s.A * s.B) %
  }, to make them comparable.

  \code{plot.feemsplithalf} provides sane defaults for \code{\link{xyplot}}
  parameters \code{xlab}, \code{ylab}, \code{as.table}, but they can be
  overridden.
}
\value{\describe{
  \item{feemsplithalf, print.feemsplithalf}{
    An object of class \code{feemsplithalf}, containing named fields:
    \describe{
      \item{factors}{
        A \code{\link{list}} of \code{\link{feemparafac}} objects containing
        the factors of the halves. The list has dimensions, the first one
        corresponding to the halves (always 2), the second to different
        numbers of factors (as many as in \code{nfac}) and the third to
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
  }
  \item{plot.feemsplithalf}{
    A \pkg{lattice} plot object. Its \code{print} or \code{plot} method
    will draw the plot on an appropriate plotting device.
  }
}}
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
  \code{\link{feemparafac}}, \code{\link{parafac}}, \code{\link{congru}},
  \code{\link{feemcube}}.
}
\examples{\donttest{
  data(feems)
  cube <- feemscale(
    feemscatter(feemcube(feems, FALSE), rep(24, 4))[1:30*6, 1:9*6,],
    na.rm = TRUE
  )
  (sh <- feemsplithalf( # takes a long time
    cube, 2:4, splits = 4, # 4 splits => S4C6T3
    # the rest is passed to multiway::parafac
    const = rep('nonneg', 3) # setting ctol and maxit is recommended
  ))
  plot(sh)
  plot(sh, 'factors')
}}
\keyword{hplot}
\keyword{method}
