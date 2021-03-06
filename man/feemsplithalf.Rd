\name{feemsplithalf}
\alias{feemsplithalf}
\alias{plot.feemsplithalf}
\alias{print.feemsplithalf}
\alias{coef.feemsplithalf}
\title{Split-half analysis of PARAFAC models}
\description{
  This function validates PARAFAC with different numbers of components
  by means of splitting the data cube in halves, fitting PARAFAC to them
  and comparing the results [1].
}
\usage{
  feemsplithalf(
    cube, nfac, splits, random, groups, \dots, progress = TRUE
  )
  \method{plot}{feemsplithalf}(x, kind = c('tcc', 'factors'), \dots)
  \method{print}{feemsplithalf}(x, \dots)
  \method{coef}{feemsplithalf}(object, kind = c('tcc', 'factors'), \dots)
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

    The number of PARAFAC models fitted is \eqn{
      \mathtt{splits} \choose {\mathtt{splits}/2}
    }{choose(splits, splits/2)}.

    Mutually incompatible with the \code{random} parameter.
  }
  \item{random}{
    Number of times to shuffle the dataset, split into non-intersecting
    halves, fit a PARAFAC model to each of the halves and compare halves
    against each other.

    The number of PARAFAC models fitted is
    \eqn3{2 \cdot \mathtt{random}}{2 * random}{
      2 &middot; <code>random</code>}.

    Mutually incompatible with the \code{splits} parameter.
  }
  \item{groups}{
    Use this argument to preserve the ratios between the groups present
    in the original dataset when constructing the halves. If specified,
    must be a factor or an integer vector of length \code{dim(cube)[3]}
    (specifying the group each sample belongs to) or a list of them,
    i.e., a valid \code{f} argument to \code{\link{split}}. By default,
    samples are considered to form a single group.

    For the split-combine method (\code{splits}), each group must have
    at least \code{splits} elements; for best results, sizes of groups
    should be close to a multiple of \code{splits}. For the randomised
    split-half method (\code{random}), each group should have at least
    \eqn{2} elements.
  }
  \item{progress}{
    Set to FALSE to disable the progress bar.
  }
  \item{x, object}{
    An object returned by \code{feemsplithalf}.
  }
  \item{kind}{
    Chooses what type of data to return or plot:
    \describe{
      \item{tcc}{
        Statistics of between-half TCCs for different numbers of
        components. The smallest TCC is chosen between emission- and
        excitation-mode values, but otherwise they are not aggregated.

        When plotting, TCC values for the component with the same number
        have the same colour.
      }
      \item{factors}{
        The resulting loading values.

        When plotting, split the plot into panels per each number of
        components and each mode (emission or excitation). Components
        with the same number have the same colour.
      }
    }
  }
  \item{\dots}{\describe{
    \item{feemsplithalf}{
      Remaining options are passed to \code{\link{feemparafac}} and,
      eventually, to \code{\link{parafac}}. It is recommended to fully
      name the parameters instead of relying on partial or positional
      matching.
    }
    \item{plot.feemsplithalf}{Passed as-is to \code{\link{xyplot}}.}
    \item{print.feemsplithalf, coef.feemsplithalf}{
      No additional options are allowed.
    }
  }}
}
\details{
  As the models (loadings \eqn3{\mathbf A}{A}{<b>A</b>},
  \eqn3{\mathbf B}{B}{<b>B</b>} and scores \eqn3{\mathbf C}{C}{<b>C</b>})
  are fitted, they are compared to the first model of the same number
  of factors (Tucker's congruence coefficient is calculated using
  \code{\link{congru}} for emission and excitation mode factors, then the
  smallest value of the two is chosen for the purposes of matching). The
  models are first reordered according to the best match by TCC value,
  then rescaled [3] by minimising \eqn3{
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

  To perform stratified sampling on a real-valued variable (e.g. salinity,
  depth), consider binning samples into groups using \code{\link{cut}},
  perhaps after histogram flattening using \code{\link{ecdf}(x)(x)}. To
  determine the number of breaks, consider \code{\link{nclass.Sturges}}.

  To conserve memory, \code{feemsplithalf} puts the user-provided
  \code{cube} in an environment and passes it via \code{envir} and
  \code{subset} options of \code{\link{feemparafac}}. This means that,
  unlike in \code{\link{feemparafac}}, the \code{cube} argument has
  to be a \code{\link{feemcube}} object and passing \code{envir} and
  \code{subset} options to \code{feemsplithalf} is not supported.

  \code{plot.feemsplithalf} plots results of the split-half procedure
  (TCC or loading values depending on the \code{kind} argument)
  using \pkg{lattice} graphics. Sane defaults are provided for
  \code{\link{xyplot}} parameters \code{xlab}, \code{ylab},
  \code{as.table}, but they can be overridden.

  \code{print.feemsplithalf} displays a very short summary of the
  analysis, currently the minimum TCC value for each number of components.

  \code{coef.feemsplithalf} returns the Tucker's congruence
  coefficients resulting from the split-half analysis.
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
        list correspond to, in order: the factors (1 to \code{nfac[i]}), the
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
  \item{coef.feemsplithalf}{
    A \code{\link{data.frame}} containing various columns, depending on
    the value of the \code{kind} argument:

    \describe{
      \item{tcc}{\describe{
        \item{factor}{
          The factor (out of \code{nfac}) under consideration.
        }
        \item{tcc}{
          Tucker's congruence coefficient between a pair of matching
          components.  Out of two possible values (TCC between excitation
          loadings or emission loadings), the minimal one is chosen,
          because the same rule is used to find which components match
          when reordering them in a pair of models.
        }
        \item{test}{
          The sequence number for each pair of models in the split-half
          test, related to the third dimension of \code{object$factors}
          or \code{object$tcc}.  May be used to group values for plotting
          or aggregation.
        }
        \item{subset}{
          Consists of two-element lists containing indices of the samples
          in each half of the original cube.
        }
        \item{nfac}{
          The number of factors in the pair of models under consideration.
        }
      }}
      \item{factors}{\describe{
        \item{wavelength}{
          Emission and excitation wavelengths.
        }
        \item{value}{
          The values of the loadings.
        }
        \item{factor}{
          Number of the factor, \eqn{1} to \code{nfac}.
        }
        \item{mode}{
          The mode the loading value belongs to, \dQuote{Emission}
          or \dQuote{Excitation}.
        }
        \item{nfac}{
          Total number of factors.
        }
        \item{test}{
          Sequence number of a split-half test, indicating a given way
          to split the dataset in a group of splits with the same numbers
          of factors.
        }
        \item{half}{
          Number of the half, \eqn{1} or \eqn{2}.
        }
        \item{subset}{
          For every row, this is an integer vector indicating the
          subset of the original data cube that the loadings have been
          obtained from.
        }
      }}
    }
  }
}}
\references{\enumerate{
  \item \bibentry{'Article',
    author = person(c('Wayne', 'S.'), 'DeSarbo'),
    title = paste(
      'An Application of PARAFAC to a Small Sample Problem,',
      'Demonstrating Preprocessing, Orthogonality Constraints,',
      'and Split-Half Diagnostic Techniques (Appendix)'
    ),
    journal = 'Research Methods for Multimode Data Analysis',
    year = 1984,
    pages = '602-642',
    url = 'https://papers.ssrn.com/abstract=2783446'
  }

  \item \bibentry{'Article',
    author = c(
      person(c('Kathleen', 'R.'), 'Murphy'),
      person(c('Colin', 'A.'), 'Stedmon'),
      person('Daniel', 'Graeber'),
      person('Rasmus', 'Bro')
    ),
    title = 'Fluorescence spectroscopy and multi-way techniques. PARAFAC',
    journal = 'Analytical Methods',
    doi = '10.1039/c3ay41160e',
    volume = 5,
    year = 2013,
    pages = '6557-6566'
  }

  \item \bibentry{'Article',
    author = c(
      person('Jordi', 'Riu'),
      person('Rasmus', 'Bro')
    ),
    title = paste(
      'Jack-knife technique for outlier detection and',
      'estimation of standard errors in PARAFAC models'
    ),
    journal = 'Chemometrics and Intelligent Laboratory Systems',
    volume = 65,
    number = 1,
    pages = '35-49',
    year = 2003,
    doi = '10.1016/S0169-7439(02)00090-4',
  }
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
  # feemparafac methods should be able to use the environment and subset
  plot(sh$factors[[1]])
  plot(sh)
  plot(sh, 'factors')
  head(coef(sh))
  head(coef(sh, 'factors'))
}}
\keyword{hplot}
\keyword{method}
