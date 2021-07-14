\encoding{UTF-8}
\name{feemscatter}
\alias{feemscatter}
\alias{feemscatter.feem}
\alias{feemscatter.list}
\alias{feemscatter.feemcube}
\title{Handle scattering signal in FEEMs}
\description{
  Remove or interpolate scattering signal in individual FEEM objects,
  FEEM cube objects, or lists of them.
}
\usage{
  feemscatter(x, \dots)
  \method{feemscatter}{list}(x, \dots, cl, progress = TRUE)
  \method{feemscatter}{feemcube}(x, \dots, cl, progress = TRUE)
  \method{feemscatter}{feem}(
    x, widths, method = c("omit", "pchip", "loess", "kriging", "whittaker"),
    add.zeroes = 30, Raman.shift = 3400, \dots
  )
}
\arguments{
  \item{x}{
    An individual FEEM object, FEEM cube object, or a list of them, to
    handle the scattering signal in.
  }
  \item{widths}{
    A numeric vector of length 4 containing the widths (in nm) of the
    scattering signal, in the following order:
    \enumerate{
      \item Rayleigh scattering
      \item Raman scattering
      \item Rayleigh scattering, \eqn3{2\lambda}{2 lambda}{2&lambda;}
      \item Raman scattering, \eqn3{2\lambda}{2 lambda}{2&lambda;}
    }
    Set a width to \eqn{0} if you don't want to handle this particular
    kind of scattering signal.
  }
  \item{method}{
    A string choosing \emph{how} to handle the scattering signal:
    \describe{
      \item{omit}{Replace it with \code{NA}.}
      \item{pchip}{
        Interpolate it line-by-line using piecewise cubic Hermitean
        polynomials (\code{\link{pchip}}). Pass a \code{by} argument
        to choose the direction of interpolation; see Details.
     }
      \item{loess}{
        Interpolate it by fitting a locally weighted polynomial surface
        (\code{\link{loess}}). In this case the remaining parameters are
        passed verbatim to \code{\link{loess}}, which may be used to set
        parameters such as \code{span}.
      }
      \item{kriging}{
        Interpolate it by means of ordinary or simple Kriging, as
        implemented in \pkg{pracma} function \code{\link{kriging}}.
        Pass a \code{type} argument to choose between the two
        methods. This method is not recommended due to its high
        CPU time and memory demands: it has to invert a dense
        \eqn3{O(N^2)}{O(N^2)}{O(N<sup>2</sup>)} matrix (which easily
        reaches multiple gigabytes for some EEMs), then take its products
        with vectors \eqn{O(N)} times.
      }
      \item{whittaker}{
        Interpolate it by minimising a weighted sum of squared residuals
        (for known part of the spectrum) and roughness penalty (squared
        central difference approximations for derivatives by
        \eqn3{\lambda_\mathrm{em}}{\lambda.em}{&lambda;<sub>em</sub>}
        and \eqn3{\lambda_\mathrm{ex}}{\lambda.ex}{&lambda;<sub>ex</sub>}).
        See Details for more information and parameters.
      }
    }
  }
  \item{add.zeroes}{
    Set intensities at \eqn3{
      \lambda_\mathrm{em} < \lambda_\mathrm{ex}
      - \mathtt{add.zeroes}\:\mathrm{nm}
    }{\lambda.em < \lambda.ex - add.zeroes (nm)}{
      &lambda;<sub>em</sub> &lt;
      &lambda;<sub>ex</sub> - <code>add.zeroes</code>
      <span style="font-style: normal;">nm</span>
    } to \eqn{0} unless they have been measured. Set to \code{NA} to
    disable this behaviour.
  }
  \item{Raman.shift}{
    Raman shift of the scattering signal of water,
    \eqn3{\textrm{cm}^{-1}}{cm^-1}{cm<sup>-1</sup>}.
  }
  \item{\dots}{
    Passed verbatim from \code{feemscatter} generics to
    \code{feemscatter.feem}.

    If \dQuote{pchip} method is selected, the \code{by} parameter
    chooses between interpolating by row, by column, or averaging both,
    see Details.

    If \dQuote{loess} method is selected, remaining options are passed
    to \code{\link{loess}} (the \code{span} parameter is of particular
    interest there).

    If \dQuote{kriging} method is selected, remaining
    options are passed to \code{\link{kriging}}.

    If \dQuote{whittaker} method is selected, available parameters
    include \code{d}, \code{lambda}, \code{nonneg} and \code{logscale},
    see Details.
  }
  \item{cl}{
    If not \code{\link{missing}}, a \pkg{parallel} cluster object to
    run the scattering correction code on or \code{NULL} for the default
    cluster object registered via \code{\link{setDefaultCluster}}.
  }
  \item{progress}{
    Set to \code{FALSE} to disable the progress bar.
  }
}
\details{
  The \dQuote{pchip} method works by default as described in [1]: each
  emission spectrum at different excitation wavelengths is considered one
  by one. Zeroes are inserted in the corners of the spectrum if they are
  undefined (\code{NA}) to prevent extrapolation from blowing up, then
  the margins are interpolated using the corner points, then the rest of
  the spectrum is interpolated line by line. Since \code{\link{pchip}}
  requires at least 3 points to interpolate, the function falls back to
  linear interpolation if it has only two defined points to work with.
  The \code{by} argument controls whether the function proceeds by rows
  of the matrix (\dQuote{emission}, default), by columns of the matrix
  (\dQuote{excitation}), or does both (\dQuote{both}) and averages the
  results to make the resulting artefacts less severe [2, see the
  \pkg{staRdom} package itself].

  The \dQuote{loess} method feeds the whole FEEM except the area to
  be interpolated to \code{\link{loess}}, then asks it to predict the
  remaining part of the spectrum. Any negative values predicted by
  \code{\link{loess}} are replaced by \eqn{0}.

  The \dQuote{kriging} method [3] is much more computationally expensive
  than the previous two, but, on some spectra, provides best results,
  not affected by artefacts resulting from line-by-line one-dimensional
  interpolation (\code{pchip}) or varying degrees of smoothness in
  different areas of the spectrum (\code{loess}). Any negative values
  returned by \code{\link{kriging}} are replaced by \eqn{0}.

  The \dQuote{whittaker} method [4] works by unfolding \code{x}
  into a vector \eqn3{\mathbf z}{z}{<b>z</b>} and looking for a vector
  \eqn3{\hat{\mathbf z}}{z.hat}{<b>z&#x0302;</b>} that is close enough to
  \eqn3{\mathbf z}{z}{<b>z</b>}, but also smooth, by minimising a sum
  of penalties:

  \deqn3{
    (\mathbf z - \hat{\mathbf z})^\top
    \mathrm{diag}(\mathbf w)
    (\mathbf z - \hat{\mathbf z})
    + \sum_k \lambda_k |\mathbf{D}_{n_k} \hat{\mathbf z}|^2
  }{sum(w * (z.hat - z)^2) + sum over k (\lambda[k] * D(n[k]) \\\%*\\\% z.hat)}{
    (<b>z</b> - <b>z&#x0302;</b>)<sup>
      <span style="font-style: normal;">T</span>
    </sup>
    <span style="font-style: normal;">diag</span>(<b>w</b>)
    (<b>z</b> - <b>z&#x0302;</b>)
    + &#x2211;<sub>k</sub> &lambda;<sub>k</sub>
    |<b>D</b><sub>n<sub>k</sub></sub> <b>z&#x0302;</b>|<sup>2</sup>
  }

  The weights \eqn3{\mathbf w}{w}{<b>w</b>} are set to \eqn{0} for
  missing (\code{NA}) points and for those to be interpolated and to
  \eqn{1} otherwise. The matrix
  \eqn3{\mathbf D_n}{D(n)}{<b>D</b><sub>n</sub>} is constructed in
  such a way that multiplying it by
  \eqn3{\hat{\mathbf z}}{z.hat}{<b>z&#x0302;</b>}
  results in a vector of \eqn{n}-th order derivative estimates in both
  directions and in all applicable points of
  \eqn3{\hat{\mathbf z}}{z.hat}{<b>z&#x0302;</b>} as a matrix. The
  wavelength grid is correctly taken into account by solving a
  Vandermonde system for every \eqn{n+1} consecutive points.

  The parameters \code{d} and \code{lambda} should be numeric vectors
  of the same length, choosing the difference orders
  (\eqn3{n_k}{n[k]}{n<sub>k</sub>}) and their weights
  (\eqn3{\lambda_k}{\lambda[k]}{&lambda;<sub>k</sub>}). It has been
  shown in [5] that a combination of first- and second-order penalty
  (\eqn3{
    2 \lambda \mathbf{D}_1 + \lambda^2 \mathbf{D}_2
  }{2 \lambda D(1) + \lambda^2 D(2)}{2&lambda;<b>D</b><sub>1</sub> +
    &lambda;<sup>2</sup><b>D</b><sub>2</sub>}) results in non-negative
  impulse response, but the resulting peak shape may be
  sub-optimal. Instead, the default penalty is
  \eqn3{
    10^{-2} \mathbf{D}_1 + 10 \mathbf{D}_2
  }{10^-2 D(1) + 10 D(2)}{10<sup>-2</sup> <b>D</b><sub>1</sub> +
    10 <b>D</b><sub>2</sub>
  }, and resulting negative values are pulled to \eqn{0} with weight
  \code{nonneg} (default \eqn{1}, same as fidelity weight) by adding a
  penalty of \eqn3{
    \mathtt{nonneg} \cdot \sum_i \mathbf{1}_{\hat{z}_i < 0} \, \hat{z}_i^2
  }{nonneg * sum(z.hat[z.hat < 0]^2)}{
    <code>nonneg</code> &middot; &#x2211;<sub>i</sub>
    <b>1</b><sub>z&#x0302;<sub>i</sub> &lt; 0</sub>
    z&#x0302;<sub>i</sub><sup>2</sup>
  } and retrying until no new penalty weights are added. Set \code{nonneg}
  to \eqn{0} to disable this behaviour.

  It is also possible to deal with resulting negative values by
  scaling and shifting the signal between \code{logscale} (typically
  \eqn3{10^{-4}}{10^-4}{10<sup>-4</sup>}) and \eqn{1}, interpolating
  the logarithm of the signal, then undoing the transformation. This
  prevents the resulting values from getting lower than \eqn3{
    \mathrm{min}(x) - (\mathrm{max}(x) - \mathrm{min}(x))
    \frac{\mathtt{logscale}}{1 - \mathtt{logscale}}
  }{min(x) - diff(range(x)) * logscale / (1 - logscale)}{
    <span style="font-style: normal;">min</span>(x) -
    (<span style="font-style: normal;">max</span>(x) -
      <span style="font-style: normal;">min</span>(x))
    <sup><code>logscale</code></sup>/<sub>1 - <code>logscale</code></sub>
  }, which is approximately \eqn3{
    -\mathtt{logscale} \cdot \mathrm{max}(x)
  }{-logscale * max(x)}{
    -<code>logscale</code> &middot;
    <span style="font-style: normal;">max</span>(x)
  } if \code{logscale} and \eqn3{\mathrm{min}(x)}{min(x)}{
    <span style="font-style: normal;">min</span>(x)
  } are both close to \eqn{0}. By default \code{logscale} is \code{NA},
  disabling this behaviour, since it may negatively affect the shape
  of interpolated signal.
}
\value{
  An object of the same kind (FEEM object / FEEM cube / list of them)
  with scattering signal handled as requested.
}
\references{\enumerate{
  \item \bibentry{'Article',
    author = c(
      person('Morteza', 'Bahram'),
      person('Rasmus', 'Bro'),
      person('Colin', 'Stedmon'),
      person('Abbas', 'Afkhami')
    ),
    title = paste(
      'Handling of Rayleigh and Raman scatter for PARAFAC modeling of',
      'fluorescence data using interpolation'
    ),
    journal = 'Journal of Chemometrics',
    volume = 20,
    number = '3-4',
    pages = '99-105',
    year = 2006,
    doi = '10.1002/cem.978',
  }

  \item \bibentry{'Article',
    author = c(
      person('Matthias', 'Pucher'),
      person('Urban',  paste0(rawToChar(as.raw(0x5c)), 'enc{Wünsch}{Wuensch}')),
      person('Gabriele', 'Weigelhofer'),
      person('Kathleen', 'Murphy'),
      person('Thomas', 'Hein'),
      person('Daniel', 'Graeber')
    ),
    title = paste(
      'staRdom: Versatile Software for Analyzing Spectroscopic Data',
      'of Dissolved Organic Matter in R'
    ),
    journal = 'Water',
    volume = 11,
    number = 11,
    year = 2019,
    pages = 2366,
    doi = '10.3390/w11112366'
  }

  \item \bibentry{'InBook',
    author = c(
      person(c('William', 'H.'), 'Press'),
      person(c('Saul', 'A.'), 'Teukolsky'),
      person(c('William', 'T.'), 'Vetterling'),
      person(c('Brian', 'P.'), 'Flannery')
    ),
    booktitle =
      'Numerical recipes: The Art of Scientific Computing (3rd Ed.)',
    publisher = 'Cambridge University Press, New York',
    chapter = '3.7.4',
    title = 'Interpolation by Kriging',
    pages = '144-147',
    year = 2007
  }

  \item \bibentry{'Article',
    author = person(c('Paul', 'H.', 'C.'), 'Eilers'),
    title = 'A Perfect Smoother',
    journal = 'Analytical Chemistry',
    volume = 75,
    number = 14,
    pages = '3631-3636',
    year = 2003,
    doi = '10.1021/ac034173t'
  }

  \item \bibentry{'Article',
    author = c(
      person(c('Paul', 'H.', 'C.'), 'Eilers'),
      person(c('Jelle', 'J.'), 'Goeman')
    ),
    title = 'Enhancing scatterplots with smoothed densities',
    journal = 'Bioinformatics',
    volume = 20,
    number = 5,
    pages = '623-628',
    year = 2004,
    doi = '10.1093/bioinformatics/btg454'
  }
}}

\seealso{
  \code{\link{feem}}, \code{\link{feemcube}}
}
\examples{
  data(feems)
  plot(x <- feemscatter(
    feems[[1]], widths = c(25, 25, 20, 20),
    method = 'whittaker', Raman.shift = 3500
  ))
}
\keyword{method}
