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
    A numeric vector or a list containing the half-widths of the
    scattering bands, in nm. Rayleigh scattering is followed by Raman
    scattering, followed by second diffraction order for Rayleigh and
    Raman, and so on. (Typically, there's no need for anything higher
    than third order, and even that is rare.) For example:

    \enumerate{
      \item Rayleigh scattering
      \item Raman scattering
      \item Rayleigh scattering, \eqn3{2\lambda}{2 lambda}{2&lambda;}
      \item Raman scattering, \eqn3{2\lambda}{2 lambda}{2&lambda;}
      \item \ldots
    }

    For higher diffraction orders, the peak widths are proportionally
    scaled, making it possible to provide the same number for all kinds
    of scattering visible in the EEM. Set a width to \eqn{0} if you
    don't want to handle this particular kind of scattering signal.

    It's possible to specify the bands asymmetrically. If the area to be
    corrected should range from \var{x} nm to the left of the scattering
    peak to \var{y} nm to the right of it, pass a list instead of a
    vector, and put a two-element vector \code{c(\var{x}, \var{y})} for
    the appropriate kind of scattering. For example, passing
    \code{widths = list(c(30, 20), 20)} means \dQuote{handle \eqn{-30}
    nm to the left and \eqn{+20} nm to the right of Rayleigh peak and
    \eqn3{\pm 20}{+/- 20}{&pm; 20} nm around Raman peak}.

    To sum up, given two half-widths \eqn3{W_1}{W[1]}{W<sub>1</sub>} and
    \eqn3{W_2}{W[2]}{W<sub>2</sub>}, the test for being inside a
    \eqn{k}th diffraction order of a scattering band is as follows:

    \deqn3{
      -W_1 <
      \frac{\lambda_\mathrm{center}}{k} - \lambda_\mathrm{em}
      < +W_2
    }{-W[1] < lambda.center/k - lambda.em < +W[2]}{
      -W<sub>1</sub> &lt;
      <sup>&lambda;<sub>center</sub></sup>/<sub>k</sub>
        - &lambda;<sub>em</sub>
      &lt; +W<sub>2</sub>
    }
    \makefigure{
      pdf('man/figures/scatter-widths.pdf', 6, 3.2, pointsize = 10)
      dev.control(displaylist = 'enable')

      em.range <- c(300, 700)
      ex.range <- c(210, 550)
      widths <- list(
        list(c(50, 20), 20),
        list(c(20, 20), 20),
        list(c(10, 10), 0)
      )
      Raman.shift <- 3400

      # 1. The rectangle itself
      image(
        matrix(0, 1, 1), col = NA,
        xlim = em.range, ylim = ex.range,
        xlab = quote(lambda[em] * ', nm'),
        ylab = quote(lambda[ex] * ', nm'),
        main = 'widths = list(c(50, 20), 20, 20, 20, 10), Raman.shift = 3400'
      )

      # 2. The scattering areas and their centers
      sc <- c(min(em.range, ex.range), max(em.range, ex.range))
      Map(function(p, k) {
        polygon(
          c(sc * k - p[[1]][[1]] * k, rev(sc * k) + p[[1]][[2]] * k),
          c(sc, rev(sc)),
          border = NA, col = '#D0D0D0'
        )
        if (p[[2]] > 0) {
          wl.ex <- seq(sc[1], sc[2], 1)
          wl.em <- 1/(1/wl.ex - Raman.shift/1e7)
          polygon(
            c(wl.em * k - p[[2]] * k, rev(wl.em * k) + p[[2]] * k),
            c(wl.ex, rev(wl.ex)),
            border = NA, col = '#D0D0D0'
          )
          lines(wl.em * k, wl.ex, lty = 3)
        }
        lines(sc * k, sc, lty = 3)
      }, widths, seq_along(widths))

      # arrow length
      dp <- 15

      # Raman shift arrow
      p <- 400
      arrows(
        p, p, x1 = 1/(1/p - Raman.shift/1e7), len = .1
      )
      text(
        1/(1/p - Raman.shift/1e7) + dp / 2, p,
        bquote(
          (lambda[ex]^-1 - lambda[em]^-1) %.% 10^7
          == .(Raman.shift) * ' cm'^-1
        ),
        adj = c(0, .6)
      )

      # scatter width arrows
      p <- 530
      segments(p, p, x1 = p - widths[[1]][[1]][1])
      arrows(
        c(p - widths[[1]][[1]][1] - dp, p + dp), p,
        x1 = c(p - widths[[1]][[1]][1], p), len = .1
      )
      text(
        p - widths[[1]][[1]][1] - dp, p,
        paste('widths[[1]][1] = ', widths[[1]][[1]][1], 'nm'),
        pos = 2
      )

      p <- 470
      segments(p, p, x1 = p + widths[[1]][[1]][2])
      arrows(
        c(p + widths[[1]][[1]][2] + dp, p - dp), p,
        x1 = c(p + widths[[1]][[1]][2], p), len = .1
      )
      text(
        p - dp, p,
        paste('widths[[1]][2] =', widths[[1]][[1]][2], 'nm'),
        pos = 2
      )

      p <- 500
      pR <- 1/(1/p - Raman.shift/1e7)
      segments(pR, p, x1 = pR + widths[[1]][[2]])
      arrows(
        c(pR - dp, pR + widths[[1]][[2]] + dp), p,
        x1 = c(pR, pR  + widths[[1]][[2]]), len = .1
      )
      text(
        pR + widths[[1]][[2]], p,
        paste('widths[[2]] =', widths[[2]][[1]][1], 'nm'),
        pos = 1
      )

      p <- 310
      segments(2 * p, p, x1 = 2 * (p - widths[[2]][[1]][1]))
      arrows(
        c(2 * (p - widths[[2]][[1]][1]) - dp, 2 * p + dp), p,
        x1 = c(2 * (p - widths[[2]][[1]][1]), 2 * p), len = .1
      )
      text(
        2 * (p - widths[[2]][[1]][1]) - dp, p,
        bquote('widths[[3]]' / 2 == .(widths[[2]][[1]][1]) * ' nm'),
        pos = 2
      )

      p <- 230
      segments(3 * p, p, x1 = 3 * (p - widths[[2]][[1]][1]))
      arrows(
        c(3 * (p - widths[[3]][[1]][1]) - dp, 3 * p + dp), p,
        x1 = c(3 * (p - widths[[3]][[1]][1]), 3 * p), len = .1
      )
      text(
        3 * (p - widths[[2]][[1]][1]), p,
        bquote('widths[[5]]' / 3 == .(widths[[3]][[1]][1]) * ' nm'),
        pos = 2
      )

      dev.print(
        svg, 'man/figures/scatter-widths.svg',
        width = 6, height = 3.2, pointsize = 10
      )
      dev.off()
    }
    \if{latex}{\figure{scatter-widths.pdf}{options: width=4.5in}}
    \if{html}{ \figure{scatter-widths.svg}{options: width="700"}}

    \if{latex,html}{
      In this example, a much larger portion of the anti-Stokes area is
      removed near the first order Rayleigh scattering signal than in
      the Stokes area. This can be useful to get rid of undesired signal
      where no fluorescence is observed on some spectrometers. The
      second and third order scattering signal areas are automatically
      scaled \eqn{2} and \eqn{3} times, respectively.
    }
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
        implemented in \CRANpkg{pracma} function \code{\link{kriging}}.
        Pass a \code{type} argument to choose between the two
        methods. This method is not recommended due to its high
        CPU time and memory demands: it has to invert a dense
        \eqn3{O(N^2)}{O(N^2)}{O(N<sup>2</sup>)} matrix (which easily
        reaches multiple gigabytes for some EEMs), then take its products
        with vectors \eqn{O(N)} times, with \eqn{N =} \code{length(x)}.
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
  \CRANpkg{staRdom} package itself].

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

  The \dQuote{whittaker} method [4] works by unfolding each FEEM into a
  vector \eqn3{\mathbf z}{z}{<b>z</b>} and looking for a vector
  \eqn3{\hat{\mathbf z}}{z.hat}{<b>z&#x0302;</b>} that is close enough
  to \eqn3{\mathbf z}{z}{<b>z</b>}, but also smooth, by minimising a sum
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
