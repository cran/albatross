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
    corrected should range from \v{x} nm to the left of the scattering
    peak to \v{y} nm to the right of it, pass a list instead of a
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
          (lambda[ex]^-1 - lambda[em]^-1) \\\%.\\\% 10^7
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
        polynomials (\code{\link[pracma]{pchip}}). Pass a \code{by}
        argument to choose the direction of interpolation; see Details.
     }
      \item{loess}{
        Interpolate it by fitting a locally weighted polynomial surface
        (\code{\link[stats]{loess}}). Extra arguments are passed
        verbatim to \code{\link[stats]{loess}}, which may be used to set
        parameters such as \code{span}.
      }
      \item{kriging}{
        Interpolate it by means of ordinary or simple Kriging, as
        implemented in \CRANpkg{pracma} function
        \code{\link[pracma]{kriging}}.  Pass a \code{type} argument to
        choose between the two methods. This method is not recommended
        due to its high CPU time and memory demands: it has to invert a
        dense \eqn3{O(N^2)}{O(N^2)}{O(N<sup>2</sup>)} matrix (which
        easily reaches multiple gigabytes for some EEMs), and compute
        its product with a vector then take scalar products \eqn{O(N)}
        times, with \eqn{N =} \code{length(x)}.
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
    } to \eqn{0} unless they have been measured in order to stabilise
    the resulting decomposition \cite{Thygesen2004}. Set to \code{NA} to
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
    to \code{\link[stats]{loess}} (the \code{span} parameter is of
    particular interest there).

    If \dQuote{kriging} method is selected, remaining options are passed
    to \code{\link[pracma]{kriging}}.

    If \dQuote{whittaker} method is selected, available parameters
    include \code{d}, \code{lambda}, \code{nonneg} and \code{logscale},
    see Details.
  }
  \item{cl}{
    If not \code{\link[base]{missing}}, a \pkg{parallel} cluster object
    to run the scattering correction code on or \code{NULL} for the
    default cluster object registered via
    \code{\link[parallel]{setDefaultCluster}}.
  }
  \item{progress}{
    Set to \code{FALSE} to disable the progress bar.
  }
}
\details{
  The \dQuote{pchip} method works by default as described in
  \cite{Bahram2006}: each emission spectrum at different excitation
  wavelengths is considered one by one. Zeroes are inserted in the
  corners of the spectrum if they are undefined (\code{NA}) to prevent
  extrapolation from blowing up, then the margins are interpolated using
  the corner points, then the rest of the spectrum is interpolated line
  by line. Since \code{\link[pracma]{pchip}} requires at least 3 points
  to interpolate, the function falls back to linear interpolation if it
  has only two defined points to work with. The \code{by} argument
  controls whether the function proceeds by rows of the matrix
  (\dQuote{emission}, default), by columns of the matrix
  (\dQuote{excitation}), or does both (\dQuote{both}) and averages the
  results to make the resulting artefacts less severe \cite{Pucher2019}
  (see the \CRANpkg{staRdom} package itself).

  The \dQuote{loess} method feeds the whole FEEM except the area to be
  interpolated to \code{\link[stats]{loess}}, then asks it to predict
  the remaining part of the spectrum. Any negative values predicted by
  \code{\link[stats]{loess}} are replaced by \eqn{0}.

  The \dQuote{kriging} method \cite{NR3-Kriging} is much more
  computationally expensive than the previous two, but, on some spectra,
  provides best results, not affected by artefacts resulting from
  line-by-line one-dimensional interpolation (\code{pchip}) or varying
  degrees of smoothness in different areas of the spectrum
  (\code{loess}). Any negative values returned by
  \code{\link[pracma]{kriging}} are replaced by \eqn{0}.

  \subsection{Whittaker smoothing}{
    The \dQuote{whittaker} method \cite{Krylov2023} works by minimising a
    sum of penalties, requiring the interpolated surface to be close to
    the original points around it and to be smooth in terms of derivatives
    by
    \eqn3{\lambda_\mathrm{em}}{emission wavelength}{&lambda;<sub>em</sub>}
    and
    \eqn3{\lambda_\mathrm{ex}}{excitation wavelength}{&lambda;<sub>ex</sub>}.

    The parameters \code{d} and \code{lambda} should be numeric vectors
    of the same length, corresponding to the derivative orders (whole
    numbers \eqn{\ge 1}) and their respective penalty weights (small
    real numbers; larger is smoother). For interpolation purposes, the
    default penalty is \eqn3{
      10^{-2} \mathbf{D}_1 + 10 \mathbf{D}_2
    }{10^-2 D(1) + 10 D(2)}{10<sup>-2</sup> <b>D</b><sub>1</sub> +
      10 <b>D</b><sub>2</sub>},
    which corresponds to \code{d = 1:2} and \code{lambda = c(1e-2, 10)}.

    Any resulting negative values are pulled towards \eqn{0} by adding
    zero-valued points with weight \code{nonneg} (default \eqn{1}) and
    retrying. Set \code{nonneg} to \eqn{0} to disable this behaviour. It
    is also possible to deal with resulting negative values by scaling and
    shifting the signal between \code{logscale} (typically) and \eqn{1},
    interpolating the logarithm of the signal, then undoing the
    transformation. By default \code{logscale} is \code{NA}, disabling
    this behaviour, since it may negatively affect the shape of
    interpolated signal.

    See the internal help page \link{whittaker2} for implementation
    details.
  }
}
\value{
  An object of the same kind (FEEM object / FEEM cube / list of them)
  with scattering signal handled as requested.
}
\references{\bibliography}
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
