\name{feemindex}
\alias{feemindex}
\alias{feemindex.feem}
\alias{feemindex.feemcube}
\alias{feemindex.list}
\title{Fluorescence indices and peak values}
\description{
  Calculate fluorescence indices or peak values for individual FEEMs or
  groups of them.
}
\usage{
  feemindex(x, \dots)
  \method{feemindex}{feem}(
    x,
    indices = c(
      "HIX", "BIX", "MFI", "CFI", "YFI", "FrI",
      "A", "B", "C", "M", "P", "T"
    ),
    tolerance = 1, interpolate = FALSE, \dots
  )
  \method{feemindex}{feemcube}(x, \dots, progress = FALSE)
  \method{feemindex}{list}(x, \dots, progress = FALSE)
}
\arguments{
  \item{x}{
    A FEEM, a FEEM cube, or a list of \code{\link{feem}} objects.
  }
  \item{indices}{
    Fluorescence indices or peaks to return. By default, all indices and
    peaks known to the function are returned. See Details for their
    meaning.
  }
  \item{tolerance}{
    A numeric scalar signifying the acceptable emission and excitation
    wavelength error in nm. For example, if a wavelength of \eqn{254} nm is
    needed to calculate an index, a value at \eqn{255} nm can be considered if
    \code{tolerance >= 1}. Defaults to \eqn{1} nm. See below for what happens
    if no matching value is found.
  }
  \item{interpolate}{
    A string specifying an interpolation method (\dQuote{whittaker},
    \dQuote{loess}, \dQuote{kriging}, \dQuote{pchip}), or \code{FALSE}
    to disable interpolation (default).

    If interpolation is disabled, an index will get an \code{NA} value
    when required points are too far from the measured grid or are
    present in the grid but set to \code{NA}.

    When interpolation is enabled, required points that are missing from
    the grid or present but set to \code{NA} will be interpolated using
    \code{\link{feemgrid}} as long as they are within the wavelength
    bounds of the FEEM. \code{NA}s may still be returned only when the
    desired value is impossible to interpolate due to it being outside
    the wavelength range.
  }
  \item{\dots}{
    Additional parameters eventually passed to interpolation methods.
    See \code{\link{feemscatter}} for details.
  }
  \item{progress}{
    Set to \code{TRUE} to enable a progress bar (implemented via
    \code{\link[utils]{txtProgressBar}}).
  }
}
\details{
  Available indices and peaks are: \describe{
    \item{HIX}{
      \deqn3{
        \mathrm{HIX} = \frac{
          \int_{435 \, \mathrm{nm}}^{480 \, \mathrm{nm}}
          I \, d\lambda_\mathrm{em}
        }{
          \int_{300 \, \mathrm{nm}}^{345 \, \mathrm{nm}}
          I \, d\lambda_\mathrm{em}
        }
        \; \mathrm{at} \; \lambda_\mathrm{ex} = 254 \, \mathrm{nm}
      }{HIX = area(em in [435, 480])/area(em in [300, 345]) at ex = 254}{
        HIX = <sup>
          &int; I(&lambda;<sub>em</sub> &in; [435, 480] nm)
        </sup>/<sub>
          &int; I(&lambda;<sub>em</sub> &in; [300, 345] nm)
        </sub> at &lambda;<sub>ex</sub> = 254 nm
      }

      Higher values of the humification index correspond to more
      condensed fluorescing molecules (higher C/H), more humified
      matter. \cite{Zsolnay1999}
    }
    \item{BIX}{
      \deqn3{
        \mathrm{BIX} =
          \frac{
            I(\lambda_\mathrm{em} = 380 \, \mathrm{nm})
          }{
            I(\lambda_\mathrm{em} = 430 \, \mathrm{nm})
          } \; \mathrm{at} \; \lambda_\mathrm{ex} = 310 \, \mathrm{nm}
      }{BIX = I(em = 380 nm)/I(em = 430 nm) at ex = 310 nm}{
        BIX = <sup>
          I(&lambda;<sub>em</sub> = 380 nm)
        </sup>/<sub>
          I(&lambda;<sub>em</sub> = 430 nm)
        </sub> at &lambda;<sub>ex</sub> = 310 nm
      }

      Index of recent autochthonous contribution determines the presence
      of the \eqn{\beta} fluorophore, characteristic of autochthonous
      biological activity in water samples. \cite{Huguet2009}
    }
    \item{MFI}{
      \deqn3{
        \mathrm{MFI} =
          \frac{
            I(\lambda_\mathrm{em} = 450 \, \mathrm{nm})
          }{
            I(\lambda_\mathrm{em} = 500 \, \mathrm{nm})
          } \; \mathrm{at} \; \lambda_\mathrm{ex} = 370 \, \mathrm{nm}
      }{MFI = I(em = 450 nm)/I(em = 500 nm) at ex = 370 nm}{
        MFI = <sup>
          I(&lambda;<sub>em</sub> = 450 nm)
        </sup>/<sub>
          I(&lambda;<sub>em</sub> = 500 nm)
        </sub> at &lambda;<sub>ex</sub> = 370 nm
      }

      The fluorescence index by \cite{McKnight2001} helps distinguish
      sources of isolated aquatic fulvic acids and may indicate their
      aromaticity.
    }
    \item{CFI}{
      \deqn3{
        \mathrm{CFI} =
          \frac{
            I(\lambda_\mathrm{em} = 470 \, \mathrm{nm})
          }{
            I(\lambda_\mathrm{em} = 520 \, \mathrm{nm})
          } \; \mathrm{at} \; \lambda_\mathrm{ex} = 370 \, \mathrm{nm}
      }{CFI = I(em = 470 nm)/I(em = 520 nm) at ex = 370 nm}{
        CFI = <sup>
          I(&lambda;<sub>em</sub> = 470 nm)
        </sup>/<sub>
          I(&lambda;<sub>em</sub> = 520 nm)
        </sub> at &lambda;<sub>ex</sub> = 370 nm
      }

      The fluorescence index by \cite{Cory2005} is correlated to
      relative contribution of microbial versus higher plant-derived
      organic matter to the DOM pool.
    }
    \item{YFI}{
      \deqn3{
        \mathrm{YFI} = \frac{
          \bar{I}(\lambda_\mathrm{em} \in [350, 400] \, \mathrm{nm})
        }{
          \bar{I}(\lambda_\mathrm{em} \in [400, 450] \, \mathrm{nm})
        } \; \mathrm{at} \; \lambda_\mathrm{ex} = 280 \, \mathrm{nm}
      }{YFI = mean(em in [350, 400])/mean(em in [400, 450]) at ex = 280}{
        YFI = <sup>
          <span style="font-style:normal;">mean</span>
          I(&lambda;<sub>em</sub> &in; [350, 400] nm)
        </sup>/<sub>
          <span style="font-style:normal;">mean</span>
          I(&lambda;<sub>em</sub> &in; [400, 450] nm)
        </sub> at &lambda;<sub>ex</sub> = 280 nm
      }

      Yeomin fluorescence index \cite{Heo2016} is lowest for humic-like
      and fulvic-like samples, higher for aminosugar-like samples and
      highest for protein-like samples.
    }
    \item{FrI}{
      \deqn3{
        \mathrm{FrI} = \frac{
          I(\lambda_\mathrm{em} = 380 \, \mathrm{nm})
        }{
          \max I(\lambda_\mathrm{em} \in [420, 435] \, \mathrm{nm})
        } \; \mathrm{at} \; \lambda_\mathrm{ex} = 310 \, \mathrm{nm}
      }{FrI = I(em = 380 nm)/max(I, em in [420, 435] nm) at ex = 310 nm}{
        FrI = <sup>
          I(&lambda;<sub>em</sub> = 380 nm)
        </sup>/<sub>
          <span style="font-style:normal;">max</span>
          I(&lambda;<sub>em</sub> &in; [420, 435] nm)
        </sub> at &lambda;<sub>ex</sub> = 310 nm
      }

      The freshness index, also known as
      \eqn3{\frac{\beta}{\alpha}}{\beta/\alpha}{
        <sup>&beta;</sup>/<sub>&alpha;</sub>}, is an indicator of
      autochthonous inputs \cite{Wilson2009} and may provide indication
      of relative contribution of microbially produced DOM.
    }
    \item{A, B, C, M, P, T}{
      Fluorophore peaks taken from \cite{Coble2007}: \tabular{cccl}{
        Peak
        \tab \eqn3{\lambda_\mathrm{ex}}{lambda.ex}{&lambda;<sub>ex</sub>}
        \tab \eqn3{\lambda_\mathrm{em}}{lambda.em}{&lambda;<sub>em</sub>}
        \tab Fluorescence \cr

        A \tab   260   \tab 400-460 \tab humic-like        \cr
        B \tab   275   \tab   305   \tab tyrosine-like     \cr
        C \tab 320-360 \tab 420-460 \tab humic-like        \cr
        M \tab 290-310 \tab 370-410 \tab marine humic-like \cr
        P \tab   398   \tab   660   \tab pigment-like      \cr
        T \tab   275   \tab   340   \tab tryptophan-like
      }
      When a range of wavelengths specified in one or both axes, the
      maximal signal value over that range is taken.
    }
  }

  \if{latex,html}{
    \makefigure{
      indices <- list(
        HIX = list(254, list(c(435, 480), c(300, 345))),
        BIX = list(310, list(380, 430)),
        MFI = list(370, list(450, 500)),
        CFI = list(370, list(470, 520)),
        YFI = list(280, list(c(350, 400), c(400, 450))),
        FrI = list(310, list(380, c(420, 435)))
      )
      peaks <- list(
        'A' = list(260, c(400, 460)),
        'B' = list(275, 305),
        'C' = list(c(320, 360), c(420, 460)),
        'M' = list(c(290, 310), c(370, 410)),
        'P' = list(398, 660),
        'T' = list(275, 340)
      )

      expand <- function(l) lapply(l, function(pair) {
        if (is.list(pair[[2]])) lapply(
          pair[[2]], function(em) list(pair[[1]], em)
        )
        else list(pair)
      })

      draw.areas <- function(
        x, ex.range = c(240, 410), em.range = c(260, 675), lwd = 1.5,
        col = c("#DF536B", "#61D04F", "#2297E6", "#28E2E5", "#E5C010", "#CD0BBC"),
        pch = c(0, 19, 2:5), sc.delta = 20, sc.Raman = 3400
      ) {
        image(
          matrix(0, 1, 1), col = NA,
          xlim = em.range, ylim = ex.range,
          xlab = quote(lambda[em] * ', nm'),
          ylab = quote(lambda[ex] * ', nm')
        )

        # draw scattering regions
        sc <- c(min(em.range, ex.range), max(em.range, ex.range))
        for (mul in 1:2) {
          polygon(
            c(mul * sc - sc.delta, rev(mul * sc) + sc.delta), c(sc, rev(sc)),
            border = NA, col = '#D0D0D0'
          )
          wl.ex <- seq(sc[1], sc[2], 1)
          wl.em <- 1/(1/wl.ex - sc.Raman/1e7)
          polygon(
            c(mul * wl.em - sc.delta, rev(mul * wl.em) + sc.delta),
            c(wl.ex, rev(wl.ex)),
            border = NA, col = '#D0D0D0'
          )
        }

        x <- expand(x)
        for (i in seq_along(x)) {
          for (p in x[[i]]) {
            if (sum(lengths(p)) == 2) {
              if (length(x[[i]]) > 1)
                points(p[[2]], p[[1]], pch = pch[i], col = col[i], lwd = lwd)
            } else if (sum(lengths(p)) == 3) {
              lines(
                rep_len(p[[2]], 2), rep_len(p[[1]], 2),
                col = col[i], lwd = lwd
              )
            } else if (sum(lengths(p)) == 4) {
              rect(
                p[[2]][1], p[[1]][1], p[[2]][2], p[[1]][2],
                border = col[i], lwd = lwd
              )
            }
          }
          if (length(x[[i]]) == 1) {
            xx <- x[[i]][[1]]
            text(
              mean(xx[[2]]), mean(xx[[1]]),
              names(x)[i], col = col[i], cex = 1.2, offset = .3,
              pos = if (sum(lengths(xx)) == 3) 3 else NULL
            )
          }
        }

        if (all(lengths(x) == 2)) {
          kinds <- lapply(expand(indices), sapply, function(p) sum(lengths(p)))
          lty. <- rep_len(1, length(kinds))
          lty.[sapply(kinds, function(x) all(x != 3))] <- NA
          pch. <- pch
          pch.[sapply(kinds, function(x) all(x != 2))] <- NA

          legend(
            'topright', legend = names(x),
            col = col, cex = .9,
            lty = lty., pch = pch., lwd = lwd
          )
        }
      }

      pdf('man/figures/feemindex.pdf', 6, 2.8, pointsize = 9)
      dev.control(displaylist = 'enable')

      par(mfrow = c(1, 2), mar = c(5, 4, 2, 2) + .1)
      draw.areas(indices)
      draw.areas(peaks)

      dev.print(
        svg, 'man/figures/feemindex.svg',
        width = 6, height = 2.8, pointsize = 9
      )
      dev.off()
    }

    \if{latex}{\figure{feemindex.pdf}{options: width=6in}}
    \if{html}{ \figure{feemindex.svg}{options: width="700"}}

    Positions of the peaks and the areas used to determine the
    fluorescence indices of an EEM. The Rayleigh and Raman scattering
    areas for both 1st and 2nd diffraction orders are shown in grey,
    assuming a width of \eqn3{\pm 20}{+/- 20}{&pm;20} nm and a Raman
    shift of \eqn3{3400 \: \mathrm{cm}^{-1}}{3400 cm^-1}{3400
      <span style="font-style:normal;">cm</span><sup>-1</sup>}. The
    tolerance interval of \eqn3{\pm 1}{+/- 1}{&pm;1} nm is invisible
    at the scale of the figure.

    Integration for HIX and YFI is done using the trapezoidal method:

    \deqn3{
      \int_a^b f(x) dx \approx (b - a) \frac{f(a) + f(b)}{2}
    }{area(f, x in [a, b]) ~ (b - a) (f(a) + f(b) / 2}{
      &int; f(x &in; [a, b]) &approx; (b - a)
        <sup>f(a) + f(b)</sup>/<sub>2</sub>
    }
  }

}
\value{
  For individual \code{\link{feem}} objects, a named numeric vector
  containing the values requested via the \code{indices} argument.

  Otherwise, a \code{\link[base]{data.frame}} containing the values from
  the vectors above and a column named \code{sample} containing the
  names of the samples (or numbers, if names were absent).
}
\author{
  With edits and suggestions by Anastasia Drozdova.
}
\references{\bibliography}

\seealso{\code{\link{feem}}}
\examples{
  data(feems)

  x <- feemscatter(feems$a, rep(25, 4), 'omit')
  feemindex(x)
  feemindex(x, interpolate = 'whittaker')

  feemindex(feems[2:3])
  feemindex(feemcube(feems[4:5], TRUE))
}
\keyword{methods}
