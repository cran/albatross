\name{marine.colours}
\alias{marine.colours}
\title{Marine colours}
\description{
  Create a perceptually contiguous palette of R colours, using hues
  typically associated with natural waters.
}
\usage{
marine.colours(
  n, chroma = 0.65, luminance = c(0.35, 1),
  alpha = 1, gamma = 1, fixup = TRUE
)
}
\arguments{
  \item{n}{
    Number of colours to return.
  }
  \item{chroma}{
    Specifies the chroma (how saturated should the colours be) for the
    palette, a real number between 0 and 1. May also be a two-element
    vector, in which case the chroma is changed smoothly from start to
    finish of the resulting palette.
  }
  \item{luminance}{
    Specifies the luminance (how bright should the colours be) of the
    colours constituting the palette. Typically, a two-element vector
    of real numbers between 0 and 1 to indicate smooth change along the
    palette, but can also be a fixed number.
  }
  \item{alpha}{
    Specifies the transparency of the colours of the palette. As above,
    can be a fixed number or a two-element vector in the range
    \eqn{[0,1]}. Typically, fully opaque (\code{alpha=1}) colours are
    used.
  }
  \item{gamma}{
    Provides the power coefficient for the hue/chroma/luminance/alpha
    growth formulae. May be useful when it is needed to sacrifice the
    perceptual linearity of the palette to provide more contrast for
    smaller or bigger values on the plot. The gamma-corrected values
    are obtained by computing
    \eqn3{x ^ \gamma}{x ^ \gamma}{x<sup>&gamma;</sup>},
    \eqn3{x \in [0;1]}{x in [0;1]}{x&nbsp;&isin;&nbsp;[0;1]}, then scaling
    the result linearly to the required range. Typically, linear growth
    (\code{gamma = 1}) is preferred.
  }
  \item{fixup}{
    Whether to correct the palette if the resulting colours happen to
    fall outside the valid RGB range (passed as-is to \code{hcl}).
    Unrepresentable colours are returned as \code{NA}s, but fixing the
    palette may make it less perceptually uniform.
  }
}
\value{
  A character vector of length \code{n} containing colour specifications
  for use with R graphics functions.
  \if{latex,html}{
    \makefigure{
      png('man/figures/marinecolours.png', 400, 32)
      par(mar = rep(0, 4))
      n <- 256
      image(matrix(seq(n), n), col = albatross::marine.colours(n))
      dev.off()
    }
    \center{
      \if{latex}{\figure{marinecolours.png}{options: width=2.5in}}
      \if{html}{ \figure{marinecolours.png}{options: width="400"}}
    }

    \center{
      The palette at the default values of \eqn3{C_{uv}^* = 0.65
      }{chroma = 0.65}{C<sub>uv</sub><sup>*</sup> = 0.65},
      \eqn3{L^* \in [0.35; 1]}{luminance in [0.35, 1]}{L<sup>*</sup> &in;
      [0.35; 1]}, \eqn3{\alpha = \gamma = 1}{\alpha = \gamma = 1}{
      &alpha; = &gamma; = 1}.
    }
  }
}
\references{
  Insired by cmocean palette called \dQuote{haline}
  (\url{https://matplotlib.org/cmocean/#haline}), but using R's
  implementation of polar CIE-LUV colour space instead of CAM02-UCS.
}

\seealso{
 \code{\link{hcl}} for the colour space used, CUBEHELIX
 (\url{https://www.mrao.cam.ac.uk/~dag/CUBEHELIX/}) for a
 similar technique using BT.601 luminance coefficients and RGB colour
 space.
}

\examples{
  image(volcano, col = marine.colours(256))
}
\keyword{color}
