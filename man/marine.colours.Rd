\name{marine.colours}
\alias{marine.colours}
\alias{diverging.colours}
\title{Perceptually uniform palettes}
\description{
  Create perceptually continuous palettes of \R colours.
}
\usage{
  marine.colours(
    n, chroma = 0.65, luminance = c(0.35, 1),
    alpha = 1, gamma = 1, fixup = TRUE
  )
  diverging.colours(
    n, chroma = c(.1, .75), luminance = c(1, .35),
    alpha = 1, gamma = 1, fixup = TRUE
  )
}
\arguments{
  \item{n}{
    Number of colours to return.
  }
  \item{chroma}{
    Specifies the chroma (how saturated should the colours be) for the
    palette, a real number between \eqn{0} and \eqn{1}. May also be a
    two-element vector, in which case the chroma is changed smoothly
    from start to finish of the resulting palette.
  }
  \item{luminance}{
    Specifies the luminance (how bright should the colours be) of the
    colours constituting the palette. Typically, a two-element vector
    of real numbers between \eqn{0} and \eqn{1} to indicate smooth
    change along the palette, but can also be a fixed number.
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
\details{
  The \code{marine.colours} palette is used by default by all
  \code{plot} methods (e.g. \code{\link{plot.feem}}) for FEEM-like data
  to show absolute values.  It is designed to retain perceptual
  uniformity even after complete desaturation.

  The \code{diverging.colours} palette is used by
  \code{\link{plot.feem.resid}} to display residual values. People with
  severe colour vision deficiency (tritanopia or monochromacy) won't be
  able to discern positive and negative branches of the palette, but
  it's supposed to be legible for people with deuteranopia and
  protanopia.
}
\value{
  A character vector of length \code{n} containing colour specifications
  for use with \R graphics functions.

  \if{latex,html}{
    \makefigure{
      png('man/figures/marinecolours.png', 400, 32)
      par(mar = rep(0, 4))
      n <- 256
      image(matrix(seq(n), n), col = albatross::marine.colours(n))
      dev.off()
    }
    \if{latex}{\figure{marinecolours.png}{options: width=2.5in}}
    \if{html}{ \figure{marinecolours.png}{options: width="400"}}

    The \code{marine.colours} palette at the default values of
    \eqn3{C_{uv}^* = 0.65}{chroma = 0.65}{C<sub>uv</sub><sup>*</sup> =
    0.65}, \eqn3{L^* \in [0.35; 1]}{luminance in [0.35, 1]}{L<sup>*</sup>
    &in; [0.35; 1]}, \eqn3{\alpha = \gamma = 1}{\alpha = \gamma = 1}{
    &alpha; = &gamma; = 1}.

    \makefigure{
      png('man/figures/divergingcolours.png', 400, 32)
      par(mar = rep(0, 4))
      n <- 256
      image(matrix(seq(n), n), col = albatross::diverging.colours(n))
      dev.off()
    }
    \if{latex}{\figure{divergingcolours.png}{options: width=2.5in}}
    \if{html}{ \figure{divergingcolours.png}{options: width="400"}}

    The \code{diverging.colours} palette at the default values of
    \eqn3{C_{uv}^* \in [0.1; 0.75]}{chroma in [0.1, 0.75]}{
    C<sub>uv</sub><sup>*</sup> &in; [0.1; 0.75]}, \eqn3{L^* \in [0.35;
    1]}{luminance in [0.35, 1]}{L<sup>*</sup> &in; [0.35; 1]},
    \eqn3{\alpha = \gamma = 1}{\alpha = \gamma = 1}{&alpha; = &gamma; = 1}.
  }
}
\references{
  Insired by cmocean palette called \dQuote{haline}
  (\url{https://matplotlib.org/cmocean/#haline}), but using \R's
  implementation of polar CIE-LUV colour space instead of CAM02-UCS.

  CUBEHELIX (\url{https://www.mrao.cam.ac.uk/~dag/CUBEHELIX/}) is a
  similar technique using BT.601 luminance coefficients and RGB colour
  space.
}
\seealso{\code{\link{plot.feem}}, \code{\link[grDevices]{hcl}}}
\examples{
  image(volcano, col = marine.colours(256))
}
\keyword{color}
