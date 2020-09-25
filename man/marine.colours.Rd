\name{marine.colours}
\alias{marine.colours}
\title{Marine colours}
\description{
  Create a perceptually contigous palette of R \code{"#RRGGBB"} colours,
  using hues typically associated with natural waters.
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
    vector, in which case the chroma is changed linearly from start to
    finish of the resulting palette.
  }
  \item{luminance}{
    Specifies the luminance (how bright should the colours be) of the
    colours constituting the palette. Typically, a two-element vector
    of real numbers between 0 and 1 to indicate linear change along the
    palette, but can also be a fixed number. See the \code{gamma}
    argument if you want to make the change non-linear.
  }
  \item{alpha}{
    Specifies the transparency of the colours of the palette. As above,
    can be a fixed number or a two-element vector in the range
    \eqn{[0,1]}. Typically, fully opaque (\code{alpha=1}) colours are
    used.
  }
  \item{gamma}{
    Provides the power coefficient for the luminance growth formula.
    May be useful when it is needed to sacrifice the perceptual
    linearity of the palette to provide more constrast for darker or
    brighter areas of the plot. The exact expression used is
    \code{luminance ^ gamma}. Typically, linear growth (\code{gamma=1})
    is preferred.
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
}
\references{
  Insired by \href{https://matplotlib.org/cmocean/#haline}{cmocean
  palette called \dQuote{haline}}, but using R's implementation of polar
  CIE-LUV colour space instead of CAM02-UCS.
}

\seealso{
 \code{\link{hcl}} for the colour space used,
 \href{http://www.mrao.cam.ac.uk/~dag/CUBEHELIX/}{CUBEHELIX} for a
 similar technique using BT.601 luminance coefficients and RGB colour
 space.
}

\examples{
  image(volcano, col = marine.colours(256))
}
\keyword{color}
