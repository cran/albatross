\name{absindex}
\alias{absindex}
\title{Functions of absorbance data}
\description{
  Calculate absorption coefficients and/or absorbance data at given
  wavelengths, spectral slopes, and their ratios.
}
\usage{
  absindex(
    x, abs.path, unit = c("log10", "m^-1"), out.A = 254,
    out.a = c(350, 355, 374, 443),
    out.a.ratio = list(c(250, 365), c(465, 665)),
    out.slope = list(c(275, 295), c(350, 400)),
    out.slope.ratio = list(c(275, 295, 350, 400)),
    out.slope.nrmse = FALSE
  )
}
\arguments{
  \item{x}{
    Absorption data, either a \code{\link[base]{list}} of two-column
    matrix-like objects with wavelengths in first column and values in
    second column, or a multi-column matrix-like object with wavelengths
    in first column and values in all other columns. Can be named in
    order to match with the \code{abs.path} values, see
    \code{\link{feemife}} for details.
  }
  \item{abs.path}{
    A numeric vector of optical path lengths for every spectrum, in
    centimetres. Defaults to \eqn{1} for all samples. If specified,
    should be either named, with names matching \code{x}, or unnamed and
    containing exactly the same number of path lengths.
  }
  \item{unit}{
    Specifies whether \code{x} contains absorbance (\verb{log10};
    \eqn3{A = \log_{10}\frac{I_0}{I}}{A = log10(I0 / I)}{
      A = log<sub>10</sub>(<sup>I<sub>0</sub></sup>/<sub>I</sub>)};
    unitless) or absorption coefficients (\verb{m^-1};
    \eqn3{
      \alpha = \frac{\ln \frac{I_0}{I}}{l}
    }{\alpha = log(I0 / I) / l}{
      &alpha; = <sup>1</sup>/<sub>l</sub>
      ln(<sup>I<sub>0</sub></sup>/<sub>I</sub>)};
    \eqn3{[\mathrm{m}^{-1}]}{[m^-1]}{[m<sup>-1</sup>]}).
  }
  \item{out.A}{
    Return absorbance values at the wavelengths given as a numeric
    vector.
  }
  \item{out.a}{
    Return absorption coefficients at the wavelengths given as a numeric
    vector.
  }
  \item{out.a.ratio}{
    Return ratios of absorption coefficients at the wavelengths given as
    a list of two-element numeric vectors. For every pair of
    wavelengths, \eqn3{\frac{\alpha(\lambda_1)}{\alpha(\lambda_2)}
    }{\alpha(\lambda[1]) / \alpha(\lambda[2])}{
      <sup>&alpha;(&lambda;<sub>1</sub>)</sup>/<sub>&alpha;(&lambda;<sub>2</sub>)</sub>
    } is returned.
  }
  \item{out.slope}{
    Return spectral slopes at wavelength ranges given as a list of
    two-element numeric vectors. See the \code{slope.method} for the
    description of how spectral slopes are computed.
  }
  \item{out.slope.ratio}{
    Return ratios of spectral slopes for pairs of wavelength ranges
    given as a list of four-element numeric vectors. For every list
    element, the value returned is \eqn3{
      \frac{S(\lambda_1, \lambda_2)}{S(\lambda_3, \lambda_4)}
    }{S(\lambda[1], \lambda[2]) / S(\lambda[3], \lambda[4])}{
      <sup>S(&lambda;<sub>1</sub>,
      &lambda;<sub>2</sub>)</sup>/<sub>S(&lambda;<sub>3</sub>,
      S(&lambda;<sub>4</sub>)</sub>}.
  }
  \item{out.slope.nrmse}{
    When computing slopes, also return the root-mean-square error for
    the models providing them, divided by the range of the response:
    \eqn3{
      \frac{1}{y_\mathrm{max} - y_\mathrm{min}}
      \sqrt{\sum_{i=1}^n\frac{(\hat{y}_i - y_i)^2}{n}}
    }{sqrt(mean((y.fitted - y.true)^2))/diff(range(y.true))}{
      <sup>1</sup>/<sub>y<sub>max</sub> - y<sub>min</sub></sub>
      &radic; &sum;
      <sup>(y&#x0302; - y)<sup>2</sup></sup>/<sub>n</sub>}.
  }
}
\details{
  Currently, the spectral slopes are calculated by fitting a linear
  model \eqn3{
      \ln \alpha = b_0 - b_1 \lambda
    }{log(\alpha) = b0 - b1 * \lambda}{
      ln &alpha; = b<sub>0</sub> - b<sub>1</sub> &lambda;
  } and returning \eqn3{b_1}{b1}{b<sub>1</sub>} as the slope. See
  \cite{Twardowski2004} for a discussion of the calculation methods for
  spectral slopes.

  Requested wavelengths missing from the original grid are interpolated
  using \code{\link[stats]{spline}}. \code{NA} values are returned
  outside the original wavelength range.
}
\value{
  A \code{\link[base]{data.frame}} with one row per sample, containing
  the following columns:

  \item{sample}{
    Names or numbers of the samples.
  }
  \item{A.<wavelength>}{
    Absorbance values, for every \v{wavelength} in \code{out.A}.
  }
  \item{a.<wavelength>}{
    Absorption coefficients, for every \v{wavelength} in \code{out.a}.
  }
  \item{aR.<wl[1]>.<wl[2]>}{
    Ratios of absorption coefficients, for every \v{wl} in
    \code{out.a.ratio}.
  }
  \item{S.<wl[1]>.<wl[2]>}{
    Spectral slopes, for every \v{wl} in \code{out.slope}.
  }
  \item{NRMSE.S.<wl[1]>.<wl[2]>}{
    If \code{out.slope.nrmse} is \code{TRUE}, root-mean-square errors
    normalised by the range of the absorption coefficients for the
    models providing spectral slopes.
  }
  \item{SR.<wl[1]>.<wl[2]>.<wl[3]>.<wl[4]>}{
    Ratios of spectral slopes, for every \v{wl} in
    \code{out.slope.ratio}.
  }
}
\references{\bibliography}

\seealso{\code{\link{feemife}}}
\examples{
  data(feems)
  absindex(absorp)
}
