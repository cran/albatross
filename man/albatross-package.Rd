\encoding{UTF-8}
\name{albatross-package}
\alias{albatross-package}
\alias{albatross}
\docType{package}
\title{
  \packageTitle{albatross}
}
\description{
  \if{html}{\out{<span style="text-align: right; font-size: small;">}}
  \if{latex}{\out{
    \begin{flushright}\{\small
  }}
  Day after day, day after day, \cr
  We stuck, nor breath nor motion; \cr
  As idle as a painted ship \cr
  Upon a painted ocean.

  Water, water, every where, \cr
  And all the boards did shrink; \cr
  Water, water, every where, \cr
  Nor any drop to drink.

  -- Samuel Taylor Coleridge, \emph{The Rime of the Ancient Mariner}

  \if{latex}{\out{
    \}\end{flushright}
  }}
  \if{html}{\out{</span>}}
  \packageDescription{albatross}
}
\details{
  \packageIndices{albatross}
}
\author{
Timur Labutin

Maintainer: Ivan Krylov
}
\references{
  K.R. Murphy, C.A. Stedmon, D. Graeber, R. Bro, Fluorescence spectroscopy and
  multi-way techniques. PARAFAC, Analytical Methods. 5 (2013) 6557.
  \doi{10.1039/c3ay41160e}

  M. Pucher, U. \enc{Wünsch}{Wuensch}, G. Weigelhofer, K. Murphy, T. Hein,
  D. Graeber, staRdom: Versatile Software for Analyzing Spectroscopic Data of
  Dissolved Organic Matter in R, Water. 11 (2019) 2366.
  \doi{10.3390/w11112366}

  J. Cleese, T. Jones, Albatross: Flavours of different sea birds, Journal
  of Flying Circus. 1.13 (1970) 7:05-7:45.
}
\keyword{package}
\seealso{
  \code{\link{feem}}, \code{\link{feemcube}}, \code{\link{feemscatter}},
  \code{\link{feemife}}, \code{\link{feemscale}}, \code{\link{feemparafac}},
  \code{\link{feemsplithalf}}, \code{\link{feemjackknife}}.
}
\examples{
  plot(x <- feem(matrix(1:42, 7), 400:406, 350:355))

  data(feems)

  dataset <- feemcube(feems, FALSE)[1:30*6, 1:9*6,]
  dataset <- feemscatter(dataset, rep(24, 4), 'pchip')
  dataset <- feemife(dataset, absorp)
  plot(dataset <- feemscale(dataset, na.rm = TRUE))

  \donttest{
    # takes a long time
    (sh <- feemsplithalf(
      dataset, nfac = 2:5, const = rep('nonneg', 3), splits = 4)
    )
    plot(sh)
    jk <- feemjackknife(dataset, nfac = 3, const = rep('nonneg', 3))
    plot(jk)
  }

  pf <- feemparafac(dataset, nfac = 2, const = rep('nonneg', 3))
  plot(pf)
}
