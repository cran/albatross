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
  \bibentry{'Article',
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

  \bibentry{'Article',
    author = c(
      person('Matthias', 'Pucher'),
      # my kingdom for a backslash!
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

  \bibentry{'Article',
    author = c(
      person('John', 'Cleese'),
      person('Terry', 'Jones')
    ),
    title = 'Albatross: Flavours of different sea birds',
    journal = 'Journal of Flying Circus',
    year = 1970,
    volume = '1.13',
    pages = '7:05-7:45'
  }
}
\keyword{package}
\seealso{
  \code{\link{feem}}, \code{\link{feemlist}}, \code{\link{feemife}},
  \code{\link{feemscatter}}, \code{\link{feemgrid}},
  \code{\link{feemcube}}, \code{\link{feemscale}},
  \code{\link{feemsplithalf}}, \code{\link{feemparafac}},
  \code{\link{feemjackknife}}.
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
