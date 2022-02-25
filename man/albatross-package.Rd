\encoding{UTF-8}
\name{albatross-package}
\alias{albatross-package}
\alias{albatross}
\docType{package}
\title{
  \packageTitle{albatross}
}
\description{
  \if{html}{\out{<div style="text-align: right; font-size: small;">}}
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
  \if{html}{\out{</div>}}
  \packageDescription{albatross}
}
\details{
  In order to work with your data, create \code{\link{feem}} and/or
  \code{\link{feemcube}} objects from files or matrix or array objects.
  Use \code{\link{feemlist}} to import files in bulk. If your files
  aren't in one of the formats supported by \code{\link{feem}} but you
  can read their contents by other means, you can supply an importer
  function to \code{\link{feemlist}}; it should take a file name and
  return the corresponding \code{\link{feem}} object.

  Operations that can be performed on the objects include plotting
  (\code{\link{plot.feem}}), calculation of fluorescence indices
  (\code{\link{feemindex}}), inner-filter effect correction
  (\code{\link{feemife}}), handling of scattering signal
  (\code{\link{feemscatter}}), changing the wavelength grid of the data
  by means of interpolation (\code{\link{feemgrid}}), and scaling
  (\code{\link{feemscale}}). Scaling may be automatically undone after
  performing the PARAFAC decomposition so that the resulting scores
  would correspond to the data as it was before the scaling.

  All processing functions can take individual \code{\link{feem}}
  objects, lists of them, or \code{\link{feemcube}} objects and return
  values of the appropriate kind. For example, \code{\link{feemscatter}}
  always returns an object of the same class but with the scattering
  signal handled, while \code{\link{feemindex}} returns named numeric
  vectors for individual \code{\link{feem}}s but
  \code{\link{data.frame}}s for collections of them. There's a slight
  memory benefit to using lists of \code{\link{feem}} objects, but the
  difference shouldn't be noticeable, so there's nothing to worry about
  if you started with a \code{\link{feemcube}}.

  In order to compute PARAFAC, you need to convert your data into a
  \code{\link{feemcube}}. Whether you perform jack-knifing, split-half
  analysis, or PARAFAC itself, a copy of the data cube is kept together
  with the results and can be extracted back using the
  \code{\link{feemcube}} function. The result objects support a
  \code{\link{plot}} method, or can give you the data as a few-column
  \code{\link{data.frame}} using the \code{\link{coef}} method.

  Once the analysis is finished, the PARAFAC model can be exported for
  the OpenFluor database (\code{\link{write.openfluor}}) or stored as an
  R object using standard R tools (\code{\link{save}} or
  \code{\link{saveRDS}}).

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
  \code{\link{feem}}, \code{\link{feemlist}}, \code{\link{feemindex}},
  \code{\link{feemife}}, \code{\link{feemscatter}},
  \code{\link{feemgrid}}, \code{\link{feemcube}},
  \code{\link{feemscale}}, \code{\link{feemsplithalf}},
  \code{\link{feemparafac}}, \code{\link{feemjackknife}}.
}
\examples{
  data(feems)

  dataset <- feemcube(feems, FALSE)
  dataset <- feemscatter(dataset, rep(24, 4), 'pchip')
  dataset <- feemife(dataset, absorp)
  plot(dataset <- feemscale(dataset, na.rm = TRUE))

  \donttest{
    # takes a long time
    (sh <- feemsplithalf(cube, nfac = 2:5, splits = 4))
    plot(sh)
    jk <- feemjackknife(cube, nfac = 3)
    plot(jk)
  }

  pf <- feemparafac(cube, nfac = 3)
  plot(pf)
}
