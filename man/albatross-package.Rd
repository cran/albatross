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
  \code{\link[base]{data.frame}}s for collections of them. There's a
  slight memory benefit to using lists of \code{\link{feem}} objects,
  but the difference shouldn't be noticeable, so there's nothing to
  worry about if you started with a \code{\link{feemcube}}.

  In order to compute PARAFAC, you need to convert your data into a
  \code{\link{feemcube}}. Whether you perform jack-knifing, split-half
  analysis, or PARAFAC itself, a copy of the data cube is kept together
  with the results and can be extracted back using the
  \code{\link{feemcube}} function. The resulting objects support a
  \code{plot} method (described in the same help page) and can give you
  the data as a few-column \code{\link[base]{data.frame}} using the
  \code{coef} method.

  Once the analysis is finished, the PARAFAC model can be exported for
  the OpenFluor database (\code{\link{write.openfluor}}) or stored as an
  \R object using standard \R tools (\code{\link[base]{save}} or
  \code{\link[base]{saveRDS}}).

  \packageIndices{albatross}
}
\author{\packageAuthor{albatross}}
\references{
  \reference{Murphy2013}

  \reference{Pucher2019}

  \reference{MPFC1970}

  \reference{Krylov2020}
}
\keyword{package}
\seealso{
  \code{\link{feem}}, \code{\link{feemlist}}, \code{\link{feemindex}},
  \code{\link{feemife}}, \code{\link{feemscatter}},
  \code{\link{feemgrid}}, \code{\link{feemcube}},
  \code{\link{feemscale}}, \code{\link{feemsplithalf}},
  \code{\link{feemparafac}}, \code{\link{feemjackknife}},
  \code{\link{feemflame}}, \code{\link{absindex}}.
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
