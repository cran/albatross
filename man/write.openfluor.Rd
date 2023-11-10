\name{write.openfluor}
\alias{write.openfluor}
\title{Export a PARAFAC model for the OpenFluor database}
\description{
  Prepares a fitted PARAFAC model for submission to OpenFluor - an online
  spectral database of fluorescence by environmental organic compounds.
}
\usage{
write.openfluor(
  model, filename, name = "?", creator = "?", doi = "?",
  reference = "?", unit = "?", toolbox =, date =, fluorometer = "?",
  constraints =, validation = "?", methods = "?", preprocess = "?",
  sources = "?", ecozones = "?", description = "",
  shift = FALSE, scale = TRUE
)
}
\arguments{
  \item{model}{A \code{\link{feemparafac}} object.}
  \item{filename}{
    Path to the text file to create from the \code{model} argument.
  }
  \item{name}{
    Short name of the model.
  }
  \item{creator}{
    Name of the creator of the model.
  }
  \item{doi}{
    Digital object identifier of the referenced source. Can also be
    \dQuote{ISBN:\ldots} for books.
  }
  \item{reference}{
    Full citation for the referenced source using the following style:
    \dQuote{Author AA, Author BB, Author CC, (year), \sQuote{Title},
    Journal Abbrev, Vol, pages}.
  }
  \item{unit}{
    Units the fluorescence was measured in. Typically, one of \dQuote{RU},
    \dQuote{QSE}, \dQuote{AU}.
  }
  \item{toolbox}{
    Defaults to
    \dQuote{albatross \v{version}, multiway \v{version}}.
  }
  \item{date}{
    Defaults to today, in \dQuote{\v{yyyy}-\v{mm}-\v{dd}} format.
  }
  \item{fluorometer}{
    The model of the instrument that produced the data.
  }
  \item{constraints}{
    Constraints applied to the PARAFAC model. Defaults to
    \code{model$const}, but please edit it to a more human-readable form.
  }
  \item{validation}{
    Validation method used for the PARAFAC model, examples include:
    \dQuote{Split-Half Analysis}, \dQuote{core-consistency}.
  }
  \item{methods}{
    The sequence of steps taken to handle the samples and to ensure
    proper fluorescence intensity measurement. Examples include: \itemize{
      \item Sampling: Filtration GF/F
      \item Sampling: Filtration \v{x} um
      \item Sampling: samples frozen
      \item Instrument spectral bias correction: Ex
      \item Instrument spectral bias correction: Em
      \item Instrument spectral bias correction: Ex & Em
      \item Inner filter effect correction: absorbance method
      \item Inner filter effect correction: dilution
      \item Inner filter effect correction: CDA
      \item Inner filter effect correction: \v{other (please describe)}
      \item Internal calibration: Raman Peak area
      \item Internal calibration: Raman Peak height
      \item Internal calibration: Blank Subtraction
      \item External calibration: Quinine Sulphate dilution series
      \item External calibration: STARNA reference standards
      \item External calibration: NIST reference standards
      \item External calibration: \v{other (please describe)}
    }
  }
  \item{preprocess}{
    PARAFAC-specific pre-processing steps applied to the dataset. Examples
    include (but are not limited to): \itemize{
      \item Outliers removed
      \item Scatter region excised (replaced with NaNs)
      \item Scatter region smoothed (replaced with interpolated values)
      \item Sample mode normalised to DOC concentration
      \item Sample mode normalised to unit variance
    }
  }
  \item{sources}{
    Should preferably include one or more of the following keywords:
    \itemize{
      \item river
      \item stream
      \item lake
      \item wetland
      \item reservoir
      \item estuary
      \item ocean - coastal and shelf seas
      \item ocean - surface off-shore
      \item ocean - deep off-shore
      \item freshwater
      \item seawater
      \item groundwater
      \item wastewater
      \item drinking water
      \item treated water
      \item recycled water
      \item ballast water
      \item sediment
      \item mudflat
      \item mangrove
      \item aquarium
      \item mesocosm
    }
  }
  \item{ecozones}{
    List all major or minor terrestrial, freshwater and
    marine ecozones and ecoregions that apply. The full set
    of possible options is too large to include here, but see
    \url{https://en.wikipedia.org/wiki/Lists_of_ecoregions} for a source
    of inspiration.
  }
  \item{description}{
    Brief description of the model and its source data in \eqn{\le 256}
    characters.
  }
  \item{shift, scale}{
    If \code{shift} is specified (default \code{FALSE}), the loadings
    are first shifted by subtracting \code{min(x)} to ensure that the
    minimal value is \eqn{0}.

    If \code{scale} is specified (default \code{TRUE}), the loadings
    are then rescaled by dividing by \code{max(x)} so that the maximal
    value is \eqn{1}.

    Note that OpenFluor clamps values outside the \eqn{[0, 1]} range
    and uses scale-invariant (but \emph{not} shift-invariant) Tucker's
    congruence coefficient to find matches.
  }
}
\details{
  Provided the \code{model} and the \code{filename} arguments, this
  function exports the loadings into a file that passes OpenFluor syntax
  check and is suitable for further editing. Alternatively, some or all
  of the fields may be specified programmatically.

  The fields \code{constraints}, \code{methods}, \code{preprocess},
  \code{sources}, \code{ecozones} can be specified as character vectors
  (to be comma-separated on output); others should be single strings.
}
\references{
  \reference{Murphy2014}

  \url{https://openfluor.lablicate.com/}
}

\seealso{\code{\link{feemparafac}}}
\examples{
  data(feems)
  cube <- feemscale(feemscatter(cube, rep(14, 4)), na.rm = TRUE)
  factors <- feemparafac(cube, nfac = 3)
  # all defaults
  write.openfluor(factors, f1 <- tempfile(fileext = '.txt'))
  if (interactive()) file.show(f1)
  unlink(f1)
  # all non-default arguments
  write.openfluor(
    factors, f2 <- tempfile(fileext = '.txt'), name = 'example',
    creator = 'J. Doe', doi = '10.1000/1', reference = paste(
      'Upper D, (1973),', "'The unsuccessful self-treatment of a case",
      "of \"writer's block\"',", 'J Appl Behav Anal, 7(3), 497'
    ), unit = 'AU', toolbox = 'all calculations done by hand',
    date = '2038-01-19', fluorometer = 'Acme Fluor-o-matic 9000',
    constraints = 'non-negative', validation = 'prior knowledge',
    methods = 'Instrument spectral bias correction: Ex & Em',
    preprocess = 'Scatter region excised (replaced with NaNs)',
    sources = 'freshwater', ecozones = 'Balkash',
    description = 'not a real model', shift = FALSE, scale = TRUE
  )
  if (interactive()) file.show(f2)
  unlink(f2)
}
\keyword{IO}
\keyword{file}
