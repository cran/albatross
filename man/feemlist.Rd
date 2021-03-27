\name{feemlist}
\alias{feemlist}
\alias{feemlist.character}
\alias{feemlist.eemlist}
\alias{feemlist.EEM}
\title{Create lists of FEEM objects}
\description{
  Convert vectors of file names or objects from other packages (such as
  \pkg{eemR} or \pkg{EEM}) into flat named lists of \code{\link{feem}}
  objects.
}
\usage{
  feemlist(x, \dots)
  \method{feemlist}{character}(
    x, format, pattern = NULL, recursive = TRUE, ignore.case = FALSE,
    simplify.names = TRUE, \dots
  )
  \method{feemlist}{eemlist}(x, \dots)
  \method{feemlist}{EEM}(x, \dots)
}
\arguments{
  \item{x}{
    A character vector containing names of files and directories to
    import using \code{\link{feem}}.

    Alternatively, an \code{eemlist} object from \pkg{eemR} package or
    an \code{EEM} object from \pkg{EEM} package.
  }
  \item{format}{
    Corresponds to the \code{format} argument of \code{\link{feem}}.
    Currently, one format is assumed for all files to be imported.
  }
  \item{pattern, recursive, ignore.case}{
    These options are passed to \code{\link{list.files}} for directories
    encountered in \code{x} and can be used to e.g. only choose files
    with a given suffix in the name. Note the non-default value for the
    \code{recursive} option.
  }
  \item{simplify.names}{
    If \code{TRUE} (default), split resulting names by the path separator
    (\code{/}, also \code{\\} on Windows) and remove leading components
    that have the same value for all samples, but leave at least one
    component. See Details on how this is related to name generation.
  }
  \item{\dots}{
    When importing files, remaining options are passed to
    \code{\link{feem}}. Otherwise, no options are allowed.
  }
}
\details{
  Names of \code{x} are preserved; if \code{x} is not named, names are
  assigned from the values of \code{x} itself, and so are empty names
  in partially-named \code{x}. Every directory in \code{x} is replaced
  with its contents (as returned by \code{\link{list.files}}), their
  names obtained by concatenating the name of the directory element with
  their paths inside the directory (with \code{.Platform$file.sep} as a
  separator). For example, when importing \code{x = c('foo' = 'bar')}
  with directory \file{bar} containing \file{baz.txt}, resulting name
  would be \file{foo/baz.txt}.

  When importing many files from the same directory, the
  \code{simplify.names} option could is useful to avoid duplication in
  resulting names. For example,
  \code{feemlist('.', simplify.names = FALSE)} results in a list with
  all names starting with \code{./}, while
  \code{feemlist('foo/bar/baz', simplify.names = TRUE)} (default) would
  shave off all three common path components and the separators.

  Mixing files and directories in \code{x} will most likely not preserve
  the order of the elements.

  \emph{Note}: Please don't rely on this mechanism behaving exactly as
  specified as it may be changed in the future versions.
}
\value{
  A flat named list of \code{\link{feem}} objects.
}
\seealso{
  \code{\link{feem}}; the packages \pkg{eemR} and \pkg{EEM}.
}
\examples{
  feemlist(
    system.file('extdata/pano2.txt', package = 'albatross'),
    'table', transpose = TRUE, na = 0
  )
  if (requireNamespace('eemR')) feemlist(eemR::eem_read(
    system.file('extdata/ho_aq.csv', package = 'albatross'),
    import_function='aqualog'
  ))
  if (requireNamespace('EEM')) feemlist(EEM::readEEM(
    system.file('extdata/ho_aq.dat', package = 'albatross')
  ))
}
\keyword{methods}
\keyword{utilities}
\keyword{IO}
\keyword{file}
