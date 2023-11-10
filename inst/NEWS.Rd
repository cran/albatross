\name{NEWS}
\title{News for package \pkg{albatross}}
\section{Changes in version 0.3-7 (2023-11-11)}{
  \subsection{Bug fixes}{\itemize{
    \item The previous version of the package was accidentally broken on
    versions of \R where \code{apply()} did not yet have a
    \code{simplify} argument. More strict release testing procedures
    will be applied in the future.
  }}
  \subsection{New features}{\itemize{
    \item New argument \code{like} to \code{reorder.feemparafac} and
    \code{rescale.feemparafac} makes it possible to order and scale the
    components to match another model. The package now re-exports the
    \code{rescale} generic from the \CRANpkg{multiway} package.
  }}
}
\section{Changes in version 0.3-6 (2023-07-17)}{
  \subsection{Potentially breaking changes}{\itemize{
    \item Fixed an error in \code{feemscatter(method='whittaker')} where
    the \code{nonneg} parameter could be ignored (turning the
    non-negativity protection off) if the wavelength grid was not sorted.
    This mistake has been present since version 0.3-0, where Whittaker
    smoothing had been introduced.

    \item Inverted the luminance of the \code{diverging.colours} palette,
    setting it to brighter for lower-magnitude values. The signs of the
    residuals are much easier to understand now.
  }}
  \subsection{Bug fixes}{\itemize{
    \item The \code{feemcube.feemflame} method is now actually exported.

    \item The \code{plot()} method for \code{feemflame} objects now
    correctly forwards its arguments to \pkg{lattice} functions for
    \code{type = 'both'}.
  }}
  \subsection{New features}{\itemize{
    \item New function \code{absindex} computes absorption coefficients,
    absorbance values, spectral slopes and their ratios.
  }}
}
\section{Changes in version 0.3-5 (2022-10-06)}{\itemize{
  \item When importing files with a non-default encoding,
  \code{feem.character} now avoids leaking a connection per file. (They
  were cleaned up anyway by the garbage collector, but with a warning.)

  \item The custom progress bar in \code{feemflame} is now re-drawn
  correctly on R GUIs.

  \item On \R \ifelse{html}{\out{&ge;}}{\eqn{\ge}{>=}} 4.2, equations in
  HTML help are typeset from LaTeX source, not manually translated HTML.

  \item \code{t()} on \code{feem} objects now results in an unclassed
  matrix, not a broken \code{feem} object.

  \item \code{feemscatter} now accepts scattering bands with asymmetric
  widths. On some fluorometers, this is needed in order to remove
  unwanted signal in anti-Stokes zone.

  \item New palette, \code{diverging.colours}, is used to plot residual
  objects.
}}
\section{Changes in version 0.3-4 (2022-02-25)}{\itemize{
  \item The \code{const} argument of \code{feemparafac} (and other
  functions that call it) now defaults to non-negativity. The previous
  default of an unconstrained model can be achieved by setting the
  argument to \code{NULL}.

  \item Fixed the bugs in \code{feemjackknife} and \code{feemsplithalf}
  which prevented them from working with 1-component models.

  \item Fixed the bugs preventing \code{feem} and \code{feemcube}
  objects from being properly indexed by their \code{dimnames}.

  \item The \code{data(feems)} dataset has been replaced with a
  synthetic one, with perfectly matching absorbance spectra.

  \item New function \code{feemflame} to model fluorescence and
  scattering signal at the same time.

  \item New function \code{feemindex} to calculate fluorescence indices
  and named peak values used in dissolved organic matter analysis.

  \item When the samples and the absorbance information don't match,
  \code{feemife} now provides an explanation of what went wrong.

  \item When an error happens during the processing of a \code{feem}
  inside a collection (a \code{feemcube} or a list of \code{feem}s) or
  during \code{feemlist} file import, the offending sample (or file) is
  now reported, and the error can be trapped without losing the
  traceback.

  \item The \code{feemcube} function can now return the original cube
  used to fit \code{feemparafac}, \code{feemjackknife},
  \code{feemsplithalf} objects.

  \item \code{feem} now supports text files from Edinburgh Instruments
  F900 software.

  \item \code{feemlist} can now take an importer function as an argument.

  \item \code{feemsplithalf} can limit the number of split-combine
  comparisons (the \code{splits} argument) or take a fixed list of
  halves to fit and compare (the \code{fixed} argument).

  \item HTML version of the manual now passes Tidy validation without
  warnings.
}}
\section{Changes in version 0.3-2 (2021-07-14)}{\itemize{
  \item Fixed a bug in numerical differentiation code for
  \code{feemscatter(method = 'whittaker')} that resulted in \eqn{d}-th
  order derivatives being \eqn{2^d} times smaller than they should have
  been. Existing optimal \code{lambda} values will have to be updated by
  dividing them by \eqn{2^{2d}}{2^(2d)}.

  \item When interpolating using \dQuote{loess} or \dQuote{kriging}
  methods, predicted negative values are replaced with zeroes. Previously,
  this happened for the whole FEEM, but only with \dQuote{loess}.

  \item \code{feemcube} now can also extract the original cube from
  \code{feemparafac} objects.
}}
\section{Changes in version 0.3-1 (2021-03-29)}{\itemize{
  \item Simplified the \code{feemlist(..., simplify.names = TRUE)} logic,
  hopefully preventing the test failure on some macOS machines.

  \item Made the \code{feem(..., 'table', fileEncoding = '...')} logic
  more robust, especially on older versions of \R where it could still
  fail when \code{fileEncoding} could not be represented in the locale
  encoding.

  \item \code{example(write.openfluor)} no longer leaves files it creates
  in the per-session temporary directory.
}}
\section{Changes in version 0.3-0 (2021-03-24)}{\itemize{
  \item Added a new \sQuote{whittaker} interpolation method for
  \code{feemscatter} and \code{feemgrid} based on Whittaker smoothing
  with configurable difference order(s) and optional iterative penalty
  against resulting negative values. This brings in a dependency on the
  recommended \pkg{Matrix} package.

  \item The \sQuote{pchip} interpolation method now offers a choice
  whether to interpolate the EEM by row, by column, or to do both and
  average the results (thanks \pkg{staRdom} for the trick). It also
  doesn't crash anymore when trying to 1-D interpolate less than 3
  defined points, falling back to linear interpolation instead.

  \item Added \code{feem} methods for strings and connections, making
  it possible to read EEMs from files in simple table formats, including
  those produced by Horiba and Panorama software.

  \item Added the \code{feemlist} generic to construct lists of
  \code{feem} objects from \pkg{eemR} or \pkg{EEM} objects or vectors
  of file and directory paths.

  \item Added the \code{groups} argument to \code{feemsplithalf},
  making it possible to stratify the resulting halves according to a
  given factor value.

  \item Added the \code{feemgrid} generic to interpolate FEEM(s) on a
  given wavelength grid.

  \item Added \code{coef} methods for \code{feemparafac},
  \code{feemsplithalf}, \code{feemjackknife} returning the coefficients
  and summary values (e.g. TCC for split-half) in \dQuote{long}
  \link{data.frame} format.

  \item Added function \code{write.openfluor} to export a fitted
  \code{feemparafac} object for the OpenFluor service.

  \item Equations in the HTML help pages are now typeset using HTML
  features instead of plain ASCII.

  \item Dependency on \code{R >= 3.3} is now acknowledged. The
  package doesn't actually pass \code{R CMD check} on prior versions
  of \R. Patches to lower the required \R version are welcome if they
  don't result in extensive rewrites.

  \item Other documentation improvements.
}}
\section{Changes in version 0.2-0 (2020-11-30)}{\itemize{
  \item Fixed a serious bug in \code{feemsplithalf} that would result
  in partially overlapping splits being compared. It has been present
  since v0.1-1. I apologise to any users of \pkg{albatross} who found
  out they have been bitten by this. Steps have been taken to make sure
  that other bugs like this don't still lurk in the code and to prevent
  similar bugs from occurring in the future.

  \item Potentially breaking change: \code{feemife} now accepts either
  named absorbance / optical path length data exactly matching unique
  names of samples, or a collection of exactly same length with either
  names missing.  Named collections of the same length as number of
  samples which don't contain all unique sample names are going to be
  rejected from now on.

  \item Saved \code{feemsplithalf} and \code{feemjackknife} results now
  should consume much less space and much less memory when loaded back.
  This is achieved by saving the FEEM cube object in an environment and
  only referencing it using vectors of indices instead of keeping copies
  of its subsets.

  \item New interpolation method: \dQuote{kriging}. Very
  resource-intensive, but results can be nicer-looking than either
  \code{loess} or \code{pchip}.

  \item More consistent progress bars for operations on collections of
  \code{feem}s. Scatter correction can now be run in parallel.

  \item Minor documentation improvements.
}}
\section{Changes in version 0.1-2 (2020-10-02)}{\itemize{
  \item Fixed a bug in \code{feemjackknife} and \code{feemsplithalf},
  which could manifest when \code{rescale} argument is not \code{3} and
  loading scales are wildly different, resulting in poor approximation
  for some of the modes.

  \item Made the documentation for PARAFAC/split-half/jack-knife plots
  more prominent.

  \item Made tests and examples much faster.

  \item Made the \code{gamma} argument in \code{marine.colours} much
  more useful.

  \item Minor improvements in documentation.
}}
