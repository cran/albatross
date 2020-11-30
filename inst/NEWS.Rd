\name{NEWS}
\title{News for package \pkg{albatross}}
\section{Changes in version 0.2-0 (2020-11-30)}{\itemize{
  \item Fixed a serious bug in \code{feemsplithalf} that would result
  in partially overlapping splits being compared. It has been present
  since v0.1-1. I apologise to any users of \pkg{albatross} who found
  out they have been bitten by this. Steps have been taken to make sure
  that other bugs like this don't still lurk in the code and to prevent
  similar bugs from occuring in the future.

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

  \item New interpolation method: \code{kriging}. Very resource-intensive,
  but results can be nicer-looking than either \code{loess} or
  \code{pchip}.

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
