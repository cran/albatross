\name{feem}
\alias{feem}
\alias{feem.matrix}
\alias{feem.data.frame}
\alias{feem.character}
\alias{feem.connection}
\newcommand{\reg}{\ifelse{latex}{\out{\textregistered}}{\ifelse{html}{\out{&reg;}}{(R)}}}
\newcommand{\regsp}{\reg\if{latex}{\out{\space}}}
\title{
  Create a fluorescence excitation-emission matrix object
}
\description{
  Functions to create fluorescence excitation-emission matrix objects from
  \R matrices coupled with excitation and emission wavelengths, three-column
  \code{data.frame}s containing \eqn3{
    (\lambda_\mathrm{em}, \lambda_\mathrm{ex}, I)
  }{(\lambda.em, \lambda.ex, intensity)}{
    (&lambda;<sub>em</sub>, &lambda;<sub>ex</sub>, I)
  } tuples or files.
}
\usage{
  feem(x, \dots)
  \method{feem}{matrix}(x, emission, excitation, scale = 1, \dots)
  \method{feem}{data.frame}(
    x, scale = 1, emission = 'emission',
    excitation = 'excitation', intensity = 'intensity', \dots
  )
  \method{feem}{character}(x, format, \dots)
  \method{feem}{connection}(x, format, \dots)
}
\arguments{
  \item{x}{
    The source of the information to create a FEEM object from: a matrix,
    a three-column \code{data.frame}, a file path as a single string,
    or a \link{connection}.

    If converting a matrix, its rows should correspond to different
    fluorescence emission wavelengths specified in the \code{emission}
    argument; conversely, its columns should correspond to excitation
    wavelengths specified in the \code{excitation} argument.

    If converting a \code{data.frame}, it should have exactly three
    columns containing emission wavelengths, excitation wavelength,
    and intensity values. The names of the columns are expected to be
    \dQuote{emission}, \dQuote{excitation}, and \dQuote{intensity},
    respectively, but can be overridden using namesake arguments.

    If reading a single file by file path or connection, the \code{format}
    argument must specify the kind of file to parse, see below.
  }
  \item{emission}{
    If converting a matrix, this should be a vector of emission wavelengths,
    each wavelength corresponding to a row of the matrix.

    If converting a \code{data.frame}, this optional argument specifies the
    name of the column containing the emission wavelengths.
  }
  \item{excitation}{
    If converting a matrix, this should be a vector of excitation wavelengths,
    each wavelength corresponding to a column of the matrix.

    If converting a \code{data.frame}, this optional argument specifies the
    name of the column containing the excitation wavelengths.
  }
  \item{intensity}{
    If converting a \code{data.frame}, this optional argument specifies the
    name of the column containing the fluorescence intensities.
  }
  \item{scale}{
    The scale value of a EEM is preserved through the analysis procedure
    to divide the resulting score values after running PARAFAC. If
    the EEM has been pre-multiplied prior to creating the FEEM object,
    you can set the multiplier here.
  }
  \item{format}{\describe{
    \item{table}{
      The FEEM is assumed to be stored as a plain text matrix, readable
      using \code{\link{read.table}}, with the first column and the
      first row containing wavelengths. For example, it is possible to
      import CSV files obtained from a HORIBA Aqualog\regsp fluorometer
      by using \code{feem(\var{file}, 'table', sep = ',')}.

      Rows are assumed to correspond to emission wavelengths, columns
      are assumed to correspond to excitation wavelengths; if that's
      not the case, set the \code{transpose} argument to \code{TRUE}.

      If there are unmeasured points in the spectrum (e.g. the anti-Stokes
      area) but they are hard to distinguish from measured values
      (e.g. stored as zeroes instead of \code{NA}), specify their
      values as the \code{na} argument (numeric vector). The function
      will check for triangles filled with these values (such that a
      threshold \eqn3{\Delta}{Delta}{&#x2206;} exists where for all \eqn3{
        \lambda_\mathrm{em} - \lambda_\mathrm{ex} > \Delta
      }{lambda.em - lambda.ex > Delta}{
        &lambda;<sub>em</sub> - &lambda;<sub>ex</sub> &gt; &#x2206;
      } or \eqn3{
        \lambda_\mathrm{em} - \lambda_\mathrm{ex} < \Delta
      }{lambda.em - lambda.ex < Delta}{
        &lambda;<sub>em</sub> - &lambda;<sub>ex</sub> &lt; &#x2206;},
      \eqn3{
        X(\lambda_\mathrm{em}, \lambda_\mathrm{ex}) \in \mathtt{na}
      }{X(lambda.em, lambda.ex) \\\%in\\\% na}{
        X(&lambda;<sub>em</sub>, &lambda;<sub>ex</sub>) &isin; <code>na</code>})
      and replace them with \code{NA}s. If the unmeasured values
      are not stored as numbers, use the \code{na.strings} argument of
      \code{\link{read.table}} to specify them.

      The \code{fileEncoding} argument is treated specially (but
      preserving the \code{\link{read.table}} semantics regarding
      numbers), preventing errors in case the session charset cannot
      represent all characters from the file: the contents are first
      converted to UTF-8, then forwarded to \code{\link{read.table}}.
      This only works when \code{file} is a file path, not a connection.
      (When passing a connection with an encoding attribute set, make
      sure that the file contents after decoding can be represented in
      the session charset.)

      All other arguments are passed to \code{\link{read.table}}, with
      \code{fill} defaulting to \code{TRUE} instead of \code{FALSE}.
    }
    \item{panorama}{
      Read a \file{.dat} file as created by \dQuote{Panorama} software
      that comes with FLUORAT\reg-02-PANORAMA fluorometer.
      Such files contain a header describing the wavelength range, e.g.:
      \preformatted{   240.0   650.0     1.0    ; Emission(columns)
   230.0   320.0     5.0    ; Excitation(rows)}
      The header is followed by the intensity data as matrix,
      whitespace-separated. Missing points (anti-Stokes area) stored as
      \eqn{0} or \eqn{100} and are automatically filtered out on import.
      No additional parameters are accepted.
    }
    \item{F900txt}{
      Read an emission map text file created by \dQuote{F900} software
      that comes with Edinburgh Instruments fluorescence spectrometers.
      Separated by empty lines, these files are assumed to contain:
      \enumerate{
        \item Sample name
        \item Metadata of each emission scan, with names in the first
        column, including the excitation wavelengths
        \item Fluorescence values with emission wavelengths in the first
        column
      }

      The \code{fileEncoding} argument specifies the encoding of the
      text file and has a safe default of \dQuote{latin1}. It should
      typically correspond to the ANSI code page on the computer running
      F900 that was used to create the file. See \code{\link{iconvlist}()}
      for a list of encodings understood by \R.

      The \code{sep} argument specifies the column separator used in the
      file. A semicolon \code{;} is typically used, though other options
      are possible.
    }
  }}
  \item{\dots}{
    When converting matrices and \code{\link{data.frame}}s, extra
    arguments besides those specified above are not allowed.

    When reading the FEEM from a file, additional arguments may be passed
    to \code{format}-specific reading functions, see above.
  }
}
\value{
  A FEEM object is a matrix with the following attributes added:
  \item{emission}{
    Fluorescence emission wavelengths corresponding to the rows of the
    matrix, nm.
  }
  \item{excitation}{
    Fluorescence excitation wavelengths corresponding to the columns of
    the matrix, nm.
  }
  \item{dimnames}{
    Dimension names, copies of information above. Used only for
    presentation purposes.
  }
  \item{scale}{
    Scale factor, preserved through the analysis, which may be used
    later to undo the scaling. Initially \eqn{1}.
  }
}

\seealso{
  FEEM methods: \code{\link{plot.feem}}, \code{\link{as.data.frame.feem}},
  \code{\link{[.feem}}, \code{\link{feemgrid}}, \code{\link{feemife}},
  \code{\link{feemscale}}, \code{\link{feemscatter}}.
}
\examples{
  feem(matrix(1:40, ncol = 8), 1:5, 1:8)
  feem(
    data.frame(x = 1:10, y = 21:30, z = 31:40),
    emission = 'x', excitation = 'y', intensity = 'z'
  )
  feem(
    system.file('extdata/ho_aq.csv', package = 'albatross'),
    'table', sep = ','
  )
  feem(
    system.file('extdata/F900.txt', package = 'albatross'), 'F900txt'
  )
}
\keyword{methods}
\keyword{utilities}
