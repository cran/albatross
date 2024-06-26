\encoding{UTF-8}
\name{feemflame}
\alias{feemflame}
\alias{fitted.feemflame}
\alias{residuals.feemflame}
\alias{coef.feemflame}
\alias{plot.feemflame}
\title{
  Fluorescence and scAttering Model Estimation
}
\description{
  Given a FEEM cube, model the fluorescence and the scattering signals
  at the same time as a sum of a PARAFAC model and a low-rank
  unfolded matrix factorisation.
}
\usage{
  feemflame(
    X, ffac, sfac, maxiter = 32, widths = rep(25, 4), Raman.shift = 3400,
    ctol = 1e-04, progress = TRUE, control.parafac, control.cmf
  )
  \method{fitted}{feemflame}(object, \dots)
  \method{residuals}{feemflame}(object, \dots)
  \method{coef}{feemflame}(
    object, type = c(
      "fluorescence",
      "scores", "loadings", "emission", "excitation", "samples",
      "scattering", "sc.scores", "sc.loadings"
    ), \dots
  )
  \method{plot}{feemflame}(
    x, type = c('both', 'fl.image', 'fl.lines'), \dots
  )
}
\arguments{
  \item{X}{A \code{\link{feemcube}} object.}
  \item{ffac}{
    The number of trilinear components used to model fluorescence,
    passed to \code{\link{feemparafac}}.
  }
  \item{sfac}{
    The number of bilinear (low-rank matrix factorisation) components
    used to model the scattering signal.
  }
  \item{maxiter}{
    Maximum number of alternating PARAFAC and constrained matrix
    factorisation iterations.
  }
  \item{widths}{
    Widths of the scattering regions, like in \code{\link{feemscatter}}:
    A numeric vector of length 4 containing the widths (in nm) of the
    scattering signal, in the following order:
    \enumerate{
      \item Rayleigh scattering
      \item Raman scattering
      \item Rayleigh scattering, \eqn3{2\lambda}{2 lambda}{2&lambda;}
      \item Raman scattering, \eqn3{2\lambda}{2 lambda}{2&lambda;}
    }
  }
  \item{Raman.shift}{
    Raman shift of the scattering signal, in
    \eqn3{\textrm{cm}^{-1}}{cm^-1}{cm<sup>-1</sup>}, like in
    \code{\link{feemscatter}}.
  }
  \item{ctol}{
    Given \eqn3{
      L = ||\mathbf X - \hat{\mathbf X}||^2
    }{L = ||X - X.hat||^2}{
      L = ||<b>X</b> - <b>X&#x0302;</b>||<sup>2</sup>},
    stop when \eqn3{
      \frac{|\Delta L|}{L} \le \mathtt{ctol}
    }{abs(diff(L))/L \le ctol}{
      <sup>|&Delta;L|</sup>/<sub>L</sub> &le; <code>ctol</code>}.
  }
  \item{progress}{
    Print progress information on the console, including the iteration
    number, relative sum of squared residuals, and relative change in
    sum of squared residuals.
  }
  \item{control.parafac, control.cmf}{
    Named lists of additional arguments to be passed to the underlying
    functions. Both default to \code{list(ctol = 1e-4, maxit = 10)},
    which makes both steps to run for 10 iterations or stop
    when the absolute change in \eqn3{R^2}{R^2}{R<sup>2</sup>} is less
    than \eqn3{10^{-4}}{1e-4}{10<sup>-4</sup>}, whichever happens
    sooner.
  }
  \item{object, x}{
    A \code{feemflame} object.
  }
  \item{type}{\describe{
    \item{coef}{
      Determines the type of coefficients to return: \describe{
        \item{fluorescence}{
          Equivalent to calling \code{\link{coef.feemparafac}} on the
          fluorescence model (default).
        }
        \item{scores, loadings, emission, excitation, samples}{
          Equivalent to calling \code{\link{coef.feemparafac}} on the
          fluorescence model and passing the respective \code{type}
          argument.
        }
        \item{sc.scores}{
          A \code{\link[base]{data.frame}} containing the following columns:
          \describe{
            \item{sample}{Sample numbers or names.}
            \item{value}{Scattering intensity value for a given factor.}
            \item{factor}{The number of the scattering component.}
          }
        }
        \item{sc.loadings}{
          A \code{\link[base]{data.frame}} containing the following columns:
          \describe{
            \item{emission, excitation}{
              The wavelengths corresponding to the value of the
              scattering profile.
            }
            \item{value}{Scattering intensity value for a given factor.}
            \item{factor}{The number of the scattering component.}
          }
        }
        \item{scattering}{
          A list with names \dQuote{scores} and \dQuote{loadings}
          containing results of \code{coef(object, 'sc.scores')} and
          \code{coef(object, 'sc.loadings')}, respectively.
        }
      }
    }
    \item{plot}{
      Describes the kind of plot to produce: \describe{
        \item{both}{
          Plot the loadings of the fluorescence and scattering models as
          false colour images.
        }
        \item{fl.image, fl.lines}{
          Equivalent to calling \code{\link{plot.feemparafac}} on the
          fluorescence model with the argument \dQuote{image} or
          \dQuote{lines}, respectively.
        }
      }
    }
  }}
  \item{\dots}{
    No other parameters are allowed.
  }
}
\details{
  FLAME models the input data as a sum of fluorescence signal (PARAFAC
  model) and scattering signal (low rank model):

  \deqn3{
    X_k(\lambda^\mathrm{em}_i, \lambda^\mathrm{ex}_j) =
    \underbrace{\sum_p A_{i,p} B_{j,p} C_{k,p}}_{\mbox{fluorescence}}
    + \underbrace{\sum_q S_{i,j,q} D_{k,q}}_{\mbox{scattering}}
  }{X[i,j,k] = sum(A[i,] * B[j,] * C[k,]) + sum(S[i,j,] * D[k,])}{
    X<sub>k</sub>(&lambda;<sup>em</sup><sub>i</sub>,
    &lambda;<sup>ex</sup><sub>j</sub>) =
    &sum;<sub>p</sub> A<sub>i,p</sub> B<sub>j,p</sub> C<sub>k,p</sub>
    + &sum;<sub>q</sub> S<sub>i,j,q</sub> D<sub>k,q</sub>
  }

  The function alternates between fitting the PARAFAC model on the
  dataset with scattering signal subtracted and fitting the low-rank
  model on the dataset with fluorescence signal subtracted. The PARAFAC
  model is fitted using the \code{\link{feemparafac}} function.  The
  low-rank model is fitted by means of unfolding the wavelength
  dimensions into one, resulting in a matrix, followed by the same
  alternating least squares procedure as done in multivariate curve
  resolution. Both models are constrained to result in non-negative
  factors.

  The low-rank model is additionally constrained to zero outside the
  scattering region. The scattering region is defined the same way as in
  \code{\link{feemscatter}}, using the \code{widths} and the
  \code{Raman.shift} arguments.

  Initial PARAFAC model is fitted with the scattering region set to
  missing. The low-rank model is initialised with truncated singular
  value decomposition forced to be non-negative.
}
\value{
  \item{feemflame}{
    An object of class \code{feemflame}, which is a list containing the
    following components: \describe{
      \item{fl}{
        A \code{\link{feemparafac}} object containing the fluorescence
        part of the model.
      }
      \item{sc}{
        An object of internal class \code{cmf}. Please don't rely on its
        structure.
      }
    }
  }
  \item{fitted.feemflame}{
    A \code{\link{feemcube}} object containing the part of \code{X}
    fitted by the model.
  }
  \item{residuals.feemparafac}{
    A \code{\link{feemcube}} object equal to \eqn3{
      \mathbf{X} - \hat{\mathbf{X}}
    }{X - fitted(X)}{<b>X</b> - <b>X&#x0302;</b>}, with an extra class
    \code{feem.resid} set. Objects of this class are plotted with a
    different default palette, see \code{\link{plot.feem.resid}}.
  }
  \item{coef.feemflame}{
    See the description of the \code{type} argument.
  }
}
% NB: old versions of R don't work well with stage=install in the \doi
% definition in newer R
\RdOpts{stage=build}
\references{
  \reference{Tauler1998}

  \reference{Krylov2021}
  \url{https://files.libs.chem.msu.ru/~ivan/SSC17/P13.pdf}.
}
\note{
  The structure of the \code{feemflame} object, the initialisation, and
  the constraints may be subject to change in a future version.
}
\seealso{\code{\link{feemparafac}}, \code{\link{feemcube}}}
\examples{\donttest{
  data(feems)
  cube <- feemscale(cube)
  factors <- feemflame(cube, ffac = 3, sfac = 1)
  str(coef(factors))
  str(coef(factors, 'scattering'))
  plot(factors)
}}
\keyword{method}
\keyword{models}
\keyword{multivariate}
