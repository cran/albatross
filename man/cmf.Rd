\name{cmf}
\alias{cmf}
\alias{wcmls}
\title{Implementation notes for constrained matrix factorisation}
\description{
  \describe{
    \item{cmf}{
      Compute a low-rank matrix factorisation \eqn3{
        \min_{\mathbf A, \mathbf B}
        || (\mathbf X - \mathbf A \mathbf{B}^\top ) \circ \mathbf W ||_\mathrm F
      }{min over A,B ||(X - A \\\%*\\\% t(B)) * W||_F}{
        <span style="font-style:normal;">min</span><sub><b>A</b>,<b>B</b></sub>
        ||(<b>X</b> -
        <b>A</b><b>B</b><sup><span style="font-style:normal;">T</span></sup>)
        &#x2218; <b>W</b>||<sub>F</sub>
      } subject to weights \eqn3{\mathbf W}{W}{<b>W</b>} (set to \eqn{0}
      where \eqn3{\mathbf X}{X}{<b>X</b>} is not defined) and constraints
      on rows of \eqn3{\mathbf{A}, \mathbf{B}}{A, B}{<b>A</b>, <b>B</b>}.
    }
    \item{wcmls}{
      Solve the weighted multivariate least squares problem \eqn3{
        \min_\mathbf{B}
        || (\mathbf X - \mathbf A \mathbf{B}^\top) \circ \mathbf W ||_\mathrm F
      }{min over B ||(X - A \\\%*\\\% t(B)) * W||_F}{
        <span style="font-style:normal;">min</span><sub><b>B</b></sub>
        ||(<b>X</b> -
        <b>A</b><b>B</b><sup><span style="font-style:normal;">T</span></sup>)
        &#x2218; <b>W</b>||<sub>F</sub>
      } subject to constraints on rows of \eqn3{\mathbf B}{B}{<b>B</b>}.
    }
  }

  \emph{This is not a public interface.} Subject to change without
  further notice. Please do not call from outside \CRANpkg{albatross}.
}
\usage{
  cmf(
    X, nfac = 1,
    const = list(list(const = "nonneg"), list(const = "nonneg")),
    start = c("svd", "random"), ctol = 1e-04, maxiter = 10
  )
  \method{fitted}{cmf}(object, \dots)
  wcmls(X, A, W, \dots, struc = NULL)
}
\arguments{
  \item{X}{
    The matrix for a low-rank approximation.
  }
  \item{nfac}{
    The rank of the factorisation; the number of columns in matrices
    \eqn3{\mathbf A, \mathbf B}{A, B}{<b>A</b>, <b>B</b>}.
  }
  \item{const}{
    Constraints on the two matrices: a list of two lists of arguments to
    pass to \code{wcmls} when computing the corresponding matrix.
  }
  \item{start}{
    A \code{cmf} object to take the starting values from.
    Alternatively, a string: \describe{
      \item{svd}{
        Compute a truncated SVD \eqn3{
          \mathbf X =
          \mathbf U \,
          \mathrm{diag}(\sigma_1, \dots, \sigma_k)
          \, \mathbf{V}^\top
        }{X = U \\\%*\\\% diag(s) \\\%*\\\% t(V)}{
          <b>X</b> = <b>U</b>
          <span style="font-style:normal;">diag</span>(<b>&sigma;</b>)
          <b>V</b><sup><span style="font-style:normal;">T</span></sup>}.
        Use
        \eqn3{
          \mathbf A = \mathbf U \,
          \mathrm{diag}(\sqrt{\sigma_1}, \dots, \sqrt{\sigma_k})
        }{A = U \\\%*\\\% diag(sqrt(s))}{
          <b>A</b> = <b>U</b>
          <span style="font-style:normal;">diag</span>(&radic;<b>&sigma;</b>)},
        \eqn3{
          \mathbf B = \mathbf V \,
          \mathrm{diag}(\sqrt{\sigma_1}, \dots, \sqrt{\sigma_k})
        }{B = V \\\%*\\\% diag(sqrt(s))}{
          <b>B</b> = <b>V</b>
          <span style="font-style:normal;">diag</span>(&radic;<b>&sigma;</b>)
        } as the starting values.
      }
      \item{random}{
        Use uniformly distributed nonnegative starting values rescaled to
        be of comparable norms.
      }
    }
  }
  \item{ctol}{
    Given \eqn3{L = ||
      (\mathbf X - \mathbf A \mathbf{B}^\top ) \circ \mathbf W
    ||_\mathrm F}{L = ||(X - A \\\%*\\\% t(B)) * W||_F}{
      L = ||(<b>X</b> -
      <b>A</b><b>B</b><sup><span style="font-style:normal;">T</span></sup>)
      &#x2218; <b>W</b>||<sub>F</sub>},
    stop when
    \eqn3{
      \frac{|\Delta L|}{L} \le \mathtt{ctol}
    }{abs(diff(L))/L \le ctol}{
      <sup>|&Delta;L|</sup>/<sub>L</sub> &le; <code>ctol</code>}.
  }
  \item{maxiter}{
    Iteration number limit.
  }
  \item{object}{
    An object of class \code{cmf}.
  }
  \item{A}{
    The predictor matrix in the weighted multivariate least squares
    problem.
  }
  \item{W}{
    The weights matrix.
  }
  \item{\dots, struc}{\describe{
    \item{wcmls}{Passed to \code{\link{cmls}}.}
    \item{fitted.cmf}{Ignored.}
  }}
}
\details{
  The \CRANpkg{CMLS} package function \code{\link{cmls}} can solve
  constrained multivariate least squares problems of the form:

  \deqn3{
    \min_\mathbf{B} || \mathbf X - \mathbf A \mathbf B ||_\mathrm F
    = L(\mathbf X, \mathbf A, \mathbf B)
  }{min over B ||X - A \\\%*\\\% t(B)||_F = L(X, A, B)}{
    <span style="font-style:normal;">min</span><sub><b>B</b></sub>
    ||<b>X</b> - <b>A</b><b>B</b>||<sub>F</sub> =
    L(<b>X</b>, <b>A</b>, <b>B</b>)
  }

  We use it to solve a weighted problem. Let
  \eqn3{\mathbf X, \mathbf W}{X, W}{<b>X</b>, <b>W</b>} be
  \eqn3{(m \times n)}{(m * n)}{(m &times; n)} matrices,
  \eqn3{\mathbf A}{A}{<b>A</b>} be an
  \eqn3{(m \times k)}{(m * k)}{(m &times; k)} matrix,
  \eqn3{\mathbf B}{B}{<b>B</b>} be an
  \eqn3{(n \times k)}{(n * k)}{(n &times; k)}
  matrix,
  \eqn3{\mathbf{J}_{p,q}}{J(p,q)}{<b>J</b><sub>p,q</sub>} be a
  \eqn3{(p \times q)}{(p * q)}{(p &times; q)} matrix of ones:

  \deqn3{
    \min_\mathbf{B}
    || \mathbf W \circ (\mathbf X - \mathbf A \mathbf B^\top) ||_\mathrm F
    = \sum_{i,j} (
      w_{i,j} x_{i,j} - w_{i,j} \mathbf{a}_{i,\cdot} \mathbf{b}_{j,\cdot}^\top
    )^2
    = {}
  }{
    min over B ||W * (X - A \\\%*\\\% t(B))||_F =

  = sum over i,j (W[i,j] * X[i,j] - W[i,j] * sum(A[i,] * B[j,]))^2 =
  }{
    <span style="font-style:normal;">min</span><sub><b>B</b></sub>
    ||<b>W</b> &#x2218; (<b>X</b> -
      <b>AB</b><sup><span style="font-style:normal;">T</span></sup>)||<sub>F</sub>
    =
    &sum;<sub>i,j</sub> (w<sub>i,j</sub> x<sub>i,j</sub> -
    w<sub>i,j</sub> <b>a</b><sub>i,&middot;</sub>
    <b>b</b><sub>j,&middot;</sub><sup><span style="font-style:normal;">T</span></sup>)<sup>2</sup>
    =
  }
  \deqn3{
    {} = \sum_j ||
      \mathbf{w}_{\cdot,j} \circ \mathbf{x}_{\cdot,j} - (
        (\mathbf{w}_{\cdot,j} \mathbf{J}_{1,k}) \circ \mathbf A
      ) \mathbf{b}_{j,\cdot}^\top
    ||_\mathrm F
    = \sum_j L(
      \mathbf{w}_{\cdot,j} \circ \mathbf{x}_{\cdot,j},
      (\mathbf{w}_{\cdot,j} \mathbf{J}_{1,k}) \circ \mathbf A,
      \mathbf{b}_{j,\cdot}^\top
    )
  }{
  = sum over j ||W[,j] * X[,j] - (W[,j] * A) \\\%*\\\% B[j,]||_F =

  = sum over j L(W[,j] * X[,j], W[,j] * A, B[j,])}{
    = &sum;<sub>j</sub> ||
    <b>w</b><sub>&middot;,j</sub> &#x2218; <b>x</b><sub>&middot;,j</sub> -
    ((<b>w</b><sub>&middot;,j</sub> <b>J</b><sub>1,k</sub>) &#x2218; <b>A</b>)
    <b>b</b><sub>j,&middot;</sub><sup><span style="font-style:normal;">T</span></sup>
    ||<sub>F</sub>
    = &sum;<sub>j</sub>
    L(<b>w</b><sub>&middot;,j</sub> &#x2218; <b>x</b><sub>&middot;,j</sub>,
      (<b>w</b><sub>&middot;,j</sub> <b>J</b><sub>1,k</sub>) &#x2218; <b>A</b>,
    <b>b</b><sub>j,&middot;</sub><sup><span style="font-style:normal;">T</span></sup>)
  }

  Here, \eqn3{\mathbf{w}_{\cdot,j}}{W[,j]}{<b>w</b><sub>&middot;,j</sub>}
  and \eqn3{\mathbf{x}_{\cdot,j}}{X[,j]}{<b>x</b><sub>&middot;,j</sub>}
  are columns of \eqn3{\mathbf W}{W}{<b>W</b>} and
  \eqn3{\mathbf X}{X}{<b>X</b>}, while
  \eqn3{\mathbf{a}_{i,\cdot}}{A[i,]}{<b>a</b><sub>i,&middot;</sub>} and
  \eqn3{\mathbf{b}_{j,\cdot}}{B[j,]}{<b>b</b><sub>j,&middot;</sub>} are
  rows of \eqn3{\mathbf A}{A}{<b>A</b>} and \eqn3{\mathbf B}{B}{<b>B</b>},
  respectively. Thus, in the weighted case, the
  \eqn3{\mathbf B}{B}{<b>B</b>} matrix is determined row by row by
  calling the \code{\link{cmls}} function for pre-processed
  \eqn3{\mathbf A}{A}{<b>A</b>} matrix and columns of
  \eqn3{\mathbf X}{X}{<b>X</b>}.

  The problem we're actually interested in is a low-rank approximation
  of \eqn3{\mathbf X}{X}{<b>X</b>}. It doesn't have a unique solution,
  especially if the rank is more than \eqn{1}, unless we apply
  constraints and some luck. We solve it by starting with (typically)
  SVD and refining the solution with alternating least squares until it
  satisfies the constraints: \eqn3{
    \min_\mathbf{B}
    || (\mathbf X - \mathbf A \mathbf{B}^\top) \circ \mathbf W ||_\mathrm F
  }{min over B ||(X - A \\\%*\\\% t(B)) * W||_F}{
    <span style="font-style:normal;">min</span><sub><b>B</b></sub>
    ||(<b>X</b> -
    <b>A</b><b>B</b><sup><span style="font-style:normal;">T</span></sup>)
    &#x2218; <b>W</b>||<sub>F</sub>
  } and \eqn3{
    \min_\mathbf{A} ||
    (\mathbf{X}^\top - \mathbf B \mathbf{A}^\top) \circ \mathbf{W}^\top
    ||_\mathrm F
  }{min over A ||(t(X) - B \\\%*\\\% t(A)) * t(W)||_F}{
    <span style="font-style:normal;">min</span><sub><b>A</b></sub>
    ||(<b>X</b><sup><span style="font-style:normal;">T</span></sup> -
    <b>B</b><b>A</b><sup><span style="font-style:normal;">T</span></sup>)
    &#x2218;
    <b>W</b><sup><span style="font-style:normal;">T</span></sup>||<sub>F</sub>
  }.
}
\value{\describe{
  \item{cmf}{
    An list of class \code{cmf} containing the
    \eqn3{\mathbf A, \mathbf B}{A, B}{<b>A</b>, <b>B</b>} matrices.
  }
  \item{wcmls}{
    The \eqn3{\mathbf B}{B}{<b>B</b>} matrix solving the constrained weighted
    multivariate least squares problem.
  }
  \item{fitted.cmf}{
    A matrix reconstructed from its \code{nfac}-rank decomposition.
  }
}}
\references{
  \bibentry{'Article',
    author = c(
      person('A.', 'de Juan'),
      person('J.', 'Jaumot'),
      person('R.', 'Tauler')
    ),
    doi = '10.1039/c4ay00571f',
    title = 'Multivariate Curve Resolution (MCR). Solving the mixture analysis problem',
    journal = 'Anal. Methods',
    year = 2014,
    volume = 6,
    issue = 14,
    pages = '4964-4976',
    publisher = 'The Royal Society of Chemistry'
  }
}

\seealso{
  \code{\link{cmls}}; the \CRANpkg{ALS} package.
}
\examples{
  data(feems)
  z <- feemscatter(feems$a, rep(25, 4), 'omit')
  str(zf <- albatross:::cmf(unclass(z)))
  str(albatross:::fitted.cmf(zf))
}
\keyword{internal}
