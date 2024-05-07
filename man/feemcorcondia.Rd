\name{feemcorcondia}
\alias{feemcorcondia}
\alias{print.feemcorcondia}
\title{Core consistency diagnostic for PARAFAC models}
\description{
  Compute the core consistency diagnostic (\dQuote{CORCONDIA}) by
  fitting a \dQuote{Tucker3} core array to the existing PARAFAC
  loadings.
}
\usage{
  feemcorcondia(
    model, divisor = c("nfac", "core"),
    kind = c('pinv', 'iterative', 'vec'), \dots
  )
  \method{print}{feemcorcondia}(x, \dots)
}
\arguments{
  \item{model}{
    A PARAFAC model returned by \code{\link{feemparafac}}.
  }
  \item{divisor}{
    The divisor used in computation of the CORCONDIA value, see Details.
  }
  \item{kind}{
    A string choosing the method used to compute the least squares
    Tucker3 core. Defaults to \dQuote{pinv} for PARAFAC models without
    missing data and \dQuote{iterative} for models where missing data is
    present. See Details.
  }
  \item{x}{
    An object returned by \code{feemcorcondia}.
  }
  \item{\dots}{\describe{
    \item{\code{feemcorcondia}}{
      For \code{kind = 'iterative'}, forwarded to
      \code{\link[stats]{optim}} (see Details). Otherwise, not allowed.
    }
    \item{\code{print.feemcorcondia}}{Ignored.}
  }}
}
\details{
  The \dQuote{Tucker3} model uses three loading matrices and a small
  three-way \dQuote{core array} to describe a larger three-way array:

  \deqn3{
    X_{i,j,k} = \sum_r \sum_s \sum_t A_{i,r} B_{j,s} C_{k,t} G_{r,s,t}
  }{X[i,j,k] = sum[r] sum[s] sum[t] (A[i,r] B[j,s] C[k,t] G[r,s,t])}{
    X<sub>i,j,k</sub> = &#x2211;<sub>r</sub> &#x2211;<sub>s</sub>
      &#x2211;<sub>t</sub>
      A<sub>i,r</sub> B<sub>j,s</sub> C<sub>k,t</sub> G<sub>r,s,t</sub>
  }

  It's easy to show that constraining \eqn3{
    G_{r,s,t} = 1_{r = s = t}}{G[r,s,t] = r == s == t}{
    G<sub>r,s,t</sub> = 1<sub>r = s = t</sub>
  } makes the Tucker3 model equivalent to a PARAFAC model. The core
  consistency diagnostic works by constraining the loading matrices of a
  Tucker3 model to the existing loading matrices from a PARAFAC model
  and estimating the core array. The closer the resulting
  \eqn3{\mathbf{G}}{G}{<b>G</b>} tensor is to a diagonal one, the
  better.

  Given the least-squares estimated core tensor
  \eqn3{\mathbf{G}}{G}{<b>G</b>}, the ideal core tensor \eqn3{
    T_{r,s,t} = 1_{r = s = t}}{T[r,s,t] = r == s == t}{
    T<sub>r,s,t</sub> = 1<sub>r = s = t</sub>} and the denominator
  \eqn{D}, the CORCONDIA metric is defined as follows:

  \deqn3{
    \left(
      1 - \frac{
        \sum_r \sum_s \sum_t (G_{r,s,t} - T_{r,s,t})^2
      }{D}
    \right) \cdot 100\\\%
  }{(1 - sum((G - T)^2)/D) * 100\\\%}{
    (1 - <sup>(&#x2211;<sub>r</sub> &#x2211;<sub>s</sub>
      &#x2211;<sub>t</sub> (G<sub>r,s,t</sub> -
        T<sub>r,s,t</sub>)<sup>2</sup>)</sup>/<sub>D</sub>)
    &middot; 100\\\%
  }

  The denominator can be chosen to be either \eqn3{
    \sum_r \sum_s \sum_t T_{r,s,t}^2
  }{sum(T^2)}{
    &#x2211;<sub>r</sub> &#x2211;<sub>s</sub>
      &#x2211;<sub>t</sub> T<sub>r,s,t</sub><sup>2</sup>}, which is
  equal to the number of factors in the model (\code{divisor = 'nfac'}),
  or \eqn3{
    \sum_r \sum_s \sum_t G_{r,s,t}^2
  }{sum(G^2)}{
    &#x2211;<sub>r</sub> &#x2211;<sub>s</sub>
      &#x2211;<sub>t</sub> G<sub>r,s,t</sub><sup>2</sup>}, which will
  avoid resulting negative values (\code{divisor = 'core'}).

  There are multiple ways how the least squares Tucker3 core can be
  computed. When no data is missing, the matricised form of the model
  can be used to derive the expression (\code{kind = 'pinv'}, the
  default in such cases):

  \deqn3{
    \mathbf{X} =
    \mathbf{A} \mathbf{G} (\mathbf{C} \otimes \mathbf{B})^\top
    + \epsilon
  }{X = A \\\%*\\\% G \\\%*\\\% (C \\\%x\\\% B)' + epsilon}{
    <b>X</b> = <b>A</b> <b>G</b> (<b>C</b> &otimes;
    <b>B</b>)<sup><span style="font-style:normal;">T</span></sup>
    + &epsilon;
  }
  \deqn3{
    \hat{\mathbf{G}} =
    \mathbf{A}^{+} \mathbf{X} (
      (\mathbf{C}^\top)^{+} \otimes
      (\mathbf{B}^\top)^{+}
    )
  }{Ghat = A^+ \\\%*\\\% X \\\%*\\\% (C'^+ \\\%x\\\% B'^+)}{
    <b>G&#x0302;</b> = <b>A</b><sup>+</sup> <b>X</b>
    ((<b>C</b><sup><span
      style="font-style:normal;">T</span></sup>)<sup>+</sup>
      &otimes;
      (<b>B</b><sup><span
      style="font-style:normal;">T</span></sup>)<sup>+</sup>)
  }

  With missing data present, a binary matrix of weights
  \eqn3{\mathbf{W}}{W}{<b>W</b>} appears:

  \deqn3{
    \min_\mathbf{G} \left\| \mathbf W \circ (
      \mathbf{A} \mathbf{G} (\mathbf{C} \otimes \mathbf{B})^\top
      - \mathbf{X}
    ) \right\|^2
  }{min over G ||W * (A \\\%*\\\% G \\\%*\\\% (C \\\%x\\\% B)' - X)||^2}{
    <span style="font-style:normal;">min</span><sub><b>G</b></sub>
    ||<b>W</b> &#x2218; (<b>A</b> <b>G</b> (<b>C</b> &otimes;
      <b>B</b>)<sup><span style="font-style:normal;">T</span></sup>
      - <b>X</b>)||<sup>2</sup>
  }

  A gradient-based method can be used to solve this problem iteratively
  without allocating too much memory, but care must be taken to ensure
  convergence. For \code{kind = 'iterative'} (which is the default for
  models with missing data), \code{\link[stats]{optim}} is used with
  parameters \code{method = 'L-BFGS-B'},
  \code{control = list(maxit = 1000, factr = 1)}. Warnings will be
  produced if the method doesn't indicate successful convergence.

  The problem can also be solved exactly by unfolding the tensor into a
  vector and skipping the elements marked as missing:

  \deqn3{
    \min_\mathbf{G} \left\|
      \mbox{vec}(
        \mathbf{A} \mathbf{G} (\mathbf{C} \otimes \mathbf{B})^\top
      )_{\mbox{non-missing}}
      - \mbox{vec}(\mathbf{X})_{\mbox{non-missing}}
    \right\|^2
  }{min over G norm(
    vec(A \\\%*\\\% G \\\%*\\\% (C \\\%x\\\% B)')[non-missing]
    - vec(X)[non-missing])
  )^2}{
    <span style="font-style:normal;">min</span><sub><b>G</b></sub>
    ||<span style="font-style:normal;">vec</span>(<b>A</b>
      <b>G</b> (<b>C</b> &otimes;
      <b>B</b>)<sup><span
      style="font-style:normal;">T</span></sup>)<sub>non-missing</sub>
    - <span
      style="font-style:normal;">vec</span>(<b>X</b>)<sub>non-missing</sub>||<sup>2</sup>
  }
  \deqn3{
    \mbox{vec}(\mathbf{C} \otimes \mathbf{B} \otimes \mathbf{A})
      _{\mbox{non-missing}}
    \times \mbox{vec}\:\mathbf{G}
    = \mbox{vec}(\mathbf{X})_{\mbox{non-missing}}
  }{vec(C \\\%x\\\% B \\\%x\\\% A)[non-missing] \\\%*\\\% vec(G) = vec(X)[non-missing]}{
    <span style="font-style:normal;">vec</span>(<b>C</b>
      &otimes; <b>B</b> &otimes; <b>A</b>)<sub>non-missing</sub> &times
    <span style="font-style:normal;">vec</span>(<b>G</b>) =
    <span
      style="font-style:normal;">vec</span>(<b>X</b>)<sub>non-missing</sub>
  }

  Unfortunately, when this method is used (\code{kind = 'vec'}), the
  left-hand side of the equation has the size of
  \eqn3{\mathbf{X}}{X}{<b>X</b>} times the number of components cubed,
  which grows very quickly for models with large numbers of components.
}
\value{
  A numeric scalar of class \code{feemcorcondia} with the following
  attributes:
  \item{divisor}{
    The \code{divisor} argument, expanded to one of the valid options.
  }
  \item{core}{
    A three-way array containing the least-squares Tucker core for the
    given PARAFAC model.
  }
}
\references{\reference{Bro2003-CORCONDIA}}
\seealso{\code{multiway::\link[multiway]{corcondia}}}
\examples{
  data(feems)
  cube <- feemscale(feemscatter(cube, c(20, 14)), na.rm = TRUE)
  # kind = 'vec' is exact but may take a lot of memory
  feemcorcondia(feemparafac(cube, nfac = 3, ctol = 1e-4), kind = 'vec')
  # kind = 'iterative' used by default for models with missing data
  feemcorcondia(feemparafac(cube, nfac = 4, ctol = 1e-4))
}
\keyword{utilities}
\keyword{models}
