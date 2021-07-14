\name{whittaker2}
\alias{whittaker2}
\alias{diffmat}
\alias{vandermonde}
\title{
  Implementation notes for Whittaker smoothing and interpolation of surfaces
}
\description{
  Smooth [1] or estimate the baseline [2] of surface measured on an
  arbitrary grid by minimising a sum of penalties. Combined difference
  orders and two different methods of preventing negative values in the
  output are supported.

  \emph{This is not a public interface.} Subject to change without
  further notice. Please do not call from outside \pkg{albatross}.
}
\usage{
  whittaker2(x, y, z, lambda, d, p, logscale, nonneg)
  diffmat(x, y, d)
  vandermonde(x0)
}
\arguments{
  \item{x}{Grid values along the rows of \code{z}.}
  \item{y}{Grid values along the columns of \code{z}.}
  \item{z}{
    Matrix containing the surface values to smooth or \code{NA}s to
    interpolate.
  }
  \item{lambda}{
    A vector of smoothness penalties, one for every difference order.
    Must be of the same length as \code{d}.
  }
  \item{d}{\describe{
    \item{\code{whittaker2}}{
      A vector of difference orders corresponding to elements of
      \code{lambda}.
    }
    \item{\code{diffmat}}{
      Difference order, an integer scalar.
    }
  }}
  \item{p}{
    If not missing, use the asymmetric penalty method [2] to estimate
    the baseline by penalising the differences with weight \eqn{p} if
    \eqn3{\hat{z} < z}{z.hat < z}{z&#x0302; &lt; z} and \eqn{1 - p}
    otherwise. Typically, values around
    \eqn3{10^{-3}}{1e-3}{10<sup>-3</sup>} are used.
  }
  \item{logscale}{
    If not \code{NA}, rescale \code{z} between \code{logscale} and
    \eqn{1} and take its logarithm before smoothing, then undo it. This
    is supposed to prevent negative numbers from appearing in the
    interpolation, but may damage the surface shape in the process.
  }
  \item{nonneg}{
    If not 0, for every resulting negative value in the interpolated
    surface, add a penalty of \code{nonneg} pulling it to \eqn{0} and
    repeat the process until no new penalties are introduced.
  }
  \item{x0}{
    A vector specifying the grid where the function to be differentiated
    is measured. Must be sorted.
  }
}
\details{
  \subsection{Whittaker smoothing}{
    Whittaker smoothing works by minimising a sum of penalties [1].
    Taking weights into account, we get: \deqn3{
      \min_{\mathbf{\hat z}} \:
      (\mathbf{\hat z} - \mathbf{z})^\top
        \mathrm{diag}(\mathbf w)
        (\mathbf{\hat z} - \mathbf{z})
      + \lambda | \mathbf D \mathbf{\hat z} |^2
    }{min over z.hat { sum((z.hat - z)^2 * w) + \lambda |D z.hat|^2 }}{
      <span style="font-style:normal;">min</span>
        <sub><b>z&#x0302;</b></sub>
      (<b>z&#x0302;</b> - <b>z</b>)<sup>
        <span style="font-style:normal;">T</span>
      </sup>
      <span style="font-style:normal;">diag</span>(<b>w</b>)
      (<b>z&#x0302;</b> - <b>z</b>)
      + &lambda; |<b>D</b> <b>z&#x0302;</b>|<sup>2</sup>
    }

    By writing down the derivatives over
    \eqn3{\mathbf{\hat z}}{z.hat}{<b>z&#x0302;</b>} and equating them to
    \eqn{0}, we get the normal equation: \deqn3{
      (
        \mathrm{diag}(\mathbf w) +
        \lambda \mathbf{D}^\top \mathbf{D}
      ) \mathbf{\hat z} = \mathrm{diag}(\mathbf w) \mathbf z
    }{(diag(w) + \lambda crossprod(D)) z.hat = diag(w) z}{
      (<span style="font-style:normal;">diag</span>(<b>w</b>)
        + &lambda; <b>D</b><sup>
          <span style="font-style:normal;">T</span>
        </sup><b>D</b>) <b>z&#x0302;</b> =
      <span style="font-style:normal;">diag</span>(<b>w</b>)
      <b>z</b>
    }

    The equation is then solved using \code{\link{solve}}.

    Equations above deal with vectors, not matrices. We can't just
    smooth the matrix row-by-row (or column-by-column): that wouldn't be
    any better than \code{\link{feemscatter}(method = 'pchip')}. One way
    to generalise this to matrices would be [3], but it only works for
    uniform grids. Instead, we unfold the \code{z} matrix into a vector,
    then strategically position the coefficients in the difference
    matrix to obtain both row-wise and column-wise differences in the
    resulting vector.
  }
  \subsection{Finite difference approximation}{
    How to differentiate a function tabulated on a fixed, potentially
    nonuniform grid before you even know its values? Use its Taylor
    series.

    First derivative is special because it's possible to use central
    differences and get a second-order accurate result, even on a
    non-uniform grid, by carefully choosing the points where the
    derivative is calculated. Let \eqn{x + \frac{h}{2}}{x + h/2} and
    \eqn{x - \frac{h}{2}}{x - h/2} be a pair of adjacent points from the
    grid. Here's the Taylor series expansion for \eqn{f} around \eqn{x},
    with the Lagrange form of the reminder:

    \ifelse{latex}{}{\strong{Honestly, just read the LaTeX version.}}

    \deqn{
      f \left(x + \frac{h}{2}\right) =
        f(x) + \frac{h}{2} f'(x) + \frac{h^2}{8} f''(x) +
        \frac{h^3}{48} f'''(\zeta)
    }{f(x + h/2) = f(x) + h/2 f'(x) + h^2/8 f''(x) + h^3/48 f'''(\zeta)}

    \deqn{
      f \left(x - \frac{h}{2}\right) =
        f(x) - \frac{h}{2} f'(x) + \frac{h^2}{8} f''(x) -
        \frac{h^3}{48} f'''(\eta)
    }{f(x - h/2) = f(x) - h/2 f'(x) + h^2/8 f''(x) - h^3/48 f'''(\eta)}

    \deqn{
      f \left(x + \frac{h}{2}\right) -
        f \left(x - \frac{h}{2}\right) = h f'(x) +
        \frac{h^3}{48} \left(f'''(\zeta) + f'''(\eta)\right)
    }{f(x + h/2) - f(x - h/2) = h f'(x) + h^3/48 (f'''(\zeta) + f'''(\eta))}

    \deqn{
      f'(x) = \frac{
        f\left(x + \frac{h}{2}\right) - f\left(x - \frac{h}{2}\right)
      }{h} - \frac{h^2}{48} \left(f'''(\zeta) + f'''(\eta)\right)
    }{f'(x) = (f(x + h/2) - f(x - h/2))/h - h^2/48 (f'''(\zeta) + f'''(\eta))}

    \deqn{
      |\delta f'(x)| \le \max_{
        \xi \in [x - \frac{h}{2}, x + \frac{h}{2}]
      } \frac{h^2}{24} f'''(\xi)
    }

    \if{latex,html}{
      \Sexpr[stage=build]{
        if (!dir.exists('man/figures')) dir.create('man/figures')
        pdf('man/figures/whittaker2points.pdf', 3, 1.3, pointsize = 9)
        dev.control(displaylist = 'enable')
        x <- cumsum(c(0, 2, 2, 5, 5))
        par(mar = c(2.5, 1, 0, 1) + .1)
        plot(0, 0, 'n', range(x), c(0, .8), axes = FALSE, xlab = '', ylab = '')
        axis(1, x[c(1, 3, 5)], expression(xi[1], xi[2], xi[3]))
        axis(
          1, x[c(2, 4)], expression(x[1], x[2]),
          tck = .07, padj = -4, lwd = NA, lwd.ticks = 1
        )
        axis(
          1, (x[2]+x[4])/2, expression(frac(x[1]+x[2],2)),
          tck = .07, padj = -1.2, lwd = NA, lwd.ticks = 1
        )
        arrows(c(x[1], x[3]), c(.5, .7), c(x[3], x[5]), length = .1, code = 3)
        text(c(x[2], x[4]), c(.5, .7), expression(h[1], h[2]), pos = c(3, 1))
        dev.print(
          svg, 'man/figures/whittaker2points.svg',
          width = 3, height = 1.3, pointsize = 9
        )
        dev.off()
        ''
      } % FIXME: with results=hide we get <pre>newlines</pre> in HTML
      \center{
        \if{latex}{
          \figure{whittaker2points.pdf}{options: width=2.5in}
        }
        \if{html}{\figure{whittaker2points.svg}}
      }
    }
    \if{text}{\preformatted{
           h1  <--------------->
       <------->       h2
             (x1+x2)/2
           x1    |     x2
       |---|---|-|-----|-------|
      xi1     xi2             xi3
    }}

    Suppose the three grid points \eqn{\xi_1 = x_1 - \frac{h_1}{2}},
    \eqn{\xi_2 = x_1 + \frac{h_1}{2} = x_2 - \frac{h_2}{2}},
    \eqn{\xi_3 = x_2 + \frac{h_2}{2}} are adjacent on the grid, and we
    know the \eqn{f'} estimates in \eqn{x_1} and \eqn{x_2}.  On the one
    hand, Taylor series expansion for \eqn{f'(x)} around
    \eqn{\frac{x_1 + x_2}{2}} gives:

    \deqn{
      f'(x_2) = f'\left(\frac{x_1 + x_2}{2}\right)
        + f''\left(\frac{x_1 + x_2}{2}\right)\frac{x_2 - x_1}{2}
        + f'''\left(\frac{x_1 + x_2}{2}\right)\frac{(x_2 - x_1)^2}{8}
        + f''''(\zeta)\frac{(x_2 - x_1)^3}{48}
    }
    \deqn{
      f'(x_1) = f'\left(\frac{x_1 + x_2}{2}\right)
        - f''\left(\frac{x_1 + x_2}{2}\right)\frac{x_2 - x_1}{2}
        + f'''\left(\frac{x_1 + x_2}{2}\right)\frac{(x_2 - x_1)^2}{8}
        - f''''(\eta)\frac{(x_2 - x_1)^3}{48}
    }
    \deqn{
      f''\left(\frac{x_1 + x_2}{2}\right) =
        \frac{f'(x_2) - f'(x_1)}{x_2 - x_1} -
        \frac{(x_2 - x_1)^2}{48}(f''''(\zeta) + f''''(\eta))
    }
    \deqn{
      |\delta f''(x)| \le \max_{
        \xi \in [x_1, x_2]
      } \frac{(x_2 - x_1)^2}{24} f''''(\xi)
    }

    On the other hand, if we substitute the estimations of \eqn{f'(x)}
    from above, we get:

    \deqn{
      f''\left(\frac{x_1 + x_2}{2}\right) = \frac{
        h_1 f(\xi_3) - h_1 f(\xi_2)
        - h_2 f(\xi_2) + h_2 f(\xi_1)
      }{h_1 h_2 (x_2 - x_1)} - \frac{
        h_2^2 f'''(\zeta_2)
        - h_1^2 f'''(\zeta_1)
      }{24(x_2 - x_1)}
      - \frac{(x_2 - x_1)^2}{24} f''''(\eta)
    }

    This is why we can't just keep using central differences and get
    second order accurate results.

    What are the general methods of finding the coefficients for the
    \eqn3{\mathbf D}{D}{<b>D</b>} matrix? Start with a system of Taylor
    series expansions for every grid point:

    \deqn{
      f(x_i) = \sum_{k=0}^{n-1} f^{(k)}(x) \frac{(x_i - x)^k}{k!}
      + f^{(n)}(\xi) \frac{(x_i - x)^{n}}{(n)!}
      \; \forall i = 1 \dots n
    }

    We can solve this system for coefficients \eqn{c_i} giving the
    needed the desired \eqn{l}-th derivative estimate with highest
    accuracy \eqn{p} possible [4]:

    \deqn{
      \sum_{i = 1}^n c_i f(x_i) = f^{(l)}(x) + o(h^p)
    }

    Substituting the approximations for \eqn{f(x_i)} into the equation,
    we get the following condition for the multiplier in front of each
    \eqn{f^{(k)}(x)}:

    \deqn{
      \frac{1}{k!} \sum_{i = 1}^n c_i (x_i - x)^k = \mathbf{1}_{k = l}
      \; \forall k = 0 \dots n-1
    }

    In the matrix form this becomes a Vandermonde system:

    \deqn{V_{k,i} = \frac{(x_i - x)^k}{k!}}
    \deqn{b_k = \mathbf{1}_{k = l}}
    \deqn{\mathbf c = \mathbf{V}^{-1} \mathbf b}

    Unfortunately, this system becomes ill-conditioned for
    \dQuote{large} numbers of points. (Experiment shows noticeably
    growing \eqn{c_i} even for third derivative from \eqn{10} points and
    no solution for \eqn{32} points on a uniform grid.) Fornberg [5]
    suggests a more numerically stable procedure, but it still breaks
    down for \eqn{1000} points.

    It is important to note that the performance of the method depends
    on the matrix \eqn{\mathbf D} being sparse. While the methods
    described above could give more accurate results, they do so at the
    cost of providing nonzero weights for a lot of points, and the
    weights get larger as the number of points increases. Therefore,
    with the knowledge that difference orders above \eqn{3} are used
    very rarely and the interest in simplicity and performance, we'll
    minimise the number of coefficients and their values by solving the
    Vandermonde system for the minimally accurate derivative
    estimations, taking exactly \eqn{k + 1} points for \eqn{k}-th
    derivative. For ever

    What is the error of such estimation? Substituting the Lagrange form
    of the remainder into \eqn{\mathbf{c}^\top f(\mathbf x)}, we get:

    \deqn{
      \sum_{i = 1}^n c_i f(x_i) = f^{(n-1)}(x) +
        \sum_{i = 1}^n c_i f^{(n)}(\xi_i) \frac{(x_i - x)^n}{n!},
      \; \xi_i \in [ x_i, x ]
    }

    Our choice of \eqn{x} (middle point for odd \eqn{n}, average of
    middle points for even \eqn{n}) lets us shave off one term from the
    sum above for odd \eqn{n} and get second order accurate results for
    \eqn{n = 2}, but other than that, the method is \eqn{n}-th order
    accurate.
  }
}
\value{
  \item{\code{whittaker2}}{
    A matrix of the same shape as \code{z}, with values smoothed,
    interpolated, or baseline estimated.
  }
  \item{\code{diffmat}}{
    A difference matrix \eqn3{\mathbf D}{D}{<b>D</b>} that, when
    multiplied by \eqn3{\mathrm{vec}(\mathbf{Z})}{vec(Z)}{
      <span style="font-style: normal;">vec</span>(<b>Z</b>)},
    returns a vector of estimated derivatives of \eqn3{\mathbf
    Z}{Z}{<b>Z</b>} of given order by \eqn{x} and \eqn{y}.
  }
  \item{\code{Vandermonde}}{
    A vector \eqn{\mathbf c} of length \eqn{n} such that
    \eqn{\mathbf{c}^\top f(\mathbf{x_0}) \approx f^{(n-1)}(x)}. The
    \eqn{x} is either the median point of \code{x0} or the average of
    two points in the middle, depending on whether \code{length(x0)} is
    odd or even.
  }
}
\references{\enumerate{
  \item \bibentry{'Article',
    title = 'A Perfect Smoother',
    volume = 75,
    doi = '10.1021/ac034173t',
    number = 14,
    journal = 'Analytical Chemistry',
    author = person(c('Paul', 'H.', 'C.'), 'Eilers'),
    year = 2003,
    pages = '3631-3636',
  }

  \item \bibentry{'Article',
    title = 'Parametric Time Warping',
    volume = 76,
    doi = '10.1021/ac034800e',
    number = 2,
    journal = 'Analytical Chemistry',
    author = person(c('Paul', 'H.', 'C.'), 'Eilers'),
    year = 2004,
    pages = '404-411',
  }

  \item \bibentry{'Article',
    title = 'Robust smoothing of gridded data in one and higher dimensions with missing values',
    volume = 54,
    doi = '10.1016/j.csda.2009.09.020',
    number = 4,
    journal = 'Computational Statistics & Data Analysis',
    author = person('Damien', 'Garcia'),
    year = 2010,
    pages = '1167-1178'
  }

  \item \bibentry{'InBook',
    author = person(c('Randall', 'J.'), 'LeVeque'),
    booktitle = 'Finite Difference Methods for Ordinary and Partial Differential Equations',
    chapter = 1,
    pages = '3-11',
    publisher = 'Society for Industrial and Applied Mathematics (SIAM)',
    title = 'Finite Difference Approximations',
    url = 'https://faculty.washington.edu/rjl/fdmbook/',
    year = 2007
  }

  \item \bibentry{'Article',
    title = 'Generation of finite difference formulas on arbitrarily spaced grids',
    volume = 51,
    doi = '10.1090/S0025-5718-1988-0935077-0',
    number = 184,
    journal = 'Mathematics of Computation',
    author = person('Bengt', 'Fornberg'),
    year = 1988,
    pages = '699-706'
  }
}}
\seealso{\code{\link{feemscatter}}}
\examples{
  data(feems)
  z <- feemscatter(feems$a, rep(25, 4), 'omit')
  str(albatross:::whittaker2(
    attr(z, 'emission'), attr(z, 'excitation'), z,
    c(1, 1e-3), 1:2, logscale = NA, nonneg = 1
  ))
}
\keyword{internal}
