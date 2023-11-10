\name{whittaker2}
\alias{whittaker2}
\alias{diffmat}
\alias{vandermonde}
\title{
  Implementation notes for Whittaker smoothing and interpolation of surfaces
}
\description{
  Smooth \cite{Eilers2003} or estimate the baseline
  \cite{Eilers2004-PTW} of a surface measured on an arbitrary grid by
  minimising a sum of penalties. Combined difference orders and two
  different methods of preventing negative values in the output are
  supported.

  \emph{This is not a public interface.} Subject to change without
  further notice. Please do not call from outside \CRANpkg{albatross}.
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
    If not missing, use the asymmetric penalty method
    \cite{Eilers2004-PTW} to estimate the baseline by penalising the
    differences with weight \eqn{p} if
    \eqn3{\hat{z} < z}{z.hat < z}{z&#x0302; &lt; z} and \eqn{1 - p}
    otherwise. Typically, values around
    \eqn3{10^{-3}}{1e-3}{10<sup>-3</sup>} are used.
  }
  \item{logscale}{
    If not \code{NA}, rescale \code{z} between \code{logscale} and
    \eqn{1} and take its logarithm before smoothing, then undo it.

    Such transformation prevents the resulting values from getting lower
    than \eqn3{
      \mathrm{min}(x) - (\mathrm{max}(x) - \mathrm{min}(x))
      \frac{\mathtt{logscale}}{1 - \mathtt{logscale}}
    }{min(x) - diff(range(x)) * logscale / (1 - logscale)}{
      <span style="font-style: normal;">min</span>(x) -
      (<span style="font-style: normal;">max</span>(x) -
        <span style="font-style: normal;">min</span>(x))
      <sup><code>logscale</code></sup>/<sub>1 -
      <code>logscale</code></sub>},
    which is approximately \eqn3{
      -\mathtt{logscale} \cdot \mathrm{max}(x)
    }{-logscale * max(x)}{
      -<code>logscale</code> &middot;
      <span style="font-style: normal;">max</span>(x)
    } if \code{logscale} and \eqn3{\mathrm{min}(x)}{min(x)}{
      <span style="font-style: normal;">min</span>(x)
    } are both close to \eqn{0}.

    A typical value would be \eqn3{10^{-4}}{10^-4}{10<sup>-4</sup>}.
    Disabled by default because it may damage the surface shape in the
    process.
  }
  \item{nonneg}{
    If not \eqn{0}, for every resulting negative value in the
    interpolated surface, add a penalty of \eqn3{
      \mathtt{nonneg} \cdot \sum_i \mathbf{1}_{\hat{z}_i < 0}
      \, \hat{z}_i^2
    }{nonneg * sum(z.hat[z.hat < 0]^2)}{
      <code>nonneg</code> &middot; &#x2211;<sub>i</sub>
      <b>1</b><sub>z&#x0302;<sub>i</sub> &lt; 0</sub>
      z&#x0302;<sub>i</sub><sup>2</sup>
    } to pull it towards \eqn{0} and repeat the process until no new
    penalties are introduced.
  }
  \item{x0}{
    A vector specifying the grid where the function to be differentiated
    is measured. Must be sorted.
  }
}
\details{
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

    \deqn3{
      f \left(x + \frac{h}{2}\right) =
        f(x) + \frac{h}{2} f'(x) + \frac{h^2}{8} f''(x) +
        \frac{h^3}{48} f'''(\zeta)
    }{f(x + h/2) = f(x) + h/2 f'(x) + h^2/8 f''(x) + h^3/48 f'''(\zeta)}{
      f (x + <sup>h</sup>/<sub>2</sub>) =
        f(x) + <sup>h</sup>/<sub>2</sub> f'(x)
        + <sup>h<sup>2</sup></sup>/<sub>8</sub> f''(x)
        + <sup>h<sup>3</sup></sup>/<sub>48</sub> f'''(&zeta;)
    }

    \deqn3{
      f \left(x - \frac{h}{2}\right) =
        f(x) - \frac{h}{2} f'(x) + \frac{h^2}{8} f''(x) -
        \frac{h^3}{48} f'''(\eta)
    }{f(x - h/2) = f(x) - h/2 f'(x) + h^2/8 f''(x) - h^3/48 f'''(\eta)}{
      f (x - <sup>h</sup>/<sub>2</sub>) =
        f(x) - <sup>h</sup>/<sub>2</sub> f'(x)
        + <sup>h<sup>2</sup></sup>/<sub>8</sub> f''(x)
        - <sup>h<sup>3</sup></sup>/<sub>48</sub> f'''(&eta;)
    }

    \deqn3{
      f \left(x + \frac{h}{2}\right) -
        f \left(x - \frac{h}{2}\right) = h f'(x) +
        \frac{h^3}{48} \left(f'''(\zeta) + f'''(\eta)\right)
    }{f(x + h/2) - f(x - h/2) = h f'(x) + h^3/48 (f'''(\zeta) + f'''(\eta))}{
      f(x + <sup>h</sup>/<sub>2</sub>) - f(x - <sup>h</sup>/<sub>2</sub>) =
      h f'(x) + <sup>h<sup>3</sup></sup>/<sub>48</sub> (f'''(&zeta;) + f'''(&eta;))
    }

    \deqn3{
      f'(x) = \frac{
        f\left(x + \frac{h}{2}\right) - f\left(x - \frac{h}{2}\right)
      }{h} - \frac{h^2}{48} \left(f'''(\zeta) + f'''(\eta)\right)
    }{f'(x) = (f(x + h/2) - f(x - h/2))/h - h^2/48 (f'''(\zeta) + f'''(\eta))}{
      f'(x) = <sup>1</sup>/<sub>h</sub> (
        f(x + <sup>h</sup>/<sub>2</sub>) - f(x - <sup>h</sup>/<sub>2</sub>)
      )
      - <sup>h<sup>2</sup></sup>/<sub>48</sub> (f'''(&zeta;) + f'''(&eta;))
    }

    \deqn3{
      |\delta f'(x)| \le \max_{
        \xi \in [x - \frac{h}{2}, x + \frac{h}{2}]
      } \frac{h^2}{24} f'''(\xi)
    }{|\delta f'(x)| \le h^2/24 max f'''(\xi) over \xi in [x - h/2; x + h/2]}{
      |&delta; f'(x)| &le;
        <sup>h<sup>2</sup></sup>/<sub>24</sub>
        <span style="font-style:normal;">max</span> f'''(&xi;)
        <span style="font-style:normal;">over</span>
          &xi; &in; [x - <sup>h</sup>/<sub>2</sub>; x + <sup>h</sup>/<sub>2</sub>]
    }

    Suppose the three grid points
    \eqn3{\xi_1 = x_1 - \frac{h_1}{2}}{\xi1 = x1 - h1/2}{
      &xi;<sub>1</sub> = x<sub>1</sub> - <sup>h<sub>1</sub></sup>/<sub>2</sub>
    },
    \eqn3{
      \xi_2 = x_1 + \frac{h_1}{2} = x_2 - \frac{h_2}{2}
    }{\xi2 = x1 + h1/2 = x2 - h2/2}{
      &xi;<sub>2</sub> = x<sub>1</sub> + <sup>h<sub>1</sub></sup>/<sub>2</sub>
      = x<sub>2</sub> - <sup>h<sub>2</sub></sup>/<sub>2</sub>
    },
    \eqn3{\xi_3 = x_2 + \frac{h_2}{2}}{\xi3 = x2 + h2/2}{
      &xi;<sub>3</sub> = x<sub>2</sub> + <sup>h<sub>2</sub></sup>/<sub>2</sub>
    } are adjacent on the grid, and we
    know the \eqn{f'} estimates in \eqn3{x_1}{x1}{x<sub>1</sub>} and
    \eqn3{x_2}{x2}{x<sub>2</sub>}:

    \if{latex,html}{
      \makefigure{
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
      }
      \if{latex}{\figure{whittaker2points.pdf}{options: width=2.5in}}
      \if{html}{ \figure{whittaker2points.svg}{options: width="400"}}
    }
    \if{text}{\preformatted{
           h1  <--------------->
       <------->       h2
             (x1+x2)/2
           x1    |     x2
       |---|---|-|-----|-------|
      xi1     xi2             xi3
    }}

    On the one hand, Taylor series
    expansion for \eqn{f'(x)} around
    \eqn3{\frac{x_1 + x_2}{2}}{(x1 + x2)/2}{
      <sup>1</sup>/<sub>2</sub>(x<sub>1</sub> + x<sub>2</sub>)
    } gives:

    \deqn3{
      f'(x_2) = f'\left(\frac{x_1 + x_2}{2}\right)
        + f''\left(\frac{x_1 + x_2}{2}\right)\frac{x_2 - x_1}{2}
        + f'''\left(\frac{x_1 + x_2}{2}\right)\frac{(x_2 - x_1)^2}{8}
        + f''''(\zeta)\frac{(x_2 - x_1)^3}{48}
    }{
      f'(x2) =    f'((x1 + x2) / 2)
             +   f''((x1 + x2) / 2) (x2 - x1)   /  2 +
             +  f'''((x1 + x2) / 2) (x2 - x1)^2 /  8 +
             + f''''(zeta)          (x2 - x1)^3 / 48
    }{
      f'(x<sub>2</sub>)
      = f'(<sup>1</sup>/<sub>2</sub> (x<sub>1</sub> + x<sub>2</sub>))
      + f''(<sup>1</sup>/<sub>2</sub> (x<sub>1</sub> + x<sub>2</sub>))
        (x<sub>2</sub> - x<sub>1</sub>) <sup>1</sup>/<sub>2</sub>
      + f'''(<sup>1</sup>/<sub>2</sub> (x<sub>1</sub> + x<sub>2</sub>))
        (x<sub>2</sub> - x<sub>1</sub>)<sup>2</sup> <sup>1</sup>/<sub>8</sub>
      + f''''(&zeta;)
        (x<sub>2</sub> - x<sub>1</sub>)<sup>3</sup> <sup>1</sup>/<sub>48</sub>
    }
    \deqn3{
      f'(x_1) = f'\left(\frac{x_1 + x_2}{2}\right)
        - f''\left(\frac{x_1 + x_2}{2}\right)\frac{x_2 - x_1}{2}
        + f'''\left(\frac{x_1 + x_2}{2}\right)\frac{(x_2 - x_1)^2}{8}
        - f''''(\eta)\frac{(x_2 - x_1)^3}{48}
    }{
      f'(x1) =    f'((x1 + x2) / 2)
             -   f''((x1 + x2) / 2) (x2 - x1)   /  2 +
             +  f'''((x1 + x2) / 2) (x2 - x1)^2 /  8 +
             - f''''(eta)           (x2 - x1)^3 / 48
    }{
      f'(x<sub>1</sub>)
      = f'(<sup>1</sup>/<sub>2</sub> (x<sub>1</sub> + x<sub>2</sub>))
      - f''(<sup>1</sup>/<sub>2</sub> (x<sub>1</sub> + x<sub>2</sub>))
        (x<sub>2</sub> - x<sub>1</sub>) <sup>1</sup>/<sub>2</sub>
      + f'''(<sup>1</sup>/<sub>2</sub> (x<sub>1</sub> + x<sub>2</sub>))
        (x<sub>2</sub> - x<sub>1</sub>)<sup>2</sup> <sup>1</sup>/<sub>8</sub>
      - f''''(&eta;)
        (x<sub>2</sub> - x<sub>1</sub>)<sup>3</sup> <sup>1</sup>/<sub>48</sub>
    }
    \deqn3{
      f''\left(\frac{x_1 + x_2}{2}\right) =
        \frac{f'(x_2) - f'(x_1)}{x_2 - x_1} -
        \frac{(x_2 - x_1)^2}{48}(f''''(\zeta) + f''''(\eta))
    }{
      f''((x1 + x2)/2) = (f'(x2) - f'(x1)) / (x2 - x1) -
         - (x2 - x1)^2 / 48 (f''''(zeta) + f''''(eta))
    }{
      f''(<sup>1</sup>/<sub>2</sub> (x<sub>1</sub> + x<sub>2</sub>)) =
      <sup>f'(x<sub>2</sub>) -
      f'(x<sub>1</sub>)</sup>/<sub>x<sub>2</sub> - x<sub>1</sub></sub>
      - <sup>1</sup>/<sub>48</sub>
        (x<sub>2</sub> - x<sub>1</sub>)<sup>2</sup>
        (f''''(&zeta;) + f''''(&eta;))
    }
    \deqn3{
      |\delta f''(x)| \le \max_{
        \xi \in [x_1, x_2]
      } \frac{(x_2 - x_1)^2}{24} f''''(\xi)
    }{|delta f''(x)| \le (x2 - x1)^2 / 24 max f''''(xi) over xi in [x1, x2]}{
      |&delta; f''(x)| &le; <sup>1</sup>/<sub>24</sub>
      (x<sub>2</sub> - x<sub>1</sub>)<sup>2</sup>
      <span style="font-style:normal;">max</span> f''''(&xi;)
      <span style="font-style:normal;">over</span> &xi; &in;
      [x<sub>1</sub>, x<sub>2</sub>]
    }

    On the other hand, if we substitute the estimations of \eqn{f'(x)}
    from above, we get:

    \deqn3{
      f''\left(\frac{x_1 + x_2}{2}\right) = \frac{
        h_1 f(\xi_3) - h_1 f(\xi_2)
        - h_2 f(\xi_2) + h_2 f(\xi_1)
      }{h_1 h_2 (x_2 - x_1)} - \frac{
        h_2^2 f'''(\zeta_2)
        - h_1^2 f'''(\zeta_1)
      }{24(x_2 - x_1)}
      - \frac{(x_2 - x_1)^2}{24} f''''(\eta)
    }{
         (x1 + x2)
      f''(-------) =
         (   2   )

        h1 f(xi3) - h1 f(xi2) - h2 f(xi2) + h2 f(xi1)
        --------------------------------------------- -
                        h1 h2 (x2 - x1)

          h2^2 f'''(xi2) - h1^2 f'''(xi1)   (x2 - x1)^2
        - ------------------------------- - ----------- f''''(eta)
                   24 (x2 - x1)                  24
    }{<!-- here comes the really painful part -->
      f''(<sup>1</sup>/<sub>2</sub>(x<sub>1</sub> + x<sub>2</sub>)) =
        (h<sub>1</sub> f(&xi;<sub>3</sub>)
        - h<sub>1</sub> f(&xi;<sub>2</sub>)
        - h<sub>2</sub> f(&xi;<sub>2</sub>)
        + h<sub>2</sub> f(&xi;<sub>1</sub>))
        (h<sub>1</sub> h<sub>2</sub> (x<sub>2</sub> - x<sub>1</sub>))
        <sup>-1</sup>
      -
        (h<sub>2</sub><sup>2</sup> f'''(&xi;<sub>2</sub>)
        - h<sub>1</sub><sup>2</sup> f'''(&xi;<sub>1</sub>))
        (24 (x<sub>2</sub> - x<sub>1</sub>)) <sup>-1</sup>
      - <sup>1</sup>/<sub>24</sub>
        (x<sub>2</sub> - x<sub>1</sub>)<sup>2</sup>
        f''''(&eta;)
    }

    This is why we can't just keep using central differences and get
    \eqn{n+1}th order accurate results.

    What are the general methods of finding the coefficients for the
    \eqn3{\mathbf D}{D}{<b>D</b>} matrix? Start with a system of Taylor
    series expansions for every grid point:

    \deqn3{
      f(x_i) = \sum_{k=0}^{n-1} f^{(k)}(x) \frac{(x_i - x)^k}{k!}
      + f^{(n)}(\xi) \frac{(x_i - x)^{n}}{n!}
      \; \forall i = 1 \dots n
    }{
                n-1 (k)  (x[i] - x)^k   (n)   (x[i] - x)^n
      f(x[i]) = sum f(x) ------------ + f(xi) ------------, i = 1..n
                k=0           k!                   n!
    }{
      f(x<sub>i</sub>) = &sum;<sub>k=0..n-1</sub>
        f<sup>(k)</sup>(x)
        <sup>(x<sub>i</sub> - x)<sup>k</sup></sup>/<sub>k!</sub>
      + f<sup>(n)</sup>(&xi;)
        <sup>(x<sub>i</sub> - x)<sup>n</sup></sup>/<sub>n!</sub>
    }

    We can solve this system for coefficients
    \eqn3{c_i}{c[i]}{c<sub>i</sub>} giving the desired \eqn{l}-th
    derivative estimate with highest accuracy \eqn{p} possible
    \cite{LeVeque2007-ch1}:

    \deqn3{
      \sum_{i = 1}^n c_i f(x_i) = f^{(l)}(x) + o(h^p)
    }{
      sum over i=1:n { c[i] f(x[i]) } = d^l/dx^l f(x[i]) + o(h^p)
    }{
      &sum;<sub>i=1..n</sub> c<sub>i</sub> f(x<sub>i</sub>) =
      f<sup>(l)</sup>(x) + o(h<sup>p</sup>)
    }

    Substituting the approximations for
    \eqn3{f(x_i)}{f(x[i])}{f(x<sub>i</sub>)} into the equation, we get the
    following condition for the multiplier in front of each
    \eqn3{f^{(k)}(x)}{d^k/dx^k f(x)}{f<sup>(k)</sup>(x)}:

    \deqn3{
      \frac{1}{k!} \sum_{i = 1}^n c_i (x_i - x)^k = \mathbf{1}_{k = l}
      \; \forall k = 0 \dots n-1
    }{1/k! sum over i=1:n { c[i] (x[i] - x)^k } = (k == l), k = 0..n-1}{
      <sup>1</sup>/<sub>k!</sub>
      &sum;<sub>i=1..n</sub>
        c<sub>i</sub> (x<sub>i</sub> - x)<sup>k</sup>
      = <b>1</b><sub>k=l</sub>
      &forall;k = 0..n-1
    }

    In the matrix form, this becomes a Vandermonde system:

    \deqn3{
      V_{k,i} = \frac{(x_i - x)^k}{k!}
    }{V[k,i] = (x[i] - x)^k/k!}{
      V<sub>k,i</sub> = <sup>(x<sub>i</sub> -
        x)<sup>k</sup></sup>/<sub>k!</sub>
    }
    \deqn3{
      b_k = \mathbf{1}_{k = l}
    }{b[k] = (k == l)}{
      b<sub>k</sub> = <b>1</b><sub>k=l</sub>
    }
    \deqn3{
      \mathbf c = \mathbf{V}^{-1} \mathbf b
    }{c = V^-1 \\\%*\\\% b}{
      <b>c</b> = <b>V</b><sup>-1</sup> <b>b</b>
    }

    Unfortunately, this system becomes ill-conditioned for \dQuote{large}
    numbers of points. (Experiment shows noticeably growing
    \eqn3{c_i}{c[i]}{c<sub>i</sub>} even for third derivative from \eqn{10}
    points and no solution for \eqn{32} points on a uniform grid.)
    Fornberg \cite{Fornberg1988} suggests a more numerically stable
    procedure, but it still breaks down for \eqn{1000} points.

    It is important to note that the performance of the method depends
    on the matrix \eqn3{\mathbf D}{D}{<b>D</b>} being sparse. While the methods
    described above could give more accurate results, they do so at the cost of
    providing nonzero weights for a lot of points, and the weights get larger
    as the number of points increases. Therefore, with the knowledge that
    difference orders above \eqn{3} are used very rarely and the interest in
    simplicity and performance, we'll minimise the number of coefficients and
    their values by solving the Vandermonde system for the minimally accurate
    derivative estimations, taking exactly \eqn{k + 1} points for \eqn{k}-th
    derivative.

    What is the error of such estimation? Substituting the Lagrange form
    of the remainder into \eqn3{
      \mathbf{c}^\top f(\mathbf x)
    }{t(c) \\\%*\\\% x}{
      <b>c</b><sup><span style="font-style:normal;">T</span></sup>
      f(<b>x</b>)}, we get:

    \deqn3{
      \sum_{i = 1}^n c_i f(x_i) = f^{(n-1)}(x) +
        \sum_{i = 1}^n c_i f^{(n)}(\xi_i) \frac{(x_i - x)^n}{n!},
      \; \xi_i \in [ x_i, x ]
    }{
       n                 (n-1)   n       (n)      (x[n] - x)^n
      sum c[i] f(x[i]) = f(x) + sum c[i] f(xi[i]) ------------,
      i=1                       i=1                     n!

                                           xi[i] in [x[i], x]
    }{
      &sum;<sub>i=1..n</sub> c<sub>i</sub> f(x<sub>i</sub>) =
      f<sup>(n-1)</sup>(x) +
      &sum;<sub>i=1..n</sub> c<sub>i</sub>
        f<sup>(n)</sup>(&xi;<sub>i</sub>)
        <sup>(x<sub>i</sub> - x)<sup>n</sup></sup>/<sub>n!</sub>,
      &xi;<sub>i</sub> &in; [x<sub>i</sub>, x]
    }

    Our choice of \eqn{x} (middle point for odd \eqn{n}, average of
    middle points for even \eqn{n}) lets us shave off one term from the
    sum above for odd \eqn{n} and get second order accurate results for
    \eqn{n = 2}, but other than that, the method is \eqn{n}-th order
    accurate.
  }
  \subsection{Whittaker smoothing}{
    Whittaker smoothing works by minimising a sum of penalties
    \cite{Eilers2003}. Interpolation can be achieved by setting weights
    \eqn3{\mathbf w}{w}{<b>w</b>} to \eqn{0} for the missing points.

    \deqn3{
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
    \eqn{0}, we get the normal equation:

    \deqn3{
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

    The equation is then solved using \code{\link[Matrix]{solve}}.

    Given a one-dimensional penalty matrix
    \eqn3{\mathbf{D}_d}{D(d)}{<b>D</b><sub>d</sub>} of order \eqn{d}
    obtained by solving a Vandermonde system for every successive group
    of \eqn{d+1} points, we can adapt it for every applicable group of
    points from a fluorescence excitation-emission matrix unfolded into
    a vector \eqn3{\mathbf z = \mathrm{vec}\, \mathbf F}{z = vec(F)}{
      <b>z</b> = <span style="font-style:normal">vec</span> <b>F</b>} by
    taking Kronecker products with unit matrices: \deqn3{
      \mathbf D = \begin{pmatrix}
        \mathbf I \otimes \mathbf{D}_\mathrm{em} \\\\
        \mathbf{D}_\mathrm{ex} \otimes \mathbf I
      \end{pmatrix}
    }{D = rbind(I \\\%x\\\% D.em, D.ex \\\%x\\\% I)}{
      <!-- didn't Chrome reintroduce MathML support recently? -->
      <div style="text-align:center;">
      <b>D</b> = <table style="line-height:1;display:inline-table;"><tr>
      <td>&#9115;<br>&#9116;<br>&#9117;</td>
      <td><table>
        <tr><td><b>I</b> &otimes; <b>D</b><sub>em</sub></td></tr>
        <tr><td><b>D</b><sub>ex</sub> &otimes; <b>I</b></td></tr>
      </table></td>
      <td>&#9118;<br>&#9119;<br>&#9120;</td>
      </tr></table></div>
    }

    Penalties of different orders are concatenated by rows in a similar
    manner (which is equivalent to adding the corresponding
    \eqn3{\mathbf{D}^\top\mathbf D}{crossprod(D)}{
      <b>D</b><sup><span style="font-style:normal;">T</span></sup><b>D</b>
    } matrices together in the normal equation).

    It has been shown in \cite{Eilers2004-density} that a combination of
    first- and second-order penalty (\eqn3{
      2 \lambda \mathbf{D}_1 + \lambda^2 \mathbf{D}_2
    }{2 \lambda D(1) + \lambda^2 D(2)}{2&lambda;<b>D</b><sub>1</sub> +
      &lambda;<sup>2</sup><b>D</b><sub>2</sub>})
    results in a non-negative impulse response, but the resulting peak
    shape may be sub-optimal.
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
    Z}{Z}{<b>Z</b>} of given order by \eqn{x} and \eqn{y} in every available
    point (see below).
  }
  \item{\code{vandermonde}}{
    A vector \eqn3{\mathbf c}{c}{<b>c</b>} of length \eqn{n =}
    \code{length(x0)} such that
    \eqn3{
      \mathbf{c}^\top f(\mathbf{x_0}) \approx f^{(n-1)}(x)
    }{t(c) \\\%*\\\% f(x0) =~ d^(n-1)/dx^(n-1) f(x)}{
      <b>c</b><sup><span style="font-style:normal;">T</span></sup>
      f(<b>x<sub>0</sub></b>) &approx; f<sup>(n-1)</sup>(x)}. The
    \eqn{x} is either the median point of \code{x0} or the average of
    two points in the middle, depending on whether \eqn{n} is odd or even.
  }
}
\references{\bibliography}
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
