# NB: different paths taken for even-length and odd-length x0
stopifnot(
	all.equal(albatross:::vandermonde(c(1, 2) / 2), c(-1, 1) * 2),
	all.equal(albatross:::vandermonde(c(1, 2, 3) * 2), c(1, -2, 1) / 2^2)
)

# just some nontrivial grid
x <- c(2, 3, 5, 8)
y <- c(10, 15, 35, 40)
# just some function with nonzero 1st and 2nd derivatives
e <- expression(2 * x^2 + 4 * y^3 + 5 * x * y)
ex <- D(e, 'x')
exx <- D(ex, 'x')
ey <- D(e, 'y')
eyy <- D(ey, 'y')
eyyy <- D(eyy, 'y')
# other non-mixed derivatives are 0

Z <- outer(x, y, function(x, y) eval(e))

m <- function(x) (x[-1] + x[-length(x)])/2

dZdx.true <- outer(m(x), y, function(x, y) eval(ex))
dZdy.true <- outer(x, m(y), function(x, y) eval(ey))

dZ <- as.vector(albatross:::diffmat(x, y, 1) %*% as.vector(Z))
dZdx.est <- matrix(dZ[1:((length(x)-1) * length(y))], ncol = length(y))
# NB: note the order in which dZ/dy is returned
dZdy.est <- matrix(dZ[-(1:((length(x)-1) * length(y)))], nrow = length(x), byrow = TRUE)

# For central difference approximation, Taylor error estimate in Lagrange form
# would be (x2-x1)^2/24 max [x1, x2] d^3f/dx^3
dZdy.err <- outer(
	seq_along(x), seq_along(y[-1]), Vectorize(function(ix, iy, x, y) {
		x <- x[ix]
		y <- y[iy + 0:1]
		abs(diff(y))^2 / 24 * max(abs(eval(eyyy)))
	}, c('ix', 'iy')), x, y
)

check.err <- function(true, est, err = 0, thr = sqrt(.Machine$double.eps)) {
	if (all(abs(true - est) - err <= thr)) return(invisible())

	message(tn <- deparse(substitute(true)))
	print(true)
	message(en <- deparse(substitute(est)))
	print(est)
	message(en, ' - ', tn)
	print(est - true)
	message(deparse(substitute(err)))
	print(err)
	stop('Error not within estimate')
}

# central difference should be exact for x (exxx = 0)
check.err(dZdx.true, dZdx.est)
# residuals for y should be within Taylor estimated error
check.err(dZdy.true, dZdy.est, dZdy.err)

d2Z <- as.vector(albatross:::diffmat(x, y, 2) %*% as.vector(Z))
d2Zdx2.est <- matrix(d2Z[1:((length(x)-2) * length(y))], ncol = length(y))
d2Zdy2.est <- matrix(d2Z[-(1:((length(x)-2) * length(y)))], nrow = length(x), byrow = TRUE)

d2Zdx2.true <- outer(x[3:length(x) - 1], y, function(x, y) rep(eval(exx), len = length(x)))
d2Zdy2.true <- outer(x, y[3:length(y) - 1], function(x, y) eval(eyy))

# again, exact for x
check.err(d2Zdx2.true, d2Zdx2.est)

# f(y) = sum(f^(k)(y0) (y - y0)^k / k!) + f^(3)(yy) (y - y0)^3 / 3!
# sum(vandermonde(y) * f(y)) =
#  = f''(y0) + sum( vandermonde(y) * f^(3)(yy) * (y - y0)^3 / 3! )
# here yy doesn't matter because eyyy is last, constant derivative
d2Zdy2.err <- outer(
	seq_along(x), seq_along(y[-(1:2)]), Vectorize(function(ix, iy, x, y) {
		x <- x[ix]
		y <- y[iy + 0:2]
		y0 <- y[2]
		abs(sum(
			albatross:::vandermonde(y) * eval(eyyy) * (y - y0)^3 / 6
		))
	}, c('ix', 'iy')), x, y
)
check.err(d2Zdy2.true, d2Zdy2.est, d2Zdy2.err)

# third derivative should be exact (0 in case of x, constant in case of y)
d3Z <- as.vector(albatross:::diffmat(x, y, 3) %*% as.vector(Z))
d3Zdx3.est <- matrix(d3Z[1:((length(x)-3) * length(y))], ncol = length(y))
d3Zdy3.est <- matrix(d3Z[-(1:((length(x)-3) * length(y)))], nrow = length(x), byrow = TRUE)
check.err(0, d3Zdx3.est)
check.err(eval(eyyy), d3Zdy3.est)

# all parameters must be correctly passed when recalling self
Zneg <- Z - 2*min(Z)
stopifnot(all.equal(
	albatross:::whittaker2(x, y, Zneg, 1e3, 3, .5, 1e-3, 2),
	albatross:::whittaker2(
		rev(x), y, Zneg[nrow(Zneg):1,],
		1e3, 3, .5, 1e-3, 2
	)[nrow(Zneg):1,]
))
