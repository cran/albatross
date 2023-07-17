# Implementation of weighted Whittaker smoothing. Currently, only used
# by .feeminterpolate in feemscatter.R. There's also an "asymmetric
# least squares" baseline estimation method, but it's not called from
# anywhere. See whittaker2.Rd for implementation details.

# Return coefficients to mutiply f(x) by to estimate the
# `length(x)-1`-th derivative. Estimate is calculated either at midpoint
# or between the two middle points, depending on whether length(x) is
# odd or even, respectively.
vandermonde <- function(x0) {
	x <- if (length(x0) %% 2 == 1) { # midpoint for odd-length vectors
		x0[(length(x0) + 1) / 2]
	} else mean(x0[length(x0)/2] + 0:1) # between midpoints otherwise
	k <- seq_along(x0) - 1 # derivative orders
	# build the Vandermonde matrix
	V <- outer(k, x0, function(k, x0) (x0 - x)^k / factorial(k))
	# find the coefficients for n-1-th derivative estimate
	solve(V, 1 * (k == length(x0) - 1))
}

# x: rows of the grid, y: columns of the grid, d: integer scalar >= 1.
# returns the matrix that, when multiplied with as.vector(z), gives
# a vector of d-th order derivatives of z by x followed by a vector
# of d-th order derivatives of z by y.
diffmat <- function(x, y, d) {
	# obtain a mapping from [i,j] to index in as.vector(z)
	idx <- matrix(seq_len(length(x) * length(y)), nrow = length(x))
	Dx <- sparseMatrix(
		# Differentiating each column gives us less points than its
		# length. The difference is exactly d.
		i = rep(seq_len((length(x) - d) * length(y)), each = d + 1),
		j = as.vector(vapply(
			1:length(y), # for every column...
			function(j) as.vector(vapply(
				1:(length(x) - d), # ...move a window of d over it
				function(i) idx[i + 0:d, j],
				integer(d+1)
			)),
			integer((d+1) * (length(x) - d))
		)),
		x = rep( # for every column...
			as.vector(vapply(
				1:(length(x) - d), # ...move a window of d over it
				function(i) vandermonde(x[i + 0:d]),
				numeric(d+1)
			)),
			times = length(y)
		)
	)
	Dy <- sparseMatrix(
		# Similarly, we differentiate each row, taking d points off each
		i = rep(seq_len((length(y) - d) * length(x)), each = d + 1),
		j = as.vector(vapply(
			1:length(x), # for every row...
			function(i) as.vector(vapply(
				1:(length(y) - d), # ...move a window of d over it
				function(j) idx[i, j + 0:d],
				integer(d+1)
			)),
			integer((d+1) * (length(y) - d))
		)),
		x = rep( # for every row...
			as.vector(vapply(
				1:(length(y) - d), # ...move a window of d over it
				function(j) vandermonde(y[j + 0:d]),
				numeric(d+1)
			)),
			times = length(x)
		)
	)
	# D must be ready to be multiplied by as.vector(z)
	rbind(Dx, Dy)
}

# x must correspond to rows, y to columns of z
# x, y must not repeat
# lambdas correspond to difference orders (could be multiple of them)
# baseline estimation is performed if p is specified, only smoothing otherwise
# any NAs in z are imputed
# for interpolation, [Eilers, 1987] (Graphics Gems 4)
# recommends D = 1st + 2nd differences
# difference matrix, premultiplied with penalty weights
whittaker2 <- function(x, y, z, lambda, d, p, logscale, nonneg) {
	stopifnot(
		!anyDuplicated(x), length(x) > max(d),
		!anyDuplicated(y), length(y) > max(d)
	)
	if (is.unsorted(x) || is.unsorted(y)) {
		# must be sorted because we use diff() in diffmat()
		perm.x <- order(x)
		perm.y <- order(y)
		return(
			Recall(
				x[perm.x], y[perm.y], z[perm.x, perm.y],
				lambda, d, p, logscale, nonneg
			)[
				# order(perm) gives an inverse permutation
				order(perm.x), order(perm.y)
			]
		)
	}

	if (!is.na(logscale)) {
		oldmin <- min(z, na.rm = TRUE)
		oldrange <- diff(range(z, na.rm = TRUE))
		z <- log(logscale + (1 - logscale) * (z - oldmin) / oldrange)
	}

	# start with equal weights but automatically impute missing values
	na <- is.na(z)
	w <- as.numeric(!na)
	z[na] <- 0 # the value shouldn't really matter but has to be defined

	# A Perfect Smoother. Paul H. C. Eilers, Analytical Chemistry,
	# 2003 75 (14), 3631-3636. doi:10.1021/ac034173t
	# let:
	# R = residuals: sum (w * (z - y)^2) = |(y-z)" diag(w) (y-z)|^2
	# S = smoothness: sum (grad z)^2 = |D z|^2
	# minimize Q = R + lambda * S over z:
	#  assume partial derivatives to be 0
	#  => (diag(w) + lambda D" D) z = diag(w) y
	#  solve for z
	lambdaDsq <- crossprod(Reduce(rbind,
		Map(
			function(lambda, d) sqrt(lambda) * diffmat(x, y, d),
			lambda, d
		)
	))

	# optional iterative penalty for negative values
	v <- rep(0, length(z))

	repeat {
		z.hat <- as.vector(
			solve(Diagonal(x = w) + lambdaDsq + Diagonal(x = v), w * as.vector(z))
		)
		if (!missing(nonneg)) {
			# pull negative results to 0 on next iteration
			v.prev <- v
			v[z.hat < 0] <- nonneg
		}
		# Parametric Time Warping. Paul H. C. Eilers, Analytical Chemistry,
		# 2004 76 (2), 404-411. doi:10.1021/ac034800e;
		# Appendix: Asymmetric Least-Squares Baseline Estimation
		if (!missing(p)) {
			w.prev <- w
			w[!na] <- ifelse(as.vector(z)[!na] > z.hat[!na], p, 1 - p)
		}
		if (
			(missing(p) || all(w.prev == w)) && # done estimating baseline
			(missing(nonneg) || all(v.prev == v)) # no new negative values
		) break
	}

	if (!is.na(logscale))
		z.hat <- oldmin + oldrange * (exp(z.hat) - logscale) / (1 - logscale)

	z[] <- z.hat
	z
}
