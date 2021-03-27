# x - rows of the grid, y - columns of the grid, d - integer scalar >= 1.
# returns the matrix that, when multiplied with as.vector(z), gives
# a vector of d-th order derivatives of z by x followed by a vector
# of d-th order derivatives of z by y.
diffmat <- function(x, y, d) {
	# obtain a mapping from [i,j] to index in as.vector(z)
	idx.x <- idx.y <- matrix(seq_len(length(x) * length(y)), nrow = length(x))
	# establish the invariant that Dx, Dy are matrices, build them iteratively
	Dx <- Dy <- Diagonal(length(x) * length(y))
	for (i in seq_len(d)) {
		r.dx <- .5 / diff(x)
		r.dy <- .5 / diff(y)
		# using the matrix of indices, build a single-order difference
		# matrix calculating dz/dx
		Dx <- sparseMatrix(
			# rows of Dx: enough for z minus one row
			i = rep(seq_len(length(r.dx) * ncol(idx.x)), times = 2),
			# columns of D: choose by index from as.vector(z)
			j = c(idx.x[-1,], idx.x[-nrow(idx.x),]),
			# values of Dx: +/- 1 / dx
			# NB: row numbers change the fastest, so we go over dx
			# as element numbers increase
			x = c(
				rep(+r.dx, times = ncol(idx.x)),
				rep(-r.dx, times = ncol(idx.x))
			)
		) %*% Dx # multiply it by D to make it higher order

		# build a matrix calculating dz/dy
		Dy <- sparseMatrix(
			# rows of Dy: enough for z minus one column
			i = rep(seq_len(length(r.dy) * nrow(idx.y)), times = 2),
			# columns of Dy: choose by index from as.vector(z)
			j = c(idx.y[,-1], idx.y[,-ncol(idx.y)]),
			# values of Dy: +/- 1 / dy
			# NB: column number changes the slowest, so we go over same dy
			# multiple times as element numbers increase
			x = c(
				rep(+r.dy, each = nrow(idx.y)),
				rep(-r.dy, each = nrow(idx.y))
			)
		) %*% Dy # increases the order of Dy
		# we use the central difference formula, so the resulting
		# derivative values lie in the middle between each pair
		# of points, losing a row and a column each time
		idx.x <- matrix(seq_along(idx.x[-1,]), nrow = length(x) - 1)
		idx.y <- matrix(seq_along(idx.y[,-1]), ncol = length(y) - 1)
		x <- (x[-1] + x[-length(x)]) / 2
		y <- (y[-1] + y[-length(y)]) / 2
	}
	# eventually we have (D %*% D %*% ...)
	# ready to be multiplied by as.vector(z)
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
				lambda, d, p, logscale
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
