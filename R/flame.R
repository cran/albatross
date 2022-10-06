# FLuorescence and scAttering Model Estimation, see ?feemflame for a
# high-level overview and ?cmf for implementation details. We use our
# own constrained matrix factorisation because it's easy enough given
# our dependency on CMLS.

wcmls <- function(X, A, W, ..., struc = NULL)
	matrix( # NB: vapply with scalar FUN.VALUE returns a vector
		vapply(seq_len(ncol(X)), function(j) {
			if (any(W[, j])) as.vector(cmls(
				A * W[, j, drop = TRUE],
				X[, j, drop = FALSE] * W[, j, drop = FALSE],
				..., struc = if (is.null(struc)) NULL
					else t(struc[j, , drop = FALSE])
			))
			else numeric(ncol(A)) # don't divide 0 by 0
		}, numeric(ncol(A))),
		nrow = ncol(X),
		byrow = TRUE # t(B) is built column by column
	)

cmf <- function(
	X, nfac = 1,
	const = list(list(const = 'nonneg'), list(const = 'nonneg')),
	start = c('svd', 'random'), ctol = 1e-4, maxiter = 10
) {
	stopifnot(length(dim(X)) == 2)

	# NAs mean zero error weight
	W <- !is.na(X)
	X[!W] <- 0 # zero weight still needs a defined value

	# we need both C and S to compute the first SSE value
	if (inherits(start, 'cmf')) {
		A <- start[[1]]
		B <- start[[2]]
	} else switch(match.arg(start),
		svd = {
			# crudely force SVD to nonnegative
			Xsvd <- svd(X, nfac, nfac)
			# rescale C and S to be of comparable norms
			Sigma <- diag(sqrt(Xsvd$d), nfac, nfac)
			A <- abs(Xsvd$u %*% Sigma)
			B <- abs(Xsvd$v %*% Sigma)
		},
		random = {
			A <- matrix(runif(nrow(X) * nfac), nrow(X), nfac)
			B <- matrix(runif(ncol(X) * nfac), ncol(X), nfac)
			# rescale C and S to be of comparable norms
			Xn <- sqrt(norm(X))
			A <- A / norm(A) * Xn
			B <- B / norm(B) * Xn
		}
	)

	SSE <- rep(sum(W * (X - tcrossprod(A, B))^2), 2)
	i <- 0

	repeat {
		# ALS update
		B <- do.call(wcmls, c(list(X, A, W), const[[1]]))
		A <- do.call(wcmls, c(list(t(X), B, t(W)), const[[2]]))

		# stop criteria update
		i <- i + 1
		SSE <- c(SSE[2], sum(W * (X - tcrossprod(A, B))^2))

		# stop criteria test
		if (i >= maxiter) break
		if (abs(diff(SSE))/SSE[2] <= ctol) break
	}

	structure(list(A, B), class = 'cmf')
}

fitted.cmf <- function(object, ...) tcrossprod(object[[1]], object[[2]])

feemflame <- function(
	X, ffac, sfac,
	maxiter = 32, widths = rep(25, 4), Raman.shift = 3400,
	ctol = 1e-4, progress = TRUE
) {
	stopifnot(inherits(X, 'feemcube'))

	fl <- feemparafac(
		feemscatter(X, widths, 'omit', progress = FALSE),
		nfac = ffac, verbose = FALSE
	)
	sc.struc <- matrix(
		.scatter.mask(X[,,1], widths, Raman.shift),
		ncol = sfac, nrow = prod(dim(X)[-3]) # use recycling
	)
	sc <- 'svd' # used as start= argument of cmf() first time

	i <- 0
	SSX <- sum(X^2, na.rm = TRUE)
	SSE <- rep(sum((X - fitted(fl))^2, na.rm = TRUE), 2)
	maxvtol <- vtol <- NULL

	# TODO: early termination by catching interrupt?
	repeat {
		# ALS step
		sc <- cmf(
			matrix(X - fitted(fl), prod(dim(X)[-3]), dim(X)[3]),
			nfac = sfac, const = list(
				list(const = 'nonneg'),
				list(const = 'nonneg', struc = sc.struc)
			),
			start = sc
		)
		fl <- feemparafac(
			Xfl <- X - array(fitted(sc), dim(X)),
			nfac = ffac, maxit = 10,
			nstart = 1, Astart = fl$A, Bstart = fl$B, Cstart = fl$C,
			verbose = FALSE
		)

		# update stop criteria
		i <- i + 1
		SSE[1] <- SSE[2]
		# Note that Xfl has fitted(MCR) already subtracted
		SSE[2] <- sum((Xfl - fitted(fl))^2, na.rm = TRUE)
		vtol <- abs(diff(SSE))/SSE[2]

		if (progress) {
			if (!exists('pb')) {
				pb <- vtolProgressBar(ctol)
				on.exit(pb$close())
			}
			pb$up(SSE[2]/SSX, vtol)
		}

		# test stop criteria
		if (i >= maxiter) break
		if (vtol <= ctol) break
	}

	structure(list(fl = fl, sc = sc), class = 'feemflame')
}

fitted.feemflame <- function(object, ...) {
	stopifnot(length(list(...)) == 0)
	# adding a vector doesn't perform dimensionality check
	fitted(object$fl) + as.vector(fitted(object$sc))
}

residuals.feemflame <- function(object, ...) {
	stopifnot(length(list(...)) == 0)
	# residuals are X - fitted(MCR) - fitted(PARAFAC)
	# PARAFAC fits X - fitted(MCR) *and* is fitted last
	# coincidentally, residuals of PARAFAC are the residuals of the model
	resid(object$fl)
}

feemcube.feemflame <- function(x, ...) {
	stopifnot(length(list(...)) == 0)
	# PARAFAC fits X - fitted(MCR) *and* is fitted last
	# therefore, the original X is the PARAFAC cube + fitted(MCR)
	feemcube(x$fl) + as.vector(fitted(x$sc))
}

coef.feemflame <- function(
	object, type = c(
		'fluorescence',
		'scores', 'loadings', 'emission', 'excitation', 'samples',
		'scattering', 'sc.scores', 'sc.loadings'
	), ...
) {
	stopifnot(length(list(...)) == 0)
	cube <- feemcube(object)
	switch(type <- match.arg(type),
		fluorescence = coef(object$fl),
		scores=, loadings=, emission=, excitation=,
		samples = coef(object$fl, type),
		scattering = lapply(
			setNames(
				c('sc.scores', 'sc.loadings'),
				c('scores', 'loadings')
			), function(type) coef(object, type)
		),
		sc.scores = data.frame(
			sample = .cubenames(cube)[row(object$sc[[2]])],
			value = as.vector(object$sc[[2]]),
			factor = as.factor(col(object$sc[[2]]))
		),
		sc.loadings = {
			loads <- array(
				object$sc[[1]], c(dim(cube)[1:2], ncol(object$sc[[1]]))
			)
			data.frame(
				emission = attr(cube, 'emission')[slice.index(loads, 1)],
				excitation = attr(cube, 'excitation')[slice.index(loads, 2)],
				value = as.vector(loads),
				factor = as.factor(slice.index(loads, 3))
			)
		}
	)
}

.flame.image.plot <- function(
	x, xlab = quote(lambda[em] * ", nm"),
	ylab = quote(lambda[ex] * ", nm"),
	cuts = 128, col.regions = marine.colours(256), as.table = TRUE, ...
) {
	cube <- feemcube(x)
	fl <- as.data.frame(feemcube(
		array(
			multiway::krprod(x$fl$B, x$fl$A),
			c(dim(cube)[1:2], ncol(x$fl$A))
		),
		attr(cube, 'emission'), attr(cube, 'excitation')
	))
	colnames(fl)[colnames(fl) == 'intensity'] <- 'value'
	colnames(fl)[colnames(fl) == 'sample'] <- 'factor'
	fl <- cbind(fl, kind = 'Fluorescence')

	sl <- cbind(coef(x, 'sc.loadings'), kind = 'Scattering')
	levelplot(
		value ~ emission + excitation | paste(kind, factor),
		rbind(fl, sl),
		xlab = xlab, ylab = ylab, col.regions = col.regions,
		cuts = cuts, as.table = as.table, ...
	)
}


plot.feemflame <- function(
	x, type = c('both', 'fl.image', 'fl.lines'), ...
) {
	switch(type <- match.arg(type),
		fl.image = plot(x$fl, 'image', ...),
		fl.lines = plot(x$fl, 'lines', ...),
		both = .flame.image.plot(x)
	)
}
