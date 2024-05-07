# (Bro, 1998): matricised Tucker model is X ~ A * G * (C %kr% B)'
parafac_tucker_core_dense <- function(X, A, B, C) {
	X <- unclass(X)
	dim(X) <- c(dim(X)[1], prod(dim(X)[-1]))
	# With no missing data, this is a matter of computing a few
	# pseudoinverses.
	# Ghat = A^+ X ((C %kr% B)')^+ =
	# (Petersen&Petersen, 2012, The Matrix Cookbook, (510))
	# = A^+ X (C' %kr% B')^+ =
	# Since pseudoinverse is super-linear in cost, it's cheaper to
	# compute two of them for smaller matrices.
	# (Petersen&Petersen, 2012, The Matrix Cookbook, (513))
	# = A^+ X (C'^+ %kr% B'^+)
	array(
		pinv(A) %*% X %*% kronecker(pinv(t(C)), pinv(t(B))),
		c(ncol(A), ncol(B), ncol(C))
	)
}

# Let D = (C %kr% B)'.
# With missing data marked as W = 0, our task is to minimise:
# F(G) = 1/2 || W * (A %*% G %*% D - X) ||^2
parafac_tucker_core_iterative <- function(
	X, A, B, C,
	method = 'L-BFGS-B',
	control = list(maxit = 1e4, factr = 1),
	...
) {
	X <- unclass(X)
	dim(X) <- c(dim(X)[1], prod(dim(X)[-1]))
	W = !is.na(X)
	X[!W] <- 0
	D <- t(C %x% B)

	# Iterative methods don't involve computing C %x% B %x% A and are
	# thus much kinder on system resources.
	start <- array(0, c(ncol(A), ncol(B), ncol(C)))
	start[cbind(1:ncol(A), 1:ncol(B), 1:ncol(C))] <- 1
	ret <- optim(
		as.vector(start),
		function(G) .5 * sum(
			(W * (A %*% matrix(G, ncol(A), ncol(B) * ncol(C)) %*% D - X))^2
		),
		function(G) t(A) %*% (
			W^2 * (A %*% matrix(G, ncol(A), ncol(B) * ncol(C)) %*% D - X)
		) %*% t(D),
		method = method, control = control, ...
	)
	if (ret$convergence == 51) warning(ret$message) else
	if (ret$convergence != 0)
		warning("Convergence not achieved: ", ret$message)
	array(ret$par, c(ncol(A), ncol(B), ncol(C)))
}

parafac_tucker_core_exact <- function(X, A, B, C) {
	# Instead of elementwise multiplication by weight matrix, we can
	# unfold everything into vectors and drop the missing points.
	mask <- !is.na(X)
	# F(G) = 1/2 || vec(A %*% G %*% D - X)[mask] ||^2 =
	# 1/2 || vec(A %*% G %*% D)[mask] - vec(X)[mask] ||^2 =
	# (Petersen&Petersen, 2012, The Matrix Cookbook, (520))
	# 1/2 || ((D' %kr% A) %*% vec G)[mask] - vec(X)[mask] ||^2 =
	# 1/2 || (D' %kr% A)[mask,] %*% vec G - vec(X)[mask] ||^2.

	# dF/dG = 0
	# ((D' %kr% A)[mask,] %*% vec G - vec(X)[mask])' %*% (D' %kr% A)[mask,] = 0
	# assuming that ((D' %kr% A)[mask,])^+ is well-behaved
	# ((D' %kr% A)[mask,] %*% vec G - vec(X)[mask])' = 0
	# (D' %kr% A)[mask,] %*% vec G = vec(X)[mask,]
	# With the exception of the mask, ths matches the formula (5) given
	# in (Bro&Kiers, 2003).
	array(
		qr.solve(
			(C %x% B %x% A)[mask,,drop = FALSE],
			as.vector(X)[mask]
		),
		c(ncol(A), ncol(B), ncol(C))
	)
}

feemcorcondia <- function(
	model, divisor = c('nfac', 'core'),
	kind = c('pinv', 'iterative', 'vec'), ...
) {
	divisor <- match.arg(divisor)
	nfac <- ncol(model$A)
	cube <- feemcube(model)
	if (missing(kind))
		kind <- if (any(is.na(cube))) 'iterative' else 'pinv'
	core <- switch(match.arg(kind),
		pinv = parafac_tucker_core_dense,
		iterative = parafac_tucker_core_iterative,
		vec = parafac_tucker_core_exact
	)(
		cube, model$A, model$B, model$C / attr(cube, 'scales'), ...
	)
	core.ideal <- array(0, rep(nfac, 3))
	core.ideal[cbind(1:nfac, 1:nfac, 1:nfac)] <- 1
	structure(
		100 * (1 -
			sum((core - core.ideal)^2) /
			sum(switch(divisor, nfac = core.ideal, core = core)^2)
		),
		divisor = divisor,
		core = core,
		class = 'feemcorcondia'
	)
}

print.feemcorcondia <- function(x, ...) {
	print(`attributes<-`(x, NULL))
	cat(
		'(Divisor: ', sQuote(attr(x, 'divisor')),
		'; see attr(,"core") for the least squares Tucker3 core)\n', sep = ''
	)
	invisible(x)
}
