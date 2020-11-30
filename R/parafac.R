feemparafac <- function(
	X, ..., rescale = 3, retries = 10, subset = TRUE, envir = NULL
) {
	cube <- if (is.null(envir)) X else get(X, envir = envir)
	stopifnot(inherits(cube, 'feemcube'))
	for (i in seq_len(retries)) {
		ret <- parafac(cube[,,subset], output = 'best', ...)
		if (ret$cflag != 2) break
	}
	if (ret$cflag == 2) stop(
		"Algorithm terminated abnormally due to",
		"a problem with the constraints"
	)
	# undo per-sample scaling
	ret$C <- ret$C * attr(cube, 'scales')[subset]
	# move variance to a given mode
	if (!is.na(rescale)) {
		stopifnot(length(rescale) == 1, rescale %in% 1:3)
		factors <- LETTERS[1:3]
		for (f in factors[-rescale])
			ret <- rescale(ret, f, absorb = factors[rescale])
	}
	structure(ret,
		class = c('feemparafac', oldClass(ret)),
		cube = X, subset = subset, envir = envir
	)
}

# extract the cube from a feemcube object, directly or by reference
.pfcube <- function(X) {
	# the cube could have been stored directly or in an environment
	cube <- attr(X, 'cube')
	if (!is.null(envir <- attr(X, 'envir')))
		cube <- get(cube, envir = envir)
	# we may have been asked to process a subset of the samples
	if (is.null(subs <- attr(X, 'subset'))) subs <- TRUE
	cube[,,subs]
}

fitted.feemparafac <- function(object, ...) {
	stopifnot(length(list(...)) == 0)
	# access and subset the cube
	cube <- .pfcube(object)
	# since scaling originally present in the cube is undone for object$C,
	# the fitted cube has to undergo an opposite transformation to match
	feemcube(
		NextMethod() / attr(cube, 'scales')[slice.index(cube, 3)],
		attr(cube, 'emission'),
		attr(cube, 'excitation'),
		attr(cube, 'scales'),
		dimnames(cube)[[3]]
	)
}

residuals.feemparafac <- function(object, ...) {
	stopifnot(length(list(...)) == 0)
	.pfcube(object) - fitted(object)
}

compplot.surf <- function(X, ...) {
	cube <- .pfcube(X)
	plot(feemcube(with(X,
		lapply(setNames(nm = 1:ncol(A)), function (i)
			feem(
				A[,i] %o% B[,i],
				attr(cube, 'emission'),
				attr(cube, 'excitation')
			)
		)
	), TRUE), ...)
}

compplot.xy <- function(
	X, xlab = quote(lambda*", nm"), ylab = "Factor value", as.table = T,
	auto.key = TRUE, type = 'l', ...
) {
	cube <- .pfcube(X)
	xyplot(
		x = factor ~ wavelength | ncomp, groups = mode,
		data = do.call(rbind, lapply(1:ncol(X$A), function (i) {
			rbind(
				data.frame(
					wavelength = attr(cube, 'emission'),
					factor = X$A[,i],
					mode = 'Emission',
					ncomp = as.character(i)
				),
				data.frame(
					wavelength = attr(cube, 'excitation'),
					factor = X$B[,i],
					mode = 'Excitation',
					ncomp = as.character(i)
				)
			)
		})),
		xlab = xlab, ylab = ylab, as.table = as.table, auto.key = auto.key,
		type = type, ...
	)
}

plot.feemparafac <- function(x, type = c('image', 'lines'), ...)
	switch(
		match.arg(type),
		image = compplot.surf(x, ...),
		lines = compplot.xy(x, ...)
	)
