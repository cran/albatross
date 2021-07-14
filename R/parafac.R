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
		"Algorithm terminated abnormally due to ",
		"a problem with the constraints"
	)
	# assign dimnames for convenience
	rownames(ret$A) <- dimnames(cube)[[1]]
	rownames(ret$B) <- dimnames(cube)[[2]]
	rownames(ret$C) <- dimnames(cube)[[3]][subset]
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
feemcube.feemparafac <- function(x, ...) {
	stopifnot(length(list(...)) == 0)
	# the cube could have been stored directly or in an environment
	cube <- attr(x, 'cube')
	if (!is.null(envir <- attr(x, 'envir')))
		cube <- get(cube, envir = envir)
	# we may have been asked to process a subset of the samples
	if (is.null(subs <- attr(x, 'subset'))) subs <- TRUE
	cube[,,subs]
}

fitted.feemparafac <- function(object, ...) {
	stopifnot(length(list(...)) == 0)
	# access and subset the cube
	cube <- feemcube(object)
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
	feemcube(object) - fitted(object)
}

coef.feemparafac <- function(
	object, type = c(
		'all', 'scores', 'loadings', 'emission', 'excitation', 'samples'
	), ...
) {
	stopifnot(length(list(...)) == 0)
	cube <- feemcube(object)
	comps <- list(
		emission = list(comp = 'A', name = 'wavelength', val = attr(cube, 'emission')),
		excitation = list(comp = 'B', name = 'wavelength', val = attr(cube, 'excitation')),
		samples = list(comp = 'C', name = 'sample', val = .feemcsamples(cube))
	)
	switch(type <- match.arg(type),
		emission =, excitation =, samples = {
			set <- comps[[type]]
			values <- object[[set$comp]]
			as.data.frame(
				list(
					set$val[row(values)],
					as.vector(values),
					as.factor(col(values))
				),
				col.names = c(set$name, 'value', 'factor')
			)
		},
		scores = coef(object, 'samples'),
		all = lapply(
			setNames(nm = c('emission', 'excitation', 'samples')),
			function(n) coef(object, n)
		),
		loadings = do.call(rbind, lapply(c('emission', 'excitation'),
			function(n) cbind(
				coef(object, n),
				# kludge: uppercase "Emission" / "Excitation"
				mode = `substr<-`(n, 1, 1, 'E')
			)
		))
	)
}

compplot.surf <- function(X, ...) {
	cube <- feemcube(X)
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
	xyplot(
		x = value ~ wavelength | factor, groups = mode,
		data = coef(X, 'loadings'),
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
