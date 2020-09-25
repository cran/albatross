feemparafac <- function(X, ..., rescale = 3, retries = 10) {
	stopifnot(inherits(X, 'feemcube'))
	for (i in seq_len(retries)) {
		ret <- parafac(X, output = 'best', ...)
		if (ret$cflag != 2) break
	}
	if (ret$cflag == 2) stop(
		"Algorithm terminated abnormally due to",
		"a problem with the constraints"
	)
	# undo per-sample scaling
	ret$C <- ret$C * attr(X, 'scales')
	# move variance to a given mode
	if (!is.na(rescale)) {
		stopifnot(length(rescale) == 1, rescale %in% 1:3)
		factors <- LETTERS[1:3]
		for (f in factors[-rescale])
			ret <- rescale(ret, f, absorb = factors[rescale])
	}
	structure(ret,
		class = c('feemparafac', attr(ret, 'class')),
		cube = X
	)
}

fitted.feemparafac <- function(object, ...) {
	stopifnot(length(list(...)) == 0)
	# undo the scaling to make it comparable to the original cube
	feemcube(
		NextMethod() /
			attr(attr(object, 'cube'), 'scales')[
				slice.index(attr(object, 'cube'), 3)
			],
		attr(attr(object, 'cube'), 'emission'),
		attr(attr(object, 'cube'), 'excitation'),
		attr(attr(object, 'cube'), 'scales'),
		dimnames(attr(object, 'cube'))[[3]]
	)
}

residuals.feemparafac <- function(object, ...) {
	stopifnot(length(list(...)) == 0)
	# since scaling originally present in the cube is undone for object$C,
	# the cube has to undergo a similar transformation to match
	attr(object, 'cube') - fitted(object)
}

compplot.surf <- function(X, ...) plot(feemcube(with(X,
	lapply(setNames(nm = 1:ncol(A)), function (i)
		feem(
			A[,i] %o% B[,i],
			attr(attr(X, 'cube'), 'emission'),
			attr(attr(X, 'cube'), 'excitation')
		)
	)
), TRUE), ...)

compplot.xy <- function(
	X, xlab = quote(lambda*", nm"), ylab = "Factor value", as.table = T,
	auto.key = TRUE, type = 'l', ...
) xyplot(
	x = factor ~ wavelength | ncomp, groups = mode,
	data = do.call(rbind, lapply(1:ncol(X$A), function (i) {
		rbind(
			data.frame(
				wavelength = attr(attr(X, 'cube'), 'emission'),
				factor = X$A[,i],
				mode = 'Emission',
				ncomp = as.character(i)
			),
			data.frame(
				wavelength = attr(attr(X, 'cube'), 'excitation'),
				factor = X$B[,i],
				mode = 'Excitation',
				ncomp = as.character(i)
			)
		)
	})),
	xlab = xlab, ylab = ylab, as.table = as.table, auto.key = auto.key,
	type = type, ...
)

plot.feemparafac <- function(x, type = c('image', 'lines'), ...)
	switch(
		match.arg(type),
		image = compplot.surf(x, ...),
		lines = compplot.xy(x, ...)
	)
