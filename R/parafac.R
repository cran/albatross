# Everything related to performing the PARAFAC decomposition itself is
# here. We wrap multiway::parafac and provide restarts on unsuccessful
# fits, rescaling (e.g. make loadings of unit norm), plots and coef()
# methods returning long-format data.frames.

feemparafac <- function(
	X, ..., const = rep('nonneg', 3), ctol = 1e-6,
	rescale = 3, retries = 10, subset = TRUE, envir = NULL
) {
	# see feemsplithalf and friends for why this is needed
	cube <- if (is.null(envir)) X else get(X, envir = envir)
	# check parameters for errors
	stopifnot(inherits(cube, 'feemcube'))
	if (!is.na(rescale)) stopifnot(length(rescale) == 1, rescale %in% 1:3)
	# run the actual PARAFAC
	time <- system.time(
		for (i in seq_len(retries)) {
			ret <- parafac(
				cube[,,subset], output = 'best', ...,
				const = const, ctol = ctol
			)
			if (ret$cflag != 2) break
		}
	)
	if (ret$cflag == 2) stop(
		'Algorithm terminated abnormally after ', retries,
		' tries due to a problem with the constraints'
	)
	# assign dimnames for convenience
	rownames(ret$A) <- dimnames(cube)[[1]]
	rownames(ret$B) <- dimnames(cube)[[2]]
	rownames(ret$C) <- dimnames(cube)[[3]][subset]
	# undo per-sample scaling
	ret$C <- ret$C * attr(cube, 'scales')[subset]
	# move variance to the given mode
	if (!is.na(rescale)) {
		factors <- LETTERS[1:3]
		for (f in factors[-rescale])
			ret <- rescale(ret, f, absorb = factors[rescale])
	}
	structure(ret,
		class = c('feemparafac', oldClass(ret)),
		cube = X, subset = subset, envir = envir, time = time
	)
}

# Turns out, doing cube[,,subset] over and over is very expensive,
# so do not subset the cube if we can help it
.metadata.feemparafac <- function(x) {
	# the cube could have been stored directly or in an environment
	cube <- attr(x, 'cube')
	if (!is.null(envir <- attr(x, 'envir')))
		cube <- get(cube, envir = envir)
	# we may have been asked to process a subset of the samples
	if (is.null(subs <- attr(x, 'subset'))) subs <- TRUE
	list(
		emission = attr(cube, 'emission'),
		excitation = attr(cube, 'excitation'),
		names = .cubenames(cube)[subs],
		cube = cube,
		subset = subs
	)
}

# extract the cube from a feemcube object, directly or by reference
feemcube.feemparafac <- function(x, ...) {
	stopifnot(length(list(...)) == 0)
	ret <- .metadata.feemparafac(x)
	ret$cube[,,ret$subset]
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
	ret <- feemcube(object) - fitted(object)
	structure(ret, class = c('feem.resid', class(ret)))
}

coef.feemparafac <- function(
	object, type = c(
		'all', 'scores', 'loadings', 'surfaces', 'emission', 'excitation', 'samples'
	), ...
) {
	stopifnot(length(list(...)) == 0)
	meta <- .metadata.feemparafac(object)
	comps <- list(
		emission = list(comp = 'A', name = 'wavelength', val = meta$emission),
		excitation = list(comp = 'B', name = 'wavelength', val = meta$excitation),
		samples = list(comp = 'C', name = 'sample', val = meta$names)
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
		)),
		surfaces = {
			do.call(rbind, lapply(1:ncol(object$A), function(i)
				cbind(
					as.data.frame(feem(
						object$A[,i] %o% object$B[,i],
						meta$emission,
						meta$excitation
					)),
					factor = ordered(i, 1:ncol(object$A))
				)
			))
		}
	)
}

compplot.surf <- function(
	X,
	xlab = pgtq("lambda[em]*', nm'", translate),
	ylab = pgtq("lambda[ex]*', nm'", translate),
	col.regions = marine.colours(256), cuts = 255, as.table = TRUE, ...,
	translate = FALSE
) levelplot(
	intensity ~ emission + excitation | factor, coef(X, 'surfaces'),
	xlab = xlab, ylab = ylab, col.regions = col.regions, cuts = cuts,
	as.table = as.table, ...
)

compplot.xy <- function(
	X, xlab = pgtq("lambda*', nm'", translate),
	ylab = pgt("Factor value", translate),
	as.table = T, auto.key = TRUE, type = 'l', ..., translate = FALSE
) {
	xyplot(
		x = value ~ wavelength | factor, groups = pgt(mode, translate),
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

reorder.feemparafac <- function(x, neworder, like, ...) {
	nfac <- ncol(x$A)
	if (missing(like)) {
		if (nfac == 1) {
			# workaround multiway <= 1.0-6 problem;
			# no point reordering 1-component models anyway
			stopifnot(length(neworder) == 1, neworder == 1)
			return(x)
		}
		return(NextMethod())
	}
	stopifnot(
		# disallow extra arguments for our special case
		missing(neworder), length(list(...)) == 0,
		# for now, only reorder like models with the same nfac
		ncol(like$A) == ncol(x$A)
	)

	# no way to reorder 1-component models
	if (nfac == 1) return(x)

	# Tucker's congruence coefficient for emission & excitation components
	tcc <- as.matrix(pmin(congru(like$A, x$A), congru(like$B, x$B)))

	perm <- integer(nfac)
	for (i in seq_along(perm)) {
		# more than one component could fit perfectly
		# (despite it usually shouldn't), so choose the first match
		next.match <- which(tcc == max(tcc), arr.ind=T)[1,]
		# record the match and the distance
		perm[next.match[1]] <- next.match[2]
		# make sure that this pair won't match with anything else
		tcc[next.match[1],] <- -Inf
		tcc[,next.match[2]] <- -Inf
	}

	reorder(x, neworder = perm)
}

# argmin over c[] ||A - diag(c) %*% B||^2
.matscale <- function(A, B) {
	stopifnot(ncol(A) == ncol(B))
	vapply(
	# for every column:
	#  min over c[j] ||A[,j] - c[j] * B[,j]||
	#  c[j] = B[,j]^T A[,j] / B[,j]^T B[,j]
		1:ncol(A), function(j)
			crossprod(B[,j,drop=FALSE], A[,j,drop=FALSE]) /
				crossprod(B[,j,drop=FALSE]),
		numeric(1)
	)
}

rescale.feemparafac <- function(x, mode, newscale, absorb, like, ...) {
	if (missing(like)) return(NextMethod())
	stopifnot(
		# disallow setting newscale when given a reference model
		missing(newscale),
		# no extra arguments
		length(list(...)) == 0,
		# make sure that the number of factors matches
		ncol(x$A) == ncol(like$A)
	)
	# our defaults are slightly different from rescale.parafac
	if (missing(mode)) mode <- c('A', 'B')
	if (missing(absorb)) absorb <- 'C'
	stopifnot(
		# must not absorb the rescaling into modes that are rescaled
		length(intersect(mode, absorb)) == 0,
		# only three-way arrays supported
		union(mode, absorb) %in% c('A', 'B', 'C'),
		!anyDuplicated(mode), length(absorb) == 1
	)

	scaling <- matrix(1, ncol(x$A), 3, dimnames = list(NULL, LETTERS[1:3]))
	unscaling <- rep(1, nrow(scaling))
	for (cn in mode) {
		scaling[,cn] <- .matscale(like[[cn]], x[[cn]])
		unscaling <- unscaling / scaling[,cn]
	}
	scaling[,absorb] <- unscaling

	x$A <- x$A %*% diag(scaling[,'A'], nrow(scaling))
	x$B <- x$B %*% diag(scaling[,'B'], nrow(scaling))
	x$C <- x$C %*% diag(scaling[,'C'], nrow(scaling))

	x
}

print.feemparafac <- function(x, ...) {
	nfac <- ncol(x$A)
	cat(
		'\nPARAFAC with ', nfac, ' factors\n\n',
		'Constraints:\n',
		sep = ''
	)
	const <- setNames(x$const, LETTERS[seq_along(x$const)])
	const[x$fixed] <- 'fixed'
	const[x$struc] <- paste(const[x$struc], 'struct', sep = '+')
	print(noquote(const))
	cat(
		'\nFit Information:\n',
		'  SSE = ', x$SSE, '\n',
		'  R^2 = ', x$Rsq, '\n',
		'  GCV = ', x$GCV, '\n',
		'  EDF = ', sum(x$edf), '\n\n',
		'Converged: ', x$cflag == 0, ' (', x$iter, ' iterations)\n',
		sep = ""
	)
	invisible(x)
}
