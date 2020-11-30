feemcube <- function(x, ...) UseMethod('feemcube')

feemcube.list <- function(x, all.wavelengths, ...) {
	# must be a list of feem objects
	stopifnot(
		length(list(...)) == 0,
		sapply(x, inherits, 'feem')
	)
	scales <- sapply(x, attr, 'scale')
	# sort because order of union(...) or intersect(...) is not guaranteed
	dimensions <- Map(sort, Reduce(
		# This works because Map(f(a,b), A, B) returns
		# list(f(A[1], B[1]), f(A[2], B[2]), ...)
		function(a, b) Map(
			if (all.wavelengths) union else intersect,
			a, b
		),
		Map(
			function(l) attributes(l)[c('emission', 'excitation')],
			x
		)
	))
	feemcube(
		sapply(x, function(eem) eem[ # extract matching wavelengths or NAs
			match(dimensions$emission, attr(eem, 'emission')),
			match(dimensions$excitation, attr(eem, 'excitation'))
		], simplify = 'array'), # discards feem classes; doesn't matter
		dimensions$emission, dimensions$excitation,
		scales, names(x)
	)
}

feemcube.array <- function(x, emission, excitation, scales, names = NULL, ...) {
	if (missing(scales)) scales <- rep(1, dim(x)[3])
	stopifnot(
		length(list(...)) == 0,
		length(dim(x)) == 3,
		dim(x)[1:2] == c(length(emission), length(excitation)),
		is.null(names) || dim(x)[3] == length(names),
		dim(x)[3] == length(scales)
	)
	structure(
		x,
		emission = emission,
		excitation = excitation,
		scales = setNames(scales, names),
		dimnames = list(
			emission = emission,
			excitation = excitation,
			sample = names
		),
		class = 'feemcube'
	)
}

`[.feemcube` <- function(x, i, j, k, drop = TRUE) {
	ret <- NextMethod()
	if (missing(i)) i <- TRUE
	if (missing(j)) j <- TRUE
	if (missing(k)) k <- TRUE
	# special case: returning a cube
	if (length(dim(ret)) == 3) return(feemcube(
		ret,
		emission = attr(x, 'emission')[i],
		excitation = attr(x, 'excitation')[j],
		scales = attr(x, 'scales')[k],
		names = dimnames(ret)[[3]]
	))
	# special case: returning a FEEM
	# Only possible when drop = TRUE and choosing a single sample
	# but maybe multiple wavelengths (or all of them)
	if (length(dim(ret)) == 2 && length(seq_len(dim(x)[3])[k]) == 1)
		return(feem(
			ret,
			attr(x, 'emission')[i],
			attr(x, 'excitation')[j],
			attr(x, 'scales')[k]
		))
	ret
}

`[<-.feemcube` <- function(x, i, j, k, value) {
	# special case: assigning a cube or a FEEM
	if (inherits(value, 'feemcube') || inherits(value, 'feem')) {
		stopifnot( # wavelengths must match
			attr(x, 'emission')[i] == attr(value, 'emission'),
			attr(x, 'excitation')[j] == attr(value, 'excitation')
		)
		# scales should match, but we will proceed anyway
		rhs.scales <- attr(value,
			if (inherits(value, 'feem')) 'scale' else 'scales'
		)
		if (any(dif <- attr(x, 'scales')[k] != rhs.scales)) {
			# gather the sample names, if any
			if (is.null(rn <- dimnames(x)[[3]])) rn <- seq_len(dim(x)[3])
			# combine the LHS and RHS (could have different lengths!);
			# leave those that differ
			warn <- rbind(attr(x, 'scales')[k], rhs.scales)[, dif, drop = FALSE]
			# format the warning table
			warn <- rbind(
				format(c('',  rn[k][dif]), justify = 'right'), ' ',
				format(c('LHS', warn[1,]), justify = 'right'), ' ',
				format(c('RHS', warn[2,]), justify = 'right'), '\n'
			)
			warning(
				'Assigning from FEEM[s] with different scales:\n', warn
			)
		}
	}
	NextMethod()
}

as.list.feemcube <- function(x, ...) {
	stopifnot(length(list(...)) == 0)
	lapply(
		setNames(1:dim(x)[3], dimnames(x)[[3]]),
		function(i) x[,,i]
	)
}

as.data.frame.feemcube <- function(x, ...) {
	samples <- dimnames(x)[[3]]
	if (is.null(samples)) samples <- as.character(1:dim(x)[3])
	if (anyDuplicated(samples)) samples <- make.unique(samples)
	data.frame(
		emission = attr(x, 'emission')[slice.index(x, 1)][!is.na(x)],
		excitation = attr(x, 'excitation')[slice.index(x, 2)][!is.na(x)],
		intensity = x[!is.na(x)],
		sample = samples[slice.index(x, 3)][!is.na(x)]
	)
}

plot.feemcube <- function(
	x, xlab = quote(lambda[em]*', nm'), ylab = quote(lambda[ex]*', nm'),
	cuts = 128, col.regions = marine.colours(256), as.table = TRUE, ...
)
	levelplot(
		x = intensity ~ emission + excitation | sample,
		data = as.data.frame(x), xlab = xlab, ylab = ylab, cuts = cuts,
		col.regions = col.regions, as.table = as.table, ...
	)
