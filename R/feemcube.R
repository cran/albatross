feemcube <- function(x, ...) UseMethod('feemcube')

feemcube.list <- function(x, all.wavelengths, ...) {
	# must be a list of feem objects
	stopifnot(
		length(list(...)) == 0,
		vapply(x, inherits, logical(1), 'feem')
	)
	scales <- vapply(x, attr, numeric(1), 'scale')
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
		vapply(
			x, function(eem) eem[ # extract matching wavelengths or NAs
				match(dimensions$emission, attr(eem, 'emission')),
				match(dimensions$excitation, attr(eem, 'excitation'))
			],
			matrix(
				numeric(),
				length(dimensions$emission),
				length(dimensions$excitation)
			)
		), # discards feem classes; doesn't matter
		dimensions$emission, dimensions$excitation,
		scales, names(x)
	)
}

feemcube.array <- function(x, emission, excitation, scales, names = NULL, ...) {
	if (missing(scales)) scales <- rep(1, dim(x)[3])
	stopifnot(
		length(list(...)) == 0,
		length(dim(x)) == 3, is.numeric(x),
		is.vector(emission, 'numeric'), is.vector(excitation, 'numeric'),
		is.vector(scales, 'numeric'),
		dim(x) == c(length(emission), length(excitation), length(scales)),
		is.null(names) || (
			dim(x)[3] == length(names) && is.vector(names) && is.atomic(names)
		)
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
	# attributes must index exactly as the array itself
	# this includes being able to index by dimnames, unfortunately
	em <- unname(setNames(attr(x, 'emission'), dimnames(x)[[1]])[i])
	ex <- unname(setNames(attr(x, 'excitation'), dimnames(x)[[2]])[j])
	sc <- unname(setNames(attr(x, 'scales'), dimnames(x)[[3]])[k])
	# special case: returning a cube
	if (length(dim(ret)) == 3) return(feemcube(
		ret, emission = em, excitation = ex, scales = sc,
		names = dimnames(ret)[[3]]
	))
	# special case: returning a FEEM
	# Only possible when drop = TRUE and choosing a single sample
	# but maybe multiple wavelengths (or all of them)
	if (length(dim(ret)) == 2 && length(sc) == 1)
		return(feem(ret, em, ex, sc))
	ret
}

`[<-.feemcube` <- function(x, i, j, k, value) {
	# special case: assigning a cube or a FEEM when the subset is a cube
	if (
		nargs() > 3 &&
		(inherits(value, 'feemcube') || inherits(value, 'feem'))
	) {
		if (missing(k)) k <- TRUE # we'll need this later
		# sanity check: wavelengths must match
		xsub <- x[
			if (missing(i)) TRUE else i,
			if (missing(j)) TRUE else j,
			k,
			drop = FALSE
		]
		stopifnot(
			attr(xsub, 'emission') == attr(value, 'emission'),
			attr(xsub, 'excitation') == attr(value, 'excitation')
		)
		# scales should match, but we will proceed anyway
		scales <- attr(xsub, 'scales')
		rhs.scales <- attr(value,
			if (inherits(value, 'feem')) 'scale' else 'scales'
		)
		if (any(dif <- scales != rhs.scales)) {
			# gather the sample names, if any
			rn <- .cubenames(x)
			# combine the LHS and RHS (could have different lengths!);
			# leave those that differ
			warn <- rbind(scales, rhs.scales)[, dif, drop = FALSE]
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
	mask <- !is.na(x)
	data.frame(
		emission = attr(x, 'emission')[slice.index(x, 1)][mask],
		excitation = attr(x, 'excitation')[slice.index(x, 2)][mask],
		intensity = x[mask],
		sample = .cubenames(x)[slice.index(x, 3)][mask],
		...
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
