.splitcombine <- function(samples, splits) {
	# splits must be a scalar or a 2-element vector of whole numbers
	# splits[1] should be even because we want to combine into halves
	# we should also have at least as many samples as splits[1]
	stopifnot(
		length(splits) <= 2, splits == round(splits),
		splits[1] %% 2 == 0,
		length(samples) >= splits[1],
		is.na(splits[2]) || splits[2] <= choose(splits[1], splits[1]/2)/2
	)
	if (is.na(splits[2])) splits[2] <- choose(splits[1], splits[1]/2)/2

	# we now have `splits` groups to recombine into halves
	groups <- split(samples, rep_len(1:splits[1], length(samples)))

	# use combn() to generate combinations of groups[[i]]
	# taken splits[1]/2 at at time i.e. build halves from them
	unique(combn(
		groups, splits[1] / 2, function(cmb) {
			# not sure where the names are coming from, but they
			# break the tests for identity in unique()
			left <- unname(sort(unlist(cmb)))
			# generate the other half right away (AB => CD)
			right <- unname(sort(setdiff(samples, left)))
			# this approach would give us both pairs AB vs CD and CD vs AB
			# so we sort both halves, place the half with the
			# lowest-numbered sample first, then eliminate duplicates
			list(left, right)[order(c(min(left), min(right)))]
		}, FALSE
	))[1:splits[2]]
}

.randomhalves <- function(samples, N) {
	# should have at least 2 samples in order to separate into halves
	stopifnot(length(samples) >= 2)

	half <- floor(length(samples) / 2)
	replicate(N, {
		samples <- sample(samples)
		list(samples[1:half], samples[-(1:half)])
	}, FALSE)
}

feemsplithalf <- function(
	cube, nfac, splits, random, groups, fixed, ..., progress = TRUE
) {
	tests <- if (
		!missing(fixed) &&
		missing(splits) && missing(random) && missing(groups)
	) {
		# must be a list of pairs
		# pair is a list of two vectors containing integer indices
		stopifnot(
			is.list(fixed), lengths(fixed) == 2,
			unlist(fixed) == round(unlist(fixed))
		)
		for (i in seq_along(fixed)) {
			if (length(s <- intersect(
				fixed[[i]][[1]], fixed[[i]][[2]]
			)) > 0)
				stop(
					'Both halves in fixed[[', i, ']] contain the ',
					'following samples: ', paste(s, collapse = ', ')
				)
		}
		fixed
	} else if (
		missing(splits) != missing(random) && # exactly one of them set
		missing(fixed)
	) {
		# if not performing stratified sampling, create one fake group
		# encasing all samples
		if (missing(groups)) groups <- rep_len(1, dim(cube)[3])
		stopifnot(
			(if (is.list(groups)) lengths else length)(groups) == dim(cube)[3]
		)
		# list(list(half, half), ...)
		Reduce(
			# given same-sized per-group lists of halves,
			function(a, b) lapply(
				# concatenate individual halves between groups
				seq_along(a), function(i) Map(c, a[[i]], b[[i]])
			),
			lapply( # list(group = list(list(half, half), ...), ...)
				split(1:dim(cube)[3], groups, drop = TRUE),
				if (!missing(splits))
					function(s) .splitcombine(s, splits)
				else if (!missing(random))
					function(s) .randomhalves(s, random)
			)
		)
	} else stop(
		'Please either request split-combine or random halves ',
		'(optionally stratified by groups), ',
		'or provide fixed halves.'
	)
	# organise two parallel arrays:
	# ((half1, half2), (half3, half4), ...) =>
	# ((half1, half2) for ncomp1, (half1, half2) for ncomp2, ...) =>
	#  (half1 for ncomp1, half2 for ncomp1, half1 for ncomp2, ...)
	slices <- unlist(rep(tests, each = length(nfac)), FALSE)
	# (ncomp1, ...) => (ncomp1, ncomp1, ncomp2, ncomp2, ...)
	facarg <- lapply(
		rep(nfac, each = 2, times = length(tests)),
		function(n) list(nfac = n)
	)

	factors <- bootparafac(
		cube, slices, ..., args = facarg, progress = progress
	)

	# okay, now we have a list of parafac results:
	# half_a_1, nfac1
	# half_a_2, nfac1
	# half_a_1, nfac2
	# half_a_2, nfac2
	# ...
	# half_b_1, nfac1
	# half_b_2, nfac1
	# half_b_1, nfac2
	# half_b_2, nfac2
	# ...

	# (with "a", "b" referring to different tests; 1 & 2 referring to
	# halves and nfac1,2,... referring to numbers of components)

	# group the list by halves / component numbers / tests
	dim(factors) <- c(2, length(nfac), length(tests))
	dimnames(factors) <- list(half = NULL, nfac = nfac, test = NULL)

	tcc <- lapply(
		setNames(seq_along(nfac), nfac), # iterating over numbers of factors,
		function(ifac) structure(
			vapply(
				1:dim(factors)[3], function(grp) # for each grouping,
					vapply( # for each mode,
						c('A','B'), function(mode) # TCC for this half and mode
							diag(as.matrix(congru(
								factors[[1,ifac,grp]][[mode]],
								factors[[2,ifac,grp]][[mode]]
							))),
							# only for matching components
						numeric(nfac[ifac])
					),
				matrix(numeric(), nfac[ifac], 2)
			),
			dimnames = list(
				factor = NULL, mode = c('Emission', 'Excitation'),
				test = NULL
			)
		)
	)

	structure(
		list(factors = factors, tcc = tcc, nfac = nfac),
		class = 'feemsplithalf'
	)
}

print.feemsplithalf <- function(x, ...) {
	stopifnot(length(list(...)) == 0)
	cat("Split-half: minimal TCC between matching components\n")
	print(vapply(x$tcc, min, numeric(1)))
	invisible(x)
}

shtccplot <- function(
	x, xlab = 'Number of components',
	ylab = 'Minimum TCC between halves', jitter.x = TRUE, ...
) {
	factor <- NULL # R CMD check vs. xyplot(groups = ...)
	xyplot(
		tcc ~ nfac, coef(x, 'tcc'), group = factor,
		xlab = xlab, ylab = ylab, jitter.x = jitter.x, ...
	)
}

shxyplot <- function(
	x, xlab = quote(lambda*", nm"), ylab = 'Factor value', as.table = T, ...
) {
	df <- coef(x, 'factors')
	factor <- test <- half <- NULL # R CMD check vs xyplot(groups = ...)
	xyplot(
		value ~ wavelength | mode + nfac, df, groups = paste(factor, test, half),
		par.settings = list(superpose.line = list(col = rep(
			trellis.par.get('superpose.line')$col, each = 2 * dim(x$factors)[3]
		))), scales = list(x = list(relation = 'free')),
		type = 'l', xlab = xlab, ylab = ylab, as.table = as.table, ...
	)
}

plot.feemsplithalf <- function(x, kind = c('tcc', 'factors'), ...) {
	switch(match.arg(kind),
		tcc = shtccplot(x, ...),
		factors = shxyplot(x, ...)
	)
}

# all loadings from all split-half objects
coefshfact <- function(x)
	do.call(rbind, Map(
		function(x, test, half) {
			ret <- cbind(
				coef(x, 'loadings'),
				nfac = as.factor(ncol(x$A)),
				test = test,
				half = half
			)
			ret$subset <- rep.int(list(attr(x, 'subset')), nrow(ret))
			ret
		}, x$factors, slice.index(x$factors, 3), slice.index(x$factors, 1)
	))

coefshtcc <- function(object)
	# concatenate by number of factors
	do.call(rbind, lapply(seq_along(object$nfac), function(i) {
		# the feemparafac objects store indices subsetting the cube
		subsets <- apply(
			object$factors[,i,, drop = FALSE], 3, lapply, attr, 'subset'
		)
		# min over mode (emission / excitation) because we use
		# the same quantity to match components
		tcc <- apply(object$tcc[[i]], c(1,3), min)
		data.frame(
			factor = as.factor(row(tcc)),
			tcc = as.vector(tcc),
			test = as.factor(col(tcc)),
			subset = I(subsets[col(tcc)]),
			nfac = as.factor(object$nfac[i])
		)
	}))

coef.feemsplithalf <- function(object, kind = c('tcc', 'factors'), ...) {
	stopifnot(length(list(...)) == 0)
	switch(match.arg(kind),
		factors = coefshfact(object),
		tcc = coefshtcc(object)
	)
}

feemcube.feemsplithalf <- function(x, ...) {
	stopifnot(length(list(...)) == 0)
	x <- x$factors[[1]]
	get(attr(x, 'cube'), envir = attr(x, 'envir'))
}
