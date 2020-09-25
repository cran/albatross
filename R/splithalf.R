feemsplithalf <- function(cube, nfac, splits, random, ...) {
	# list(list(half, half), ...)
	tests <- if (!missing(splits) && missing(random)) {
		# should be even because we want to combine into halves
		stopifnot(splits %% 2 == 0)

		# a list of length(.) == splits
		groups <- split(1:dim(cube)[3], rep_len(1:splits, dim(cube)[3]))

		# use combn() to generate combinations of groups[[i]]
		# taken splits/2 at at time i.e. build halves from them
		unique(combn(
			groups, splits / 2, function(cmb) {
				# not sure where the names are coming from, but break
				# the tests for identity in unique()
				left <- unname(sort(unlist(cmb)))
				# generate the other half right away (AB => CD)
				right <- unname(sort((1:dim(cube)[3])[-left]))
				# this approach would give us both pairs AB vs CD and CD vs AB
				# so we sort both halves, place the half with the
				# lowest-numbered sample first, then eliminate duplicates
				list(left, right)[order(c(min(left), min(right)))]
			}, FALSE
		))
	} else if (!missing(random) && missing(splits)) {
		half <- floor(dim(cube)[3] / 2)
		replicate(random, {
			indices <- sample.int(dim(cube)[3])
			list(indices[1:half], indices[-(1:half)])
		}, FALSE)
	} else stop(
		"Please provide either splits = n.groups or random = n.shuffles"
	)

	# organise two parallel arrays:
	# ((half1, half2), ...) =>
	#  (half1 for ncomp1, half2 for ncomp1, half1 for ncomp2, ...)
	slices <- rep(unlist(tests, FALSE), each = length(nfac))
	# (ncomp1, ...) => (ncomp1, ncomp1, ncomp2, ncomp2, ...)
	facarg <- lapply(
		rep(nfac, each = 2, times = length(tests)),
		function(n) list(nfac = n)
	)

	factors <- bootparafac(cube, slices, ..., args = facarg)

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

	# (with "a", "b" referring to different tests; 1 & 2 referring to
	# halves and nfac1,2,... referring to numbers of components)

	# group the list by halves / component numbers / tests
	dim(factors) <- c(2, length(nfac), length(tests))
	dimnames(factors) <- list(half = NULL, nfac = nfac, test = NULL)

	tcc <- lapply(
		setNames(seq_along(nfac), nfac), # iterating over numbers of factors,
		function(ifac) structure(
			sapply(
				1:dim(factors)[3], function(grp) # for each grouping,
					sapply( # for each mode,
						c('A','B'), function(mode) # TCC for this half and mode
							diag(congru(
								factors[[1,ifac,grp]][[mode]],
								factors[[2,ifac,grp]][[mode]]
							)),
							# only for matching components
						simplify = 'array'
					),
				simplify = 'array'
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
	print(sapply(x$tcc, min))
	invisible(x)
}

shtccplot <- function(
	x, xlab = 'Number of components', ylab = 'Minimum TCC between halves', ...
) {
	# concatenate by number of factors
	df <- do.call(rbind, lapply(seq_along(x$nfac), function(i) {
		# min over mode (emission / excitation) because we use
		# the same quantity to match components
		tcc <- apply(x$tcc[[i]], c(1,3), min)
		data.frame(
			fac = as.vector(row(tcc)),
			tcc = as.vector(tcc),
			test = as.vector(col(tcc)),
			nfac = as.factor(x$nfac[i])
		)
	}))
	fac <- NULL # R CMD check vs xyplot(groups = ...)
	xyplot(
		tcc ~ nfac, df, jitter.x = T, group = fac,
		xlab = xlab, ylab = ylab, ...
	)
}

shxyplot <- function(
	x, xlab = quote(lambda*", nm"), ylab = 'Factor value', as.table = T, ...
) {
	df <- do.call(rbind, Map(
		function(x, test, half)
			do.call(rbind, lapply(1:ncol(x$A), function(i) cbind(
				rbind(
					data.frame(
						wavelength = attr(attr(x, 'cube'), 'emission'),
						value = x$A[,i],
						mode = 'Emission'
					),
					data.frame(
						wavelength = attr(attr(x, 'cube'), 'excitation'),
						value = x$B[,i],
						mode = 'Excitation'
					)
				),
				fac = as.character(i),
				nfac = as.character(ncol(x$A)),
				test = test,
				half = half
			))),
		x$factors, slice.index(x$factors, 3), slice.index(x$factors, 1)
	))
	fac <- test <- half <- NULL # R CMD check vs xyplot(groups = ...)
	xyplot(
		value ~ wavelength | mode + nfac, df, groups = paste(fac, test, half),
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
