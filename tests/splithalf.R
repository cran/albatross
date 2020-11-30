library(albatross)
data(feems)

cube <- feemscale(
	feemscatter(
		feemcube(feems, FALSE)[(1:36)*5, (1:11)*5, ],
		rep(24, 4), 'pchip'
	),
	na.rm = T
)

# must make sense for only one number of factors / one split
(sh <- feemsplithalf(cube, 2, random = 1))
stopifnot(inherits(sh, 'feemsplithalf'))

# must handle non-even numbers of samples
(sh <- feemsplithalf(cube[,,1:11], 2:3, splits = 4, const = rep('nonneg', 3)))

# must not compare splits containing same samples
for (nfac in seq_len(dim(sh$factors)[2]))
	lapply(
		apply(sh$factors[,nfac,], 2, lapply, attr, 'subset'), # extract subsets
		function(p)
			stopifnot(length(intersect(p[[1]], p[[2]])) == 0)
	)

# must handle progress argument
(sh <- feemsplithalf(cube, 2, random = 1, progress = FALSE))
