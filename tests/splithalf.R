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
for (pair in coef(sh)$subset)
	stopifnot(length(intersect(pair[[1]], pair[[2]])) == 0)

# must handle progress argument
(sh <- feemsplithalf(cube, 2, random = 1, progress = FALSE))

# must work correctly when there's only one pair of splits
coef(feemsplithalf(cube, 2, splits = 2))

# must work with a groups of length() == splits
groups <- c(rep(1, 4), rep(2, 8))
for (pair in coef(feemsplithalf(
	cube, 2, splits = 4, groups = groups
))$splits) stopifnot(
	# NB: for odd numbers of samples results are less strict but close
	table(groups[pair[[1]]]) * 2 == table(groups),
	table(groups[pair[[2]]]) * 2 == table(groups)
)

# must be able to create halves from a group of length() == 2
groups <- c(rep(1, 2), rep(2, 10))
for (pair in coef(feemsplithalf(
	cube, 2, random = 3, groups = groups
))$splits) stopifnot(
	table(groups[pair[[1]]]) * 2 == table(groups),
	table(groups[pair[[2]]]) * 2 == table(groups)
)

# must understand lists of factors
groups1 <- list(
	c(rep(1, 8), rep(2, 4)),
	c(rep(1, 4), rep(2, 8))
)
groups2 <- c(rep(1, 4), rep(2, 4), rep(3, 4))
for (pair in coef(feemsplithalf(
	cube, 2, splits = 2, groups = groups1
))$splits) stopifnot(
	table(groups2[pair[[1]]]) * 2 == table(groups2),
	table(groups2[pair[[2]]]) * 2 == table(groups2)
)

# must return the correct columns
stopifnot(colnames(coef(sh, 'tcc')) == c(
	'factor', 'tcc', 'test', 'subset', 'nfac'
))
stopifnot(colnames(coef(sh, 'factors')) == c(
	'wavelength', 'value', 'factor', 'mode',
	'nfac', 'test', 'half', 'subset'
))

# #fac, #test, Nfac must identify a point in df of TCCs
stopifnot(1 == aggregate(
	coef(sh, 'tcc')[,'tcc',drop = FALSE],
	coef(sh, 'tcc')[,c(
		'factor', 'nfac', 'test'
	)],
	FUN = length
)$tcc)

# wavelength + #fac + mode + Nfac + #test + #half must identify point
stopifnot(1 == aggregate(
	coef(sh, 'factors')[,'value',drop = FALSE],
	coef(sh, 'factors')[,c(
		'wavelength', 'factor', 'mode', 'nfac', 'test', 'half'
	)],
	FUN = length
)$value)
