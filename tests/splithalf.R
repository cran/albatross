library(albatross)
library(tools)
data(feems)

cube <- feemscale(feemscatter(cube, rep(24, 4), 'pchip'), na.rm = TRUE)

# must make sense for only one number of factors / one split
(sh <- feemsplithalf(cube, 2, random = 1, ctol = 1e-4))
stopifnot(inherits(sh, 'feemsplithalf'))

# must handle non-even numbers of samples
(sh <- feemsplithalf(cube[,,1:11], 2, splits = 4, ctol = 1e-4))
# feemsplithalf(splits=s) makes choose(s, s/2) halves then combines them
# resulting in choose(s, s/2)/2 comparisons
stopifnot(dim(sh$factors)[3] == choose(4, 2)/2)
# must be able to limit the number of comparisons
(sh <- feemsplithalf(cube, 2, splits = c(10, 5), ctol = 1e-4))
stopifnot(dim(sh$factors)[3] == 5)

# must not compare splits containing same samples
for (pair in coef(sh)$subset)
	stopifnot(length(intersect(pair[[1]], pair[[2]])) == 0)

# must handle progress argument
(sh <- feemsplithalf(cube, 2, random = 1, progress = FALSE, ctol = 1e-4))

# must work correctly when there's only one pair of splits
coef(feemsplithalf(cube, 2, splits = 2, ctol = 1e-4))

# must work with groups of length() == splits
groups <- c(rep(1, 4), rep(2, 8))
for (pair in coef(feemsplithalf(
	cube, 2, splits = 4, groups = groups, ctol = 1e-4
))$splits) stopifnot(
	# NB: for odd numbers of samples results are less strict but close
	table(groups[pair[[1]]]) * 2 == table(groups),
	table(groups[pair[[2]]]) * 2 == table(groups)
)

# must be able to create halves from a group of length() == 2
groups <- c(rep(1, 2), rep(2, 10))
for (pair in coef(feemsplithalf(
	cube, 2, random = 3, groups = groups, ctol = 1e-4
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
	cube, 2, splits = 2, groups = groups1, ctol = 1e-4
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
	tcc ~ factor + nfac + test, coef(sh, 'tcc'),
	FUN = length
)$tcc)

# wavelength + #fac + mode + Nfac + #test + #half must identify point
stopifnot(1 == aggregate(
	value ~ wavelength + factor + mode + nfac + test + half,
	coef(sh, 'factors'),
	FUN = length
)$value)

fixed <- list(list(
	1:(round(dim(cube)[3]/2)),
	(round(dim(cube)[3]/2)+1):dim(cube)[3]
))
sh <- feemsplithalf(cube, 1, fixed = fixed, ctol = 1e-4)
# NOTE: testing one-factor model so that coef(sh) will have as many rows
# as elements in `fixed`
# NOTE: unclass() removes the "AsIs" class from the data.frame column
stopifnot(all.equal(fixed, unclass(coef(sh)$subset)))
# testing intersecting splits is an error
fixed[[1]][[2]][1] <- fixed[[1]][[1]][1]
assertError(feemsplithalf(cube, 2, fixed = fixed), verbose = TRUE)
# providing groups at the same time as fixed splits, or asking for
# split-combine or random splits, is an error
assertError(
	feemsplithalf(cube, 2, fixed = fixed, groups = groups), verbose = TRUE
)
assertError(
	feemsplithalf(cube, 2, fixed = fixed, splits = 2), verbose = TRUE
)
assertError(
	feemsplithalf(cube, 2, fixed = fixed, random = 2), verbose = TRUE
)
# asking for both split-combine and random halves is an error
assertError(feemsplithalf(cube, 2, splits = 2, random = 2), verbose = TRUE)

stopifnot(all.equal(cube, feemcube(sh)))

# exercise some plot types not used otherwise
plot(sh, 'factors')
# "subset" takes some care to forward NSE bits properly
plot(sh, 'bandfactors', subset = nfac == 1 & mode == 'Emission')

# Parallel `bootparafac` is a whole separate thing now
library(parallel)
cl <- makeCluster(2)
sh <- feemsplithalf(cube, 2, random = 2, parallel = TRUE, cl = cl, ctol = 1e-4)
stopCluster(cl)
head(coef(sh, 'aggtcc'))
head(coef(sh, 'bandfactors'))
plot(sh, 'aggtcc')
stopifnot(is.environment(attr(sh$factors[[1]], 'envir')))
