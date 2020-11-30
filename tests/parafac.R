library(albatross)
data(feems)

# absolute minimum required; skip IFE correction
cube <- feemscale(
	feemscatter(
		feemcube(feems, FALSE)[(1:60)*3, (1:18)*3, ],
		rep(24, 4), 'pchip'
	),
	na.rm = T
)

# check subsetting
subs <- c(1, 2, 4, 8)
# checking values of components would be hard,
# so instead we check that dimnames are consistent
stopifnot(all.equal(
	dimnames(cube[,,subs]),
	# resid() uses fitted(), both need to account for subset
	dimnames(resid(feemparafac(cube, nfac = 3, subset = subs)))
))

factors <- feemparafac(cube, nfac = 3, const = rep('nonneg', 3))

# must return wrapped multiway::parafac object
stopifnot(inherits(factors, 'feemparafac'), inherits(factors, 'parafac'))

# need these methods
fitted(factors)
resid(factors)

# should still work when subset is missing (though it exists by default)
attr(factors, 'subset') <- NULL

# fitted / residuals must be of the same kind as original cube
stopifnot(is.null(attr.all.equal(cube, fitted(factors))))
stopifnot(is.null(attr.all.equal(cube, residuals(factors))))

# check environment access
env <- new.env(parent = emptyenv())
env$blablabla <- cube
factors <- feemparafac(
	'blablabla', nfac = 3, const = rep('nonneg', 3), envir = env
)
fitted(factors)
resid(factors)
