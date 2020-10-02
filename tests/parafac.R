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

factors <- feemparafac(cube, nfac = 3, const = rep('nonneg', 3))

# must return wrapped multiway::parafac object
stopifnot(inherits(factors, 'feemparafac'), inherits(factors, 'parafac'))

# need these methods
fitted(factors)
resid(factors)

# fitted / residuals must be of the same kind as original cube
stopifnot(is.null(attr.all.equal(cube, fitted(factors))))
stopifnot(is.null(attr.all.equal(cube, residuals(factors))))
