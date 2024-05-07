library(albatross)
data(feems)

# absolute minimum required; skip IFE correction
cube <- feemscale(feemscatter(cube, rep(24, 4), 'pchip'), na.rm = T)

# check subsetting
subs <- c(1, 2, 4, 8)
pf <- feemparafac(cube, nfac = 3, subset = subs)
# checking values of components would be hard,
# so instead we check that dimnames are consistent
stopifnot(all.equal(
	dimnames(cube[,,subs]),
	# resid() uses fitted(), both need to account for subset
	dimnames(resid(pf))
))
# also check the equivalence of the cube
stopifnot(all.equal(cube[,,subs], feemcube(pf)))

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
# residuals() sets an additional class of feem.resid which only affects plots
stopifnot(is.null(attr.all.equal(
	# Used to call class<- here, but on R-3.3.0 it overwrites the argument
	structure(cube,
		class = c('feem.resid', class(cube))
	),
	residuals(factors)
)))

# check environment access
env <- new.env(parent = emptyenv())
env$blablabla <- cube
factors <- feemparafac(
	'blablabla', nfac = 3, const = rep('nonneg', 3),
	envir = env
)
fitted(factors)
resid(factors)
stopifnot(all.equal(cube, feemcube(factors)))

# dimnames should be assigned
stopifnot(
	dimnames(cube)[[1]] == rownames(factors$A),
	dimnames(cube)[[2]] == rownames(factors$B),
	dimnames(cube)[[3]] == rownames(factors$C)
)

# coef must return data.frames or lists with correct contents
coefnames <- list(
	emission = c('wavelength', 'value', 'factor'),
	excitation = c('wavelength', 'value', 'factor'),
	samples = c('sample', 'value', 'factor'),
	scores = c('sample', 'value', 'factor'),
	loadings = c('wavelength', 'value', 'factor', 'mode')
)
for (n in names(coefnames))
	stopifnot(all.equal(colnames(coef(factors, n)), coefnames[[n]]))
allnames <- c('emission', 'excitation', 'samples')
stopifnot(all.equal(names(coef(factors, 'all')), allnames))
for (n in allnames)
	stopifnot(all.equal(
		colnames(coef(factors, 'all')[[n]]),
		coefnames[[n]]
	))

# reorder, rescale without like= should still work
stopifnot(
	inherits(reorder(pf, 3:1), 'feemparafac'),
	inherits(rescale(pf, mode = 'C', newscale = 2, absorb = 'A'), 'feemparafac')
)

stopifnot(
	inherits(
		# rescale with not all arguments specified should still work
		rescale(pf, 'A', absorb = 'C'),
		'feemparafac'
	)
)

# reorder(like=...) must order the components to match
stopifnot(
	identical(
		reorder(
			reorder(pf, 3:1),
			like = pf
		),
		pf
	)
)

# round-trip reorder() and rescale() in both forms
all.equal(
	factors,
	rescale(
		reorder(
			reorder(
				rescale(factors, mode = 'C', newscale = 2, absorb = 'A'),
				3:1
			),
			like = factors
		),
		like = factors
	),
	tolerance = .Machine$double.eps^(1/3)
)

pf.1 <- feemparafac(cube, nfac = 1, subset = subs)
tools::assertError(reorder(pf.1, neworder = 2))
pf.1 <- reorder(pf.1, neworder = 1)
stopifnot(is.matrix(pf.1$A))
