library(albatross)
data(feems)
f <- feemflame(
	feemscale(cube), ffac = 3, sfac = 1,
	control.cmf = list(maxit = 4), control.parafac = list(maxit = 2)
)
stopifnot(inherits(fitted(f), 'feemcube'))
stopifnot(inherits(resid(f), 'feemcube'))
stopifnot(inherits(feemcube(f), 'feemcube'))
