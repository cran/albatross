library(albatross)
data(feems)
# must be able to specify methods
# must be able to disable progress bar
# must be able to handle lists of FEEMs
feemgrid(
	feemscatter(feems, rep(0, 4), progress = FALSE),
	method = 'pchip', progress = FALSE
)

# must be able to handle FEEM cubes
feemgrid(
	feemscatter(feemcube(feems, TRUE), rep(0, 4), progress = FALSE),
	method = 'pchip'
)

# when requested only one wavelength along either axis, should return
# a FEEM, not a named vector
# a.k.a. yet another drop = FALSE bug
x <- feems$a
stopifnot(inherits(
	feemgrid(x, emission = attr(x, 'emission'), excitation = 254),
	'feem'
))
stopifnot(inherits(
	feemgrid(x, emission = 300.5, excitation = attr(x, 'excitation')),
	'feem'
))
