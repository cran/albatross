library(albatross)
library(tools)

# .integ must work about as well as linear interpolation
data(feems)
x <- albatross:::.extract(
	feems$a, c(290, 310), 254, TRUE, c(1,1), FALSE
)
(trapz <- albatross:::.integ(x))
with(
	integrate(
		approxfun(attr(x, 'emission'), as.vector(x)),
		min(attr(x, 'emission')), max(attr(x, 'emission'))
	),
	if (abs(trapz - value) >= abs.error) stop(paste(
		'Trapezoid value', trapz, '- approxfun value', value,
		'=', trapz - value, '|>=| error estimate', abs.error
	))
)
(ztrap <- albatross:::.integ(x[nrow(x):1,, drop = FALSE]))
# .integ must handle unsorted data
stopifnot(all.equal(trapz, ztrap))

# exact wavelength extraction
x <- feem(matrix(1:4, 2), c(300, 450), c(310, 320))
# - must return NAs for points out of range
(y <- albatross:::.extract(x, c(380, 430), 310, FALSE, c(60, 60), FALSE))
stopifnot(any(is.na(y)))
# - if multiple points match due to tolerance being too high...
assertWarning(
	y <- albatross:::.extract(x, c(380, 430), 310, FALSE, c(70, 70), FALSE),
	verbose = TRUE
)
(y) # => must warn and set them to NA anyway
stopifnot(any(is.na(y)))
# - must be able to interpolate if requested
(y <- albatross:::.extract(x, c(380, 430), 310, FALSE, c(2, 2), 'whittaker', d = 1, lambda = 1e-3))
stopifnot(!is.na(y))
# - but not if it's out of range
(y <- albatross:::.extract(x, 290, 310, FALSE, c(2, 2), 'whittaker', d = 1, lambda = 1e-3))
stopifnot(is.na(y))
# - must interpolate NAs present in the spectrum if enabled
xx <- replace(x, matrix(c(1,1), 1), NA)
(y <- albatross:::.extract(x, 300, 310, FALSE, c(2, 2), 'whittaker', d = 1, lambda = 1e-3))
stopifnot(!is.na(y))

# range extraction
x <- feem(matrix(1:12, 4), seq(300, 330, len = 4), seq(270, 360, len = 3))
# - no duplicates even if border touches grid points
(y <- albatross:::.extract(x, c(309, 321), c(314, 316), TRUE, c(2, 2), FALSE))
stopifnot(!anyDuplicated(attr(y, 'emission')), !anyDuplicated(attr(y, 'excitation')))
# - must be able to interpolate whole ranges
(y <- albatross:::.extract(x, c(305, 325), c(300, 350), TRUE, c(2, 2), 'whittaker'))
stopifnot(!is.na(y))
