library(albatross)
library(tools) # assert*

# array constructor
z <- feemcube(array(1:(11*13*3), c(11, 13, 3)), 1:11, 20 + 1:13, 2:4, names = letters[1:3])

# when not dropping dimensions, should return a cube
stopifnot(inherits(z[1:10, 1:10, 1:2], 'feemcube'))
# when dropping the sample dimension, should return a feem with correct scale
stopifnot(
	inherits(z[1:10, c(TRUE, FALSE), 1], 'feem'),
	attr(z[1:10, c(TRUE, FALSE), 1], 'scale') == 2
)

# assign from FEEM matching wavelengths and scale
z[3 + 1:3, 5 + 1:5, 1] <- feem(matrix(15:1, 3), 3 + 1:3, 20 + 5 + 1:5, 2)

# disallow assignment from non-matching wavelengths
assertError(z[,,2] <- feem(matrix(1:(11*13), 11), 1:11, 1:13, 3))
assertError(z[,,] <- structure(z, emission = attr(z, 'emission') + 10))
# no wavelength check with 1- or 0- argument form of [<-
z[,,2][] <- feem(matrix(1:(11*13), 11), 1:11, 1:13, 3)

# warn about assignment from non-matching scales
assertWarning(z[,,] <- structure(z, scales = c(1e-2, 1, 1e+3)), verbose = TRUE)
assertWarning(
	z[3 + 1:3, 5 + 1:5, 1] <- feem(matrix(15:1, 3), 3 + 1:3, 20 + 5 + 1:5, 1),
	verbose = TRUE
)

stopifnot(
	c('emission', 'excitation', 'intensity', 'sample') %in%
	colnames(as.data.frame(z[,1, 1, drop = F]))
)

z <- feemscale(z, na.rm = TRUE)
stopifnot(all.equal(z, feemcube(as.list(z), TRUE)))

# "sample" column should be factor or character, but not integer
for (cube in list(z, feemcube(unname(as.list(z)), FALSE)))
	with(as.data.frame(cube),
		stopifnot(is.factor(sample) || is.character(sample))
	)

# sub-assignment with unset indices must work
z[] <- z
z[,,] <- z

# must correctly index by dimnames
stopifnot(all.equal(
	z,
	z[dimnames(z)[[1]], dimnames(z)[[2]], dimnames(z)[[3]]]
))
zz <- z
zz[dimnames(z)[[1]], dimnames(z)[[2]], dimnames(z)[[3]]] <- z[]
stopifnot(all.equal(z, zz))

# emission, excitation, scales must be numeric vectors
# names must be atomic vector
# x must be numeric 3-way array
# sizes must match
assertError(
	feemcube(array(1:(11*13*3), c(11, 13, 3)), matrix(1:11, 1), 20 + 1:13),
	verbose = TRUE
)
assertError(
	feemcube(array(1:(11*13*3), c(11, 13, 3)), 1:11, matrix(20 + 1:13)),
	verbose = TRUE
)
assertError(
	feemcube(array(1:(11*13*3), c(11, 13, 3)), as.list(1:11), 20 + 1:13),
	verbose = TRUE
)
assertError(
	feemcube(array(1:(11*13*3), c(11, 13, 3)), 1:11, letters[1:13]),
	verbose = TRUE
)
assertError(
	feemcube(array(1:(11*13*3), c(11, 13, 3)), 1:11, 1:13, matrix(2:4)),
	verbose = TRUE
)
assertError(
	feemcube(array(1:(11*13*3), c(11, 13, 3)), 1:11, 1:13, 2:4, as.list(letters[2:4])),
	verbose = TRUE
)
assertError(
	feemcube(as.character(array(1:(11*13*3), c(11, 13, 3))), 1:11, 20 + 1:13),
	verbose = TRUE
)
