library(albatross)
library(tools) # assert*

sortfeem <- function(x)
	x[order(attr(x, 'emission')), order(attr(x, 'excitation'))]

# round-trip a matrix via data.frame
(z <- feem(matrix(1:40, ncol = 8), 66 + 1:5, 99 + 1:8))
stopifnot(
	all.equal(
		sortfeem(z),
		sortfeem(feem(as.data.frame(feem(as.data.frame(z)))))
	)
)

# extraction operator must return FEEM objects unless dropping dimensions
stopifnot(inherits(z[2:4, 3:7], 'feem'))

# replacement operator must verify wavelengths when assigning from FEEM object
assertError(
	z[1:2, 3:4] <- feem(matrix(1:4, 2), 1:2, 1:2)
)
# but the check is disabled if we use 1- or 0- argument form of [<-
z[1:2,3:4][] <- feem(matrix(1:4, 2), 1:2, 1:2)

# replacement operator must warn about scale differences
assertWarning(
	z[2:3, 4:5] <- feem(matrix(1:4, 2), 66 + 2:3, 99 + 4:5, 2),
	verbose = TRUE
)

# sub-assignment with unset subscripts must work
z[] <- z
z[,] <- z

# must keep wavelengths when subscripting with dimnames
stopifnot(all.equal(
	z,
	z[dimnames(z)[[1]], dimnames(z)[[2]]]
))
zz <- z
zz[dimnames(z)[[1]], dimnames(z)[[2]]] <- z[]
stopifnot(all.equal(z, zz))

sortdf <- function(x) x[do.call(order, x),]

# round-trip a sparse data.frame via feem
(z <- data.frame(emission = 1:10, excitation = 21:30, intensity = 31:40))
stopifnot(
	all.equal(
		sortdf(z),
		sortdf(as.data.frame(feem(as.data.frame(feem(z)))))
	)
)

# must not allow matrices of correct length:
# other code assumes that emission and excitation doesn't have dimensions
assertError(feem(matrix(1:4, 2), matrix(1:2), 1:2), verbose = TRUE)
assertError(feem(matrix(1:4, 2), 1:2, matrix(1:2, 1)), verbose = TRUE)
assertError(feem(matrix(1:4, 2), letters[1:2], matrix(1:2)), verbose = TRUE)
# non-scalar/non-numeric scale
assertError(feem(matrix(1:4, 2), 1:2, 1:2, 1:2), verbose = TRUE)
assertError(feem(matrix(1:4, 2), 1:2, 1:2, matrix(1)), verbose = TRUE)
assertError(feem(matrix(1:4, 2), 1:2, 1:2, 'x'), verbose = TRUE)
# non-numeric x
# NB: character x gets dispatched to feem.character and fails for a
# different reason
assertError(feem(as.complex(matrix(1:4, 2)), 1:2, 1:2), verbose = TRUE)
