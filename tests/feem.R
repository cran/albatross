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

# replacement operator must warn about scale differences
assertWarning(
	z[2:3, 4:5] <- feem(matrix(1:4, 2), 66 + 2:3, 99 + 4:5, 2)
)

sortdf <- function(x) x[do.call(order, x),]

# round-trip a sparse data.frame via feem
(z <- data.frame(emission = 1:10, excitation = 21:30, intensity = 31:40))
stopifnot(
	all.equal(
		sortdf(z),
		sortdf(as.data.frame(feem(as.data.frame(feem(z)))))
	)
)
