library(albatross)
data(feems)

# assuming that all wavelengths in absorp are the same,
stopifnot(0 == sum(sapply(2:length(absorp), function(i)
	sum(abs(absorp[[i]][,1] - absorp[[1]][,1]))
)))
# build a matrix of absorbance spectra with named columns
absmat <- do.call(cbind, c(list(absorp[[1]][,1]), lapply(absorp, `[`,, 2)))
matplot(absmat[,1], absmat[,-1])

# simple one-spectrum correction
plot(feemcube(list(
	before = feemscatter(feems$a, rep(24, 4)),
	after = feemscatter(feemife(feems$a, absorp$a), rep(24, 4))
), TRUE))

# list-list correction
# corr[[i]] should match the result of IFE-correcting
# feems[[name[i]]] with absorp[[name[i]]
check.ife <- function(corr, name) stopifnot(simplify2array(Map(
	function(x, n) all.equal(feemife(feems[[n]], absorp[[n]]), x),
	corr, name
)))
check.ife(feemife(feems, absorp), names(feems))
subn <- c('a','e','j','g')
check.ife(feemife(feems[subn], absorp), subn)
check.ife(feemife(unname(feems[subn]), unname(absorp[subn])), subn)

# list-matrix correction
check.ife(feemife(feems, absmat), names(feems))
check.ife(feemife(feems[subn], absmat), subn)
check.ife(
	feemife(
		unname(feems[subn]),
		unname(absmat[, c(1, match(subn, colnames(absmat)))])
	), subn
)

# feemcube-list correction
fcube <- feemcube(feems, TRUE)
check.ife2 <- function(corr, name) stopifnot(sapply(
	seq_along(name), function(i) {
		reference <- feemife(feems[[name[i]]], absorp[[name[i]]])
		all.equal(
			reference,
			corr[
				match(attr(reference, 'emission'), attr(corr, 'emission')),
				match(attr(reference, 'excitation'), attr(corr, 'excitation')),
				i
			]
		)
	}
))
check.ife2(feemife(fcube, absorp), dimnames(fcube)[[3]])
subn <- c('a','e','j','g')
check.ife2(
	feemife(fcube[,,subn], absorp), subn
)
check.ife2(
	feemife(unname(fcube[,,subn]), unname(absorp[subn])), subn
)

# feemcube-matrix correction
check.ife2(feemife(fcube, absmat), dimnames(fcube)[[3]])
check.ife2(
	feemife(fcube[,,subn], absmat), subn
)
check.ife2(
	feemife(
		unname(fcube[,,subn]),
		unname(absmat[, c(1, match(subn, colnames(absmat)))])
	), subn
)

# plot some correction coefficients
plot(z <- feemcube(lapply(setNames(nm = c(1, 2, 5, 10)), function (i)
	feemscatter(feemife(feems$a, absorp$a, i), rep(24, 4))
), TRUE))
plot(feemcube(list(
	corr.mat = z[,,1] / feems$a,
	one.to.ten = z[,,1] / z[,,4]
), TRUE))
