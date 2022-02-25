library(tools)
library(albatross)
data(feems)

# assuming that all wavelengths in absorp are the same,
stopifnot(0 == sum(sapply(2:length(absorp), function(i)
	sum(abs(absorp[[i]][,1] - absorp[[1]][,1]))
)))
# build a matrix of absorbance spectra with named columns
absmat <- do.call(cbind, c(list(absorp[[1]][,1]), lapply(absorp, `[`,, 2)))

# simple one-spectrum correction
feemcube(list(
	before = feemscatter(feems$a, rep(24, 4)),
	after = feemscatter(feemife(feems$a, absorp$a), rep(24, 4))
), TRUE)

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
# correction should fail if lengths are the same but names don't match
assertError(
	feemife(feems[c('a', 'b', 'c')], absorp[c('c', 'a', 'd')]),
	verbose = TRUE
)
# should fail if names of feems are duplicated
assertError(
	feemife(feems[c('a', 'a', 'c')], absorp[c('c', 'a', 'd')]),
	verbose = TRUE
)
# should fail if no names but lengths differ
assertError(
	feemife(unname(feems[c('a', 'a', 'c')]), unname(absorp[c('c', 'a')])),
	verbose = TRUE
)

# list-matrix correction
check.ife(feemife(feems, absmat), names(feems))
check.ife(feemife(feems[subn], absmat), subn)
check.ife(
	feemife(
		unname(feems[subn]),
		unname(absmat[, c(1, match(subn, colnames(absmat)))])
	), subn
)
assertError(feemife(
	feems[c('a', 'b', 'c')],
	absmat[, c(1, match(c('c', 'a', 'd'), colnames(absmat)))]
), verbose = TRUE)

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
check.ife2(
	feemife(fcube[,,subn], absorp), subn
)
check.ife2(
	feemife(unname(fcube[,,subn]), unname(absorp[subn])), subn
)
assertError(feemife(
	fcube[,,c('a', 'b', 'c')], absorp[c('c', 'a', 'd')]
), verbose = TRUE)

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
assertError(feemife(
	fcube[,,c('a', 'b', 'c')],
	absmat[, c(1, match(c('c', 'a', 'd'), colnames(absmat)))]
), verbose = TRUE)

feemife(feems, absorp, progress = TRUE)
# arrange must work (and report errors) for paths the same way as for
# absorption spectra themselves
feemife(feems[c('b','a','c')], absorp, setNames(rep(1, 12), letters[1:12]))
stopifnot(all.equal(
	feemife(feems, absorp, 1:12),
	feemife(feems, absorp, setNames(12:1, rev(letters[1:12])))
))
assertError(feemife(
	feems[c('a','a','b','c')], absorp,
	setNames(rep(1, 12), letters[1:12])
), verbose = TRUE)
assertError(feemife(
	feems[c('a','b','c')], absorp,
	setNames(rep(1, 3), letters[1:3 + 1])
), verbose = TRUE)
assertError(
	feemife(feems[c('a','b','c')], absorp, rep(1, 4)),
	verbose = TRUE
)

# also test cubeapply error wrapping
ab <- absorp
ab$e <- ab$e[ab$e[,1] > min(attr(feems$e, 'emission')) + 10,]
(assertCondition(feemife(feems, ab), 'feem.wrapped.error', verbose = TRUE))
