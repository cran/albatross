library(albatross)
data(feems)

z <- feemscatter(feems$a, rep(25, 4))
stopifnot(inherits(z, 'feem'))

feemscatter(feems$a, rep(25, 4), Raman.shift = 4200)
feemscatter(feems$a, rep(25, 4), 'pchip', 35)
feemscatter(feems$a, rep(25, 4), 'loess', NA)
feemscatter(feems$a, rep(-1, 4), add.zeroes = NA)
feemscatter(feems$a[-(1:90),19:25], rep(15, 4), 'o')
z <- feemscatter(feems, rep(20, 4))
stopifnot(
	inherits(z, 'list'), names(z) == names(feems)
)
(z <- feemscatter(feemcube(feems, TRUE), rep(20, 4)))
stopifnot(
	inherits(z, 'feemcube'), dimnames(z)[[3]] == names(feems)
)
