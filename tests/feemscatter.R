library(albatross)
data(feems)

plot(feemscatter(feems$a, rep(25, 4)))
plot(feemscatter(feems$a, rep(25, 4), Raman.shift = 4200))
plot(feemscatter(feems$a, rep(25, 4), 'pchip', 35))
plot(feemscatter(feems$a, rep(25, 4), 'loess', NA))
plot(feemscatter(feems$a, rep(-1, 4), add.zeroes = NA))
plot(feemscatter(feems$a[-(1:90),19:25], rep(15, 4), 'o'))
z <- feemscatter(feems, rep(20, 4))
stopifnot(
	inherits(z, 'list'), names(z) == names(feems)
)
plot(z <- feemscatter(feemcube(feems, TRUE), rep(20, 4)))
stopifnot(
	inherits(z, 'feemcube'), dimnames(z)[[3]] == names(feems)
)
