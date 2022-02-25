library(albatross)
data(feems)
library(parallel)

z <- feemscatter(feems$a, rep(25, 4))
stopifnot(inherits(z, 'feem'))

feemscatter(feems$a, rep(25, 4), Raman.shift = 4200)
feemscatter(feems$a, rep(25, 4), 'pchip', 35)
feemscatter(feems$a, rep(25, 4), 'loess', NA)
feemscatter(
	feems$a[seq(1, nrow(feems$a), 2), seq(1, ncol(feems$a), 2)],
	rep(25, 4), 'kriging', NA, type = 'simple'
)
feemscatter(
	feems$a, rep(25, 4), 'w',
	d = c(1, 2), lambda = c(2e-1, 1e-2), logscale = NA
)
feemscatter(
	feems$a, rep(25, 4), 'w',
	d = c(1, 2), lambda = c(2e-1, 1e-2), logscale = 1e-3
)
feemscatter(feems$a, rep(-1, 4), add.zeroes = NA)
feemscatter(
	feems$a[seq(nrow(feems$a)*.75, nrow(feems$a)),19:25],
	rep(15, 4), 'o'
)
z <- feemscatter(feems, rep(20, 4))
stopifnot(
	inherits(z, 'list'), names(z) == names(feems)
)
(z <- feemscatter(feemcube(feems, TRUE), rep(20, 4)))
stopifnot(
	inherits(z, 'feemcube'), dimnames(z)[[3]] == names(feems)
)

cl <- makeCluster(1)
feemscatter(feems, rep(20, 4), cl = cl, progress = FALSE)
feemscatter(feems, rep(20, 4), cl = cl, progress = TRUE)
stopCluster(cl)

# pchip must fall back to linear interpolation when it doesn't have more
# than two points
z <- feem(
	matrix(c(
		1,  NA, 1,
		NA, NA, 1,
		NA,  1, 1
	), 3, byrow = TRUE),
	c(310, 320, 330), c(300, 310, 320)
)
feemscatter(z, rep(20, 4), 'pchip')
