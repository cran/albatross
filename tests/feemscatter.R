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
setDefaultCluster(cl)
feemscatter(feems, rep(20, 4), cl = NULL, progress = TRUE)
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
# this previously failed because t() silently produced botched feems
feemscatter(feems$a, rep(25, 4), 'pchip', by = 'both')

z <- feem(
	matrix(0, 401, 69),
	300:700, seq(210, 550, 5)
)

for (w in list(
	list(
		# Rayleigh +/- 10 nm
		10,
		quote(abs(emission - excitation) < 10)
	),
	list(
		# Rayleigh & Raman +/- 10 nm
		c(10, 10),
		quote(
			abs(emission - excitation) < 10 |
			abs(emission - 1/(1/excitation - 3400/1e+07)) < 10
		)
	),
	list(
		# Rayleigh -20 +10, Raman +/- 10 nm
		list(c(20, 10), 10),
		quote(
			((emission - excitation > -20) & (emission - excitation < 10)) |
			abs(emission - 1/(1/excitation - 3400/1e+07)) < 10
		)
	),
	list(
		# orders 1&2 +/- 10 nm; auto-scale 2nd and 3rd order; 3rd order Rayleigh +/- 5 nm but no Raman
		c(10, 10, 10, 10, 5),
		quote(
			abs(emission - excitation) < 10 |
			abs(emission - 1/(1/excitation - 3400/1e+07)) < 10 |
			abs(emission/2 - excitation) < 10 |
			abs(emission/2 - 1/(1/excitation - 3400/1e+07)) < 10 |
			abs(emission/3 - excitation) < 5
		)
	)
)) {
	d <- merge(
		as.data.frame(z),
		as.data.frame(feemscatter(z, w[[1]], Raman.shift = 3400)),
		by = c('emission', 'excitation'),
		all = TRUE
	)
	mask <- eval(w[[2]], d)
	stopifnot(
		all(is.na(d$intensity.y[mask])),
		all(!is.na(d$intensity.y[!mask]))
	)
}
