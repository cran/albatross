library(albatross)
data(feems)

cube <- feemscale(
	feemscatter(
		feemcube(feems, FALSE)[(1:36)*5, (1:11)*5, ],
		rep(24, 4), 'pchip'
	),
	na.rm = T
)

jk <- feemjackknife(cube, nfac = 2, const = rep('nonneg', 3))
plot(jk)
plot(jk, 'estimations')
plot(jk, 'RIP')
plot(jk, 'IMP')
