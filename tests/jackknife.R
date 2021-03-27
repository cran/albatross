library(albatross)
data(feems)
cube <- feemscale(
	feemscatter(feemcube(feems, FALSE), rep(24, 4))[1:30*6, 1:9*6,],
	na.rm = TRUE
)
# check that the progress argument works
jk <- feemjackknife(cube, nfac = 2, const = rep('nonneg', 3), progress = FALSE)

cols <- list(
	estimations = c('loading', 'mode', 'wavelength', 'factor', 'omitted'),
	RIP = c('msq.resid', 'Emission', 'Excitation', 'omitted'),
	IMP = c('score.overall', 'score.predicted', 'factor', 'omitted')
)
for (n in names(cols)) stopifnot(cols[[n]] == colnames(coef(jk, n)))
