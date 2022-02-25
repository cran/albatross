library(albatross)
data(feems)
cube <- feemscale(feemscatter(cube, rep(24, 4)), na.rm = TRUE)
# check that the progress argument works
jk <- feemjackknife(cube, nfac = 2, const = rep('nonneg', 3), progress = FALSE)

cols <- list(
	estimations = c('loading', 'mode', 'wavelength', 'factor', 'omitted'),
	RIP = c('msq.resid', 'Emission', 'Excitation', 'omitted'),
	IMP = c('score.overall', 'score.predicted', 'factor', 'omitted')
)
for (n in names(cols)) stopifnot(cols[[n]] == colnames(coef(jk, n)))

# 1-component jack-knife should also work
jk1 <- feemjackknife(cube, nfac = 1, const = rep('nonneg', 3))
stopifnot(is.matrix(jk1$leaveone[[1]]$A))

stopifnot(all.equal(cube, feemcube(jk)))
