library(albatross)
data(feems)
cube <- feemscale(feemscatter(cube, rep(24, 4)), na.rm = TRUE)
# progress argument must work, 1-component model must work
system.time(
	jk <- feemjackknife(
		cube, nfac = 1, ctol = 1e-4, progress = FALSE
	)
)
stopifnot(is.matrix(jk$leaveone[[1]]$A))

cols <- list(
	estimations = c('loading', 'mode', 'wavelength', 'factor', 'omitted'),
	RIP = c('msq.resid', 'Emission', 'Excitation', 'omitted'),
	IMP = c('score.overall', 'score.predicted', 'factor', 'omitted')
)
for (n in names(cols)) stopifnot(cols[[n]] == colnames(coef(jk, n)))

stopifnot(all.equal(cube, feemcube(jk)))

# Also test parallel `bootparafac` with postprocessing
library(parallel)
cl <- makeCluster(2)
system.time(jk <- feemjackknife(
	cube, nfac = 1, ctol = 1e-4,
	parallel = TRUE, cl = cl, .scheduling = 'static'
))
stopCluster(cl)
plot(jk)
stopifnot(is.environment(attr(jk$leaveone[[1]], 'envir')))
