library(albatross)
data(feems)
cube <- feemscale(
	feemscatter(feemcube(feems, FALSE), rep(24, 4))[1:30*6, 1:9*6,],
	na.rm = TRUE
)
# check that the progress argument works
jk <- feemjackknife(cube, nfac = 2, const = rep('nonneg', 3), progress = FALSE)
