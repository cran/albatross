library(albatross)
data(feems)
# must be able to specify methods
# must be able to disable progress bar
# must be able to handle lists of FEEMs
feemgrid(
	feemscatter(feems, rep(0, 4), progress = FALSE),
	method = 'pchip', progress = FALSE
)

# must be able to handle FEEM cubes
feemgrid(
	feemscatter(feemcube(feems, TRUE), rep(0, 4), progress = FALSE),
	method = 'pchip'
)
