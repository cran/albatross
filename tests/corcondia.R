library(albatross)
data(feems)
# Need a cube without missing data for mode='pinv'
cube <- feemscale(feemscatter(cube, c(20, 14), 'whittaker'), na.rm = TRUE)
feemcorcondia(feemparafac(cube, nfac = 3, ctol = 1e-4))
