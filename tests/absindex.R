library(albatross)
data(feems)

absindex(absorp, abs.path = rep(2, length(absorp)))
absindex(absorp, unit = 'm^-1')
absindex(absorp, out.slope.nrmse = TRUE)
absindex(absorp, out.A = NULL, out.a = NULL, out.a.ratio = NULL, out.slope = NULL, out.slope.ratio = NULL, out.slope.nrmse = TRUE)
