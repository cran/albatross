library(albatross)

x <- feem(matrix(1:42, 6), 1:6, 1:7)
feemscale(x)

y <- replicate(2, x, FALSE)
feemscale(y)

z <- feemcube(y, TRUE)
feemscale(z)

z[1, 1, ] <- NA
feemscale(z, multiway::sumsq, na.rm = TRUE)
