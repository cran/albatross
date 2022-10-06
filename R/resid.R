`[.feem.resid` <- function(x, ...) {
	x <- NextMethod()
	if (inherits(x, c('feem', 'feemcube')))
		class(x) <- c('feem.resid', class(x))
	x
}

plot.feem.resid <- function(x, ..., at, col.regions = diverging.colours(256)) {
	if (missing(at)) {
		r <- max(abs(range(x, na.rm = TRUE)))
		at <- seq(-r, r, length.out = 255)
	}
	NextMethod(x, ..., at = at, col.regions = col.regions)
}
