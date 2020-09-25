feemscale <- function(x, ...) UseMethod('feemscale')

feemscale.list <- function(x, ...) lapply(x, feemscale, ...)

feemscale.feemcube <- function(x, ...) feemcube(feemscale(as.list(x), ...), TRUE)

feemscale.feem <- function(x, norm = sd, remember = TRUE, ...) {
	factor <- norm(x, ...)
	structure(
		x / factor,
		scale = attr(x, 'scale') * if (remember) factor else 1
	)
}
