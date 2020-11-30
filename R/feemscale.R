feemscale <- function(x, ...) UseMethod('feemscale')

feemscale.list <- feemscale.feemcube <- function(x, ..., progress = FALSE)
	cubeapply(x, feemscale, ..., progress = progress)

feemscale.feem <- function(x, norm = sd, remember = TRUE, ...) {
	factor <- norm(x, ...)
	structure(
		x / factor,
		scale = attr(x, 'scale') * if (remember) factor else 1
	)
}
