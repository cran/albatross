# Scaling is the simplest of preprocessing steps we do, but we don't
# have a notion of preprocessing pipeline, so any of those steps gets a
# separate function.
# In particular, scaling factors are remembered by default in order to
# automatically "undo" scaling of the scores.

feemscale <- function(x, ...) UseMethod('feemscale')

feemscale.list <- feemscale.feemcube <- function(x, ..., progress = FALSE)
	cubeapply(x, feemscale, ..., progress = progress)

feemscale.feem <- function(x, norm = sd, remember = TRUE, ...) {
	factor <- norm(x, ...)
	stopifnot(`Scale function returned NA: does it need na.rm = TRUE?` = !is.na(factor))
	structure(
		x / factor,
		scale = attr(x, 'scale') * if (remember) factor else 1
	)
}
