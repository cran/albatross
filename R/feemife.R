# abs2list: convert absorbance data into a list of 2-column matrices
abs2list <- function(x) UseMethod('abs2list')
abs2list.matrix <- abs2list.data.frame <- function(x) lapply(
	setNames(2:ncol(x), colnames(x)[-1]), function(n)
		cbind(x[, 1], x[, n])
)
abs2list.list <- identity

# arrange: return x[n] if names match
# or just x if x is not named but the lengths match
arrange <- function(x, n, m) if (
	!is.null(names(x)) && !is.null(n) &&
	!anyDuplicated(n) && all(n %in% names(x))
) {
	x[n]
} else if (length(x) == m && (is.null(names(x)) || is.null(n))) {
	x
} else {
	stop(
		'Names of ', deparse(substitute(x)), ' must exactly match ',
		'(unique) names of samples. Otherwise, length must match the ',
		'number of samples while some of the names are missing.'
	)
}

feemife <- function(x, ...) UseMethod('feemife')

feemife.list <- function(x, absorbance, abs.path, ..., progress = FALSE) {
	if (missing(abs.path)) abs.path <- rep(1, length(x))
	stopifnot(
		length(list(...)) == 0,
		length(abs.path) == length(x)
	)
	cubeapply(
		x, feemife, arrange(abs2list(absorbance), names(x), length(x)),
		arrange(abs.path, names(x), length(x)),
		progress = progress, .recycle = TRUE
	)
}

feemife.feemcube <- function(x, absorbance, abs.path, ..., progress = FALSE)
	cubeapply(x, feemife, absorbance, abs.path, ..., progress = progress)

feemife.feem <- function(x, absorbance, abs.path = 1, ...) {
	stopifnot(
		length(list(...)) == 0,
		min(absorbance[,1]) <= min(attr(x, 'emission')),
		max(absorbance[,1]) >= max(attr(x, 'emission')),
		min(absorbance[,1]) <= min(attr(x, 'excitation')),
		max(absorbance[,1]) >= min(attr(x, 'excitation')),
		ncol(absorbance) == 2
	)
	od <- splinefun(absorbance)
	x * outer(
		attr(x, 'emission'), attr(x, 'excitation'),
		function(em, ex) 10^((od(em) + od(ex)) / (2 * abs.path))
	)
}

