# abs2list: convert absorbance data into a list of 2-column matrices
abs2list <- function(x) UseMethod('abs2list')
abs2list.matrix <- abs2list.data.frame <- function(x) lapply(
	setNames(2:ncol(x), colnames(x)[-1]), function(n)
		cbind(x[, 1], x[, n])
)
abs2list.list <- identity

# arrange: return x[n] if names match or length(x) == length(n)
arrange <- function(x, n, m) if (
	!is.null(names(x)) && !is.null(n) &&
	!anyDuplicated(n) && all(n %in% names(x))
) {
	x[n]
} else if (length(x) == if (is.null(n)) m else length(n)) {
	x
} else {
	stop(
		'Need either exactly N absorbance spectra or named ',
		'absorbance spectra exactly matching (unique) names of samples'
	)
}

feemife <- function(x, ...) UseMethod('feemife')

feemife.list <- function(x, absorbance, abs.path, ...) {
	if (missing(abs.path)) abs.path <- rep(1, length(x))
	stopifnot(
		length(list(...)) == 0,
		length(abs.path) == length(x)
	)
	Map(
		feemife, x, arrange(abs2list(absorbance), names(x), length(x)),
		arrange(abs.path, names(x), length(x))
	)
}

feemife.feemcube <- function(x, absorbance, abs.path, ...) {
	if (missing(abs.path)) abs.path <- rep(1, dim(x)[3])
	stopifnot(
		length(list(...)) == 0,
		length(abs.path) == dim(x)[3]
	)
	feemcube(Map(
		feemife, as.list(x),
		arrange(abs2list(absorbance), dimnames(x)[[3]], dim(x)[3]),
		arrange(abs.path, dimnames(x)[[3]], dim(x)[3])
	), TRUE)
}

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

