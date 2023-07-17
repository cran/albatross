# Inner-filter effect correction and matching absorbance spectra to
# fluorescence spectra happens here.

# abs2list: convert absorbance data into a list of 2-column matrices
abs2list <- function(x) UseMethod('abs2list')
abs2list.matrix <- abs2list.data.frame <- function(x) lapply(
	setNames(2:ncol(x), colnames(x)[-1]), function(n)
		cbind(x[, 1], x[, n])
)
abs2list.list <- identity

# arrange: return x[n] if names match
# or just x if either is not named but the lengths match
arrange <- function(x, n, m, kind) if (
	!is.null(names(x)) && !is.null(n) &&
	!anyDuplicated(n) && all(n %in% names(x))
) {
	x[n]
} else if (length(x) == m && (is.null(names(x)) || is.null(n))) {
	x
} else stop(
	'Failed to look up ', kind, ' corresponding to the ',
	'spectra:\n',
	if (is.null(names(x)) || is.null(n)) c(
		if (is.null(names(x))) ' * No names for ', kind, '\n',
		if (is.null(n)) ' * No names for the spectra\n',
		if (length(x) != m) paste(
			' * Have', m, ' spectra but', length(x),
			kind, '\n'
		)
	) else c(
		if (anyDuplicated(n)) paste(
			' * The following spectra have duplicated names:',
			paste(dQuote(unique(n[duplicated(n)])), collapse = ', '), '\n'
		),
		if (any(!n %in% names(x))) paste0(
			' * The following spectra don\'t have ',
			'corresponding ', kind, ': ',
			paste(dQuote(setdiff(n, names(x))), collapse = ', '), '\n'
		),
		if (length(x) == m) paste(
			' * Since both spectra and', kind,
			'are named, refusing to match them despite having', m,
			'spectra in both datasets\n'
		)
	),
	call. = FALSE
)

feemife <- function(x, ...) UseMethod('feemife')

feemife.list <- function(x, absorbance, abs.path, ..., progress = FALSE) {
	if (missing(abs.path)) abs.path <- rep(1, length(x))
	stopifnot(length(list(...)) == 0)
	cubeapply(
		x, feemife,
		arrange(abs2list(absorbance), names(x), length(x), 'absorbance data'),
		arrange(abs.path, names(x), length(x), 'cell lengths'),
		progress = progress, .recycle = TRUE
	)
}

# NB: this works because cubeapply.feemcube(x, f) immediately calls
# f(as.list(x)), landing us in feemife.list above
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

