# Fluorescence indices are functions of FEEMs. The code that extracts
# the wavelengths in a mostly uniform way is somewhat complicated, but
# it prevents boring repetitive code later.

feemindex <- function(x, ...) UseMethod('feemindex')

# must *always* return a data.frame
feemindex.feemcube <- feemindex.list <- function(x, ..., progress = FALSE) {
	ret <- cubeapply(as.list(x), feemindex, ..., progress = progress)
	cbind(
		sample = .cubenames(x),
		as.data.frame(matrix(
			# NB: simplify2array doesn't do the right thing here
			unlist(ret), length(ret), byrow = TRUE,
			dimnames = list(NULL, names(ret[[1]]))
		))
	)
}

# extract the single wavelenghts or ranges (HIX, YFI, FrI, C) from x
# must handle following cases:
# - desired wavelength exists in attr(x, 'em' or 'ex')
# - desired wavelength exists within +/- tolerance[1 or 2] (should be 1 or
#   2 nm in both em and ex)
# - desired wavelength can be interpolated => should we warn? provide
#   an option to return NA instead
# - desired wavelength can't be interpolated (outside the range) =>
#   return NA or stop()?
# - desired wavelength exists but is NA => behave as if it's not in the
#   wavelength grid at all and handle it as above
# conditions:
# - diff(em and ex) must not be less than tolerance; otherwise same points could be
#   chosen for different wavelengths, leading to havoc
.extract <- function(x, em, ex, range, tolerance, method, ...) {
	kind <- c('emission', 'excitation')
	axes <- Map(function(need, avail, tolerance, kind) {
		if (range) # all wavelengths between requested
			need <- sort(unique(c(
				range(need),
				avail[avail >= min(need) & avail <= max(need)]
			)))
		distances <- abs(outer(need, avail, `-`))
		matches <- apply(distances, 1, which.min)
		distances <- distances[cbind(seq_along(need), matches)]
		# offer interpolation if allowed
		matches[distances > tolerance] <- if (is.character(method)) 0 else NA
		# but don't extrapolate
		matches[need < min(avail) | need > max(avail)] <- NA
		dup <- duplicated(matches, incomparables = c(0, NA))
		if (range) {
			# duplicates arise naturally when requested range boundary is
			# close to a grid point; can be safely dropped
			need <- need[!dup]
			matches <- matches[!dup]
		} else if (any(dup)) {
			# We matched same grid wavelength for multiple requested
			# wavelengths. This can't be useful.
			matches[dup] <- NA
			warning(
				'The same grid wavelength matches as closest to',
				' different wavelengths I need. Please decrease ', kind,
				' tolerance to at most ', min(diff(need)) / 2, ' and',
				' enable interpolation.'
			)
		}
		list(matches = matches, need = need)
	}, list(em, ex), attributes(x)[kind], tolerance, kind)
	if (any(unlist(lapply(axes, `[[`, 'matches')) == 0, na.rm = TRUE)) {
		# interpolation requested and needed
		new.em <- c(attr(x, 'emission'), subset(axes[[1]]$need, axes[[1]]$matches == 0))
		new.ex <- c(attr(x, 'excitation'), subset(axes[[2]]$need, axes[[2]]$matches == 0))
		# match again on spectrum with certain wls interpolated
		Recall(
			x = feemgrid(x, new.em, new.ex, method = method, ...),
			em, ex, range, tolerance, FALSE
		)
	} else x[axes[[1]]$matches, axes[[2]]$matches, drop = FALSE]
}

# integrate a piece of feem along the non-degenerate axis
.integ <- function(feem) {
	stopifnot(sum(dim(feem) != 1) == 1)
	x <- attributes(feem)[c('emission', 'excitation')][dim(feem) != 1][[1]]
	y <- as.vector(feem)
	if (any(is.na(x)) || any(is.na(y))) return(NA)
	# Trapezoidal method
	# \int_a^b f(x) dx \approx (b - a) \frac{f(a) + f(b)}{2}
	if (is.unsorted(x)) {
		perm <- order(x)
		x <- x[perm]
		y <- y[perm]
	}
	sum(diff(x) * (y[-1] + y[-length(y)]) / 2)
}

.ave <- function(feem) {
	stopifnot(sum(dim(feem) != 1) == 1)
	x <- attributes(feem)[c('emission', 'excitation')][dim(feem) != 1][[1]]
	.integ(feem) / diff(range(x))
}

.div <- function(x) {
	stopifnot(length(x) == 2)
	x[1]/x[2]
}

# must return named vector of requested fluorescence index values
feemindex.feem <- function(
	x, indices = c(
		# indices
		'HIX', 'BIX', 'MFI', 'CFI', 'YFI', 'FrI',
		# named peaks
		'A', 'B', 'C', 'M', 'P', 'T'
	),
	tolerance = 1, interpolate = FALSE, ...
) {
	# if interpolation is enabled, get rid of NAs
	if (is.character(interpolate))
		x <- .feeminterpolate(x, interpolate, is.na(x), ...)
	else
		stopifnot(length(list(...)) == 0)

	extract <- function(em, ex, range = FALSE) .extract(
		em = em, ex = ex, range = range,
		x = x, tolerance = tolerance, method = interpolate, ...
	)

	vapply(
		match.arg(indices, several.ok = TRUE),
		function(kind) switch(kind,
			HIX = .integ(extract(c(435, 480), 254, TRUE)) /
				.integ(extract(c(300, 345), 254, TRUE)),
			BIX = .div(extract(c(380, 430), 310)),
			MFI = .div(extract(c(450, 500), 370)),
			CFI = .div(extract(c(470, 520), 370)),
			YFI = .ave(extract(c(350, 400), 280, TRUE)) /
				.ave(extract(c(400, 450), 280, TRUE)),
			FrI = extract(380, 310) / max(extract(c(420, 435), 310, TRUE)),
			A = max(extract(c(400, 460), 260, TRUE)),
			B = extract(305, 275),
			C = max(extract(c(420, 460), c(320, 360), TRUE)),
			M = max(extract(c(370, 410), c(290, 310), TRUE)),
			P = extract(660, 398),
			T = extract(340, 275),
			stop("Don't know how to calculate ", dQuote(kind))
		),
		numeric(1)
	)
}
