# Use interpolation to change the wavelength grid of a FEEM or a group
# of them. Not very reliable, can be easily made to return nonsense if
# there's large areas of NAs.

feemgrid <- function(x, ...) UseMethod('feemgrid')

.cliprange <- function(wl, ranges)
	wl[wl >= max(ranges[1,]) & wl <= min(ranges[2,])]

feemgrid.feemcube <- feemgrid.list <- function(
	x, emission, excitation, ..., progress = TRUE
) {
	if (missing(emission) || missing(excitation)) {
		# need to find out the overall grid
		overall <- if (inherits(x, 'feemcube'))
			attributes(x)[c('emission', 'excitation')]
		else # unfortunately we copy the code from feemcube()
			Map(sort, Reduce(
				# union(em.a, em.b), union(ex.a, ex.b)
				function(a, b) Map(union, a, b),
				Map(
					function(l) attributes(l)[c('emission', 'excitation')],
					x
				)
			))

		# need the individual wavelength ranges to choose the strictest
		wlranges <- vapply(
			as.list(x),
			function(x) {
				x <- as.data.frame(x) # subsets !is.na(x)
				vapply(
					setNames(nm = c('emission', 'excitation')),
					function(n) range(x[,n]), numeric(2)
				)
			},
			array(NA_real_, dim = c(2, 2))
		)

		# now clip the overall range to the subranges
		if (missing(emission)) emission <- .cliprange(
			overall$emission, wlranges[,'emission',]
		)
		if (missing(excitation)) excitation <- .cliprange(
			overall$excitation, wlranges[,'excitation',]
		)
	}
	cubeapply(
		x, feemgrid, emission = emission, excitation = excitation,
		..., progress = progress
	)
}

feemgrid.feem <- function(
	x, emission, excitation,
	method = c('whittaker', 'loess', 'kriging', 'pchip'), ...
) {
	# temporarily construct a union of old and new scales to interpolate
	temp.x <- feem(
		matrix(
			NA_real_,
			length(temp.em <- union(attr(x, 'emission'), emission)),
			length(temp.ex <- union(attr(x, 'excitation'), excitation))
		),
		temp.em, temp.ex, attr(x, 'scale')
	)
	# assign all known data
	temp.x[
		match(attr(x, 'emission'), temp.em),
		match(attr(x, 'excitation'), temp.ex)
	] <- x[]
	# interpolate anything unknown
	temp.x <- .feeminterpolate(temp.x, match.arg(method), is.na(temp.x), ...)
	# only return the requested grid
	temp.x[match(emission, temp.em), match(excitation, temp.ex), drop = FALSE]
}
