# Functions typically computed on absorbance data alone (like
# fluorescence indices on fluorescence data) live here.
# See also: feemife.R for absorbance-related helpers.

# Return normalised RMSE for a model. Works for both lm and nls models.
.nrmse <- function(m) if (is.null(m)) NA else
	sqrt(mean(resid(m)^2)) / diff(range(fitted(m) - resid(m)))

# Fit y ~ exp(-slope * x) subject to a certain range of x and return
# the slope value and NRMSE (see above). The "log-linear" method is the
# one everybody uses (despite it lies when y gets close to 0) and the
# "exp-offset" method is the "more correct" one, taking a constant
# offset into account. The "hyperbolic" method uses a different model (a
# power function), also with a constant offset.
.spectral.slope <- function(
	x, y, range, method
) {
	mask <- x >= range[1] & x <= range[2]
	# NB: nonsensical situations include whole range falling between
	# grid points *and* asking for more than available
	if (!any(mask) || range[1] < min(x) || range[2] > max(x))
		return(list(slope = NA, model = NULL))
	x <- x[mask]
	y <- y[mask]
	tryCatch({
		slope <- switch(method,
			'log-linear' = {
				m <- lm(log(y) ~ x)
				-unname(coef(m)['x'])
			}
		)
		list(slope = slope, model = m)
	}, error = function(e) {
		warning(simpleWarning(conditionMessage(e), conditionCall(e)))
		list(slope = NA, model = NULL)
	})
}

.do.absindex <- function(
	X, len, unit, out.A, out.a, out.a.ratio, out.slope,
	out.slope.ratio, out.slope.nrmse, slope.method
) {
	wl <- X[,1]
	val <- X[,2]
	switch(unit,
		'log10' = {
			A <- val
			a <- log(10) * A / (len / 100)
		},
		'm^-1' = {
			a <- val
			A <- a / log(10) * (len / 100)
		}
	)
	# NB: sapply() won't work for two reasons. For slope ratios, it
	# would combine the vectors into matrices and lose their names.  For
	# everything else, it would return a list() when given a length-0
	# argument. Instead, we combine the lists as we get them, then make
	# the result into a numeric vector.
	unlist(c(
		if (length(out.A) > 0) with(
			spline(wl, A, xout = replace(
				out.A, out.A < min(wl) | out.A > max(wl), NA
			)),
			setNames(y, paste0('A.', x))
		),
		if (length(out.a) > 0) with(
			spline(wl, a, xout = replace(
				out.a, out.a < min(wl) | out.a > max(wl), NA
			)),
			setNames(y, paste0('a.', x))
		),
		lapply(out.a.ratio, function(r) {
			ar <- spline(
				wl, a, xout = replace(r, r < min(wl) | r > max(wl), NA)
			)$y
			setNames(ar[1] / ar[2], paste(c('aR', r), collapse = '.'))
		}),
		lapply(out.slope, function(r) {
			ans <- .spectral.slope(wl, a, r, slope.method)
			setNames(
				c(ans$slope, if (out.slope.nrmse) .nrmse(ans$model)),
				paste0(
					c('', if (out.slope.nrmse) 'NRMSE.'),
					paste(c('S', r), collapse = '.')
				)
			)
		}),
		lapply(out.slope.ratio, function(r) setNames(
			.spectral.slope(wl, a, r[1:2], slope.method)$slope /
				.spectral.slope(wl, a, r[3:4], slope.method)$slope,
			paste(c('SR', r), collapse = '.')
		))
	))
}

absindex <- function(
	x, abs.path, unit = c('log10', 'm^-1'),
	out.A = 254,
	out.a = c(350, 355, 374, 443),
	out.a.ratio = list(c(250, 365), c(465, 665)),
	out.slope = list(c(275, 295), c(350, 400)),
	out.slope.ratio = list(c(275, 295, 350, 400)),
	out.slope.nrmse = FALSE
) {
	x <- abs2list(x)
	if (missing(abs.path)) abs.path <- rep(1, length(x))
	else if (length(abs.path) == 1) abs.path <- rep(abs.path, length(x))
	stopifnot(
		lengths(out.a.ratio) == 2, lengths(out.slope) == 2,
		lengths(out.slope.ratio) == 4
	)
	ret <- do.call(rbind, .mapply(
		.do.absindex,
		list(
			X = x,
			len = arrange(
				abs.path, names(x), length(x), 'cell lengths'
			)
		),
		MoreArgs = list(
			unit = match.arg(unit), out.A = out.A, out.a = out.a,
			out.a.ratio = out.a.ratio, out.slope = out.slope,
			out.slope.ratio = out.slope.ratio,
			out.slope.nrmse = out.slope.nrmse,
			slope.method = 'log-linear'
		)
	))
	# workaround in case nothing is requested
	if (length(ret) == 0) ret <- matrix(numeric(), length(x), 0)
	cbind(sample = .fixnames(names(x), length(x)), as.data.frame(ret))
}
