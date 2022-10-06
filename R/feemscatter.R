# Here, we locate the scattering region and wrap interpolation methods
# to present them as a common interface. We implement Whittaker
# smoothing ourselves in a separate file.

feemscatter <- function(x, ...) UseMethod('feemscatter')

feemscatter.list <- feemscatter.feemcube <- function(x, ..., cl, progress = TRUE)
	cubeapply(x, feemscatter, ..., cl = cl, progress = progress)

# pchip requires sorted x0 and at least 3 defined points
.pchip <- function(x0, y0, xi) if (length(x0) >= 3) {
	pchip(sort(x0), y0[order(x0)], xi)
} else { # fall back to linear interpolation
	interp1(sort(x0), y0[order(x0)], xi)
}

do.interpolate.pchip <- function(feem, mask, l.row, l.col) {
	# interpolate from everything defined not in target set
	src <- !mask & is.finite(feem)

	minmax.row <- c(which.min(l.row), which.max(l.row))
	minmax.col <- c(which.min(l.col), which.max(l.col))

	# need to provide a pair of fully defined columns on border
	# wavelengths (unless they are already defined) to avoid
	# extrapolation

	# 1. provide zeroes as a last resort
	feem[minmax.row, minmax.col] <- ifelse(
		src[minmax.row, minmax.col],
		feem[minmax.row, minmax.col], 0
	)
	src[minmax.row, minmax.col] <- TRUE

	# 2. pchip-interpolate anything missing
	for (i in 1:2) { # min and max
		feem[minmax.row[i], !src[minmax.row[i],]] <- .pchip(
			l.col[src[minmax.row[i],]], feem[minmax.row[i], src[minmax.row[i],]],
			l.col[!src[minmax.row[i],]]
		)
		src[minmax.row[i],] <- TRUE
	}

	# 3. ready to work column by column
	for (i.col in seq_along(l.col)) {
		.pchip(
			l.row[src[,i.col]], feem[src[,i.col], i.col], l.row[mask[,i.col]]
		) -> feem[mask[,i.col], i.col]
	}

	feem
}

interpolate.pchip <- function(
	feem, mask, by = c('emission', 'excitation', 'both')
) switch(match.arg(by),
	emission = do.interpolate.pchip(
		feem, mask, attr(feem, 'emission'), attr(feem, 'excitation')
	),
	excitation = t(do.interpolate.pchip(
		# easier to reuse existing method than to reimplement
		# it in terms of the other axis
		t(feem), t(mask), attr(feem, 'excitation'), attr(feem, 'emission')
	)),
	both = (
		Recall(feem, mask, 'emission') +
		Recall(feem, mask, 'excitation')
	)/2
)

interpolate.loess <- function(feem, mask, span = .05, ...) {
	l.em <- attr(feem, 'emission')
	l.ex <- attr(feem, 'excitation')

	# loess understands three-column format
	xx <- l.em[row(feem)][!mask]
	yy <- l.ex[col(feem)][!mask]
	zz <- feem[!mask]

	x0 <- l.em[row(feem)][mask]
	y0 <- l.ex[col(feem)][mask]
	z0.hat <- predict(
		loess(
			zz ~ xx + yy, data.frame(xx = xx, yy = yy, zz = zz),
			span = span, ...
		),
		data.frame(xx = x0, yy = y0)
	)
	z0.hat[z0.hat < 0] <- 0 # LOESS does sometimes return negative values

	feem[mask] <- z0.hat

	feem
}

interpolate.kriging <- function(feem, mask, ...) {
	l.em <- attr(feem, 'emission')
	l.ex <- attr(feem, 'excitation')

	src <- !mask & is.finite(feem)
	xx <- l.em[row(feem)][src]
	yy <- l.ex[col(feem)][src]
	zz <- feem[src]

	x0 <- l.em[row(feem)][mask]
	y0 <- l.ex[col(feem)][mask]
	z0.hat <- kriging(cbind(xx, yy), zz, cbind(x0, y0), ...)
	z0.hat[z0.hat < 0] <- 0 # last resort against negative values

	feem[mask] <- z0.hat

	feem
}

interpolate.whittaker <- function(
	feem, mask, d = 1:2, lambda = c(1e-2, 10), logscale = NA, nonneg = 1
) {
	feem[mask] <- NA # request the points to be interpolated
	feem[mask] <- whittaker2(
		attr(feem, 'emission'), attr(feem, 'excitation'),
		feem, lambda = lambda, d = d, logscale = logscale, nonneg = nonneg
	)[mask] # substitute only requested zone
	feem
}

.feeminterpolate <- function(x, method, mask, ...) {
	# some interpolation methods don't handle empty mask well
	if (all(!mask)) return(x)
	switch(method,
		pchip     = interpolate.pchip,
		loess     = interpolate.loess,
		kriging   = interpolate.kriging,
		whittaker = interpolate.whittaker
	)(x, mask, ...)
}

.scatter.mask <- function(x, widths, Raman.shift) {
	# interpret single wavelength widths as +/- width
	widths <- lapply(widths, function(p) abs(if (length(p) == 1) rep(p, 2) else p))
	# group by Rayleigh, Raman
	widths <- split(widths, ceiling(seq_along(widths)/2))
	outer(
		attr(x, 'emission'), attr(x, 'excitation'),
		function(em, ex) Reduce(`|`, Map(
			function(w, k) {
				Ray <- em/k - ex
				Ram <- em/k - 1/(1/ex - Raman.shift/1e7)
				(Ray > -w[[1]][1] & Ray < w[[1]][2]) | (
					if (length(w) >= 2)
						(Ram > -w[[2]][1] & Ram < w[[2]][2])
					else FALSE
				)
			}, widths, as.numeric(names(widths))
		))
	)
}

feemscatter.feem <- function(
	x, widths, method = c('omit', 'pchip', 'loess', 'kriging', 'whittaker'),
	add.zeroes = 30, Raman.shift = 3400, ...
) {
	scatter <- .scatter.mask(x, widths, Raman.shift)
	if (!is.na(add.zeroes)) x[
		is.na(x) & outer(
			attr(x, 'emission'), attr(x, 'excitation') - add.zeroes,
			`<`
		)
	] <- 0
	switch(method <- match.arg(method),
		omit = { x[scatter] <- NA; x },
		.feeminterpolate(x, method, scatter, ...)
	)
}
