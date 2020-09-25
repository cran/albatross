feemscatter <- function(x, ...) UseMethod('feemscatter')

feemscatter.list <- function(x, ...) lapply(x, feemscatter, ...)

feemscatter.feemcube <- function(x, ...) feemcube(feemscatter(as.list(x), ...), TRUE)

omit.mask <- function(feem, mask) {
	feem[mask] <- NA
	feem
}

# pchip requires sorted x0
.pchip <- function(x0, y0, xi) pchip(sort(x0), y0[order(x0)], xi)

interpolate.pchip <- function(feem, mask) {
	# interpolate from everything defined not in target set
	src <- !mask & is.finite(feem)

	l.em <- attr(feem, 'emission')
	l.ex <- attr(feem, 'excitation')
	minmax.em <- c(which.min(l.em), which.max(l.em))
	minmax.ex <- c(which.min(l.ex), which.max(l.ex))

	# 1. provide a row of zeroes at minimal emission wavelength if needed
	feem[minmax.em[1],][!src[minmax.em[1],]] <- 0
	src[minmax.em[1],] <- TRUE

	# provide a fully defined excitation spectrum at maximal wavelength,
	# unless it's already defined
	# 2. provide zeroes as a last resort, since pchip doesn't like NAs
	feem[minmax.em[2], minmax.ex] <- ifelse(
		src[minmax.em[2], minmax.ex],
		feem[minmax.em[2], minmax.ex], 0
	)
	src[minmax.em[2], minmax.ex] <- TRUE

	# 3. pchip-interpolate anything missing
	feem[minmax.em[2], !src[minmax.em[2],]] <- .pchip(
		l.ex[src[minmax.em[2],]], feem[minmax.em[2], src[minmax.em[2],]],
		l.ex[!src[minmax.em[2],]]
	)
	src[minmax.em[2],] <- TRUE

	for (i.ex in seq_along(l.ex)) {
		.pchip(
			l.em[src[,i.ex]], feem[src[,i.ex], i.ex], l.em[mask[,i.ex]]
		) -> feem[mask[,i.ex], i.ex]
	}
	feem
}

interpolate.loess <- function(feem, mask, span = .05, ...) {
	l.em <- attr(feem, 'emission')
	l.ex <- attr(feem, 'excitation')

	# loess understands three-column format
	xx <- l.em[row(feem)][!mask]
	yy <- l.ex[col(feem)][!mask]
	zz <- feem[!mask]

	x0 <- l.em[row(feem)][mask]
	y0 <- l.ex[col(feem)][mask]

	feem[mask] <- predict(
		loess(
			zz ~ xx + yy, data.frame(xx = xx, yy = yy, zz = zz),
			span = span, ...
		),
		data.frame(xx = x0, yy = y0)
	)

	feem[feem < 0] <- 0 # LOESS does sometimes return negative values

	feem
}

feemscatter.feem <- function(
	x, widths, method = c('omit', 'pchip', 'loess'),
	add.zeroes = 30, Raman.shift = 3400, ...
) {
	scatter <- outer(
		attr(x, 'emission'), attr(x, 'excitation'),
		function(em, ex)
			# Rayleigh, Raman, 2*Rayleigh, 2*Raman
			abs(em - ex) < widths[1] |
			abs(em - 1/(1/ex - Raman.shift/1e7)) < widths[2] |
			abs(em/2 - ex) < widths[3] |
			abs(em/2 - 1/(1/ex - Raman.shift/1e7)) < widths[4]
	)
	if (!is.na(add.zeroes)) x[
		is.na(x) & outer(
			attr(x, 'emission'), attr(x, 'excitation') - add.zeroes,
			`<`
		)
	] <- 0
	switch(match.arg(method),
		omit = omit.mask,
		pchip = interpolate.pchip,
		loess = interpolate.loess
	)(x, scatter, ...)
}
