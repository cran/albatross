feemjackknife <- function(cube, ..., progress = TRUE) {
	slices <- c(list(1:dim(cube)[3]), lapply(1:dim(cube)[3], `-`))
	models <- bootparafac(
		cube, slices, ...,
		postprocess = function(fac, cube, slice, args, ...) {
			if (length((1:dim(cube)[3])[-slice]) == 1) {
				# X[i,j] = sum(A[i,] * B[j,] * C[]). C = ?
				# combine D[i,j,r] <- A[i,r] * B[j,r], unfold (i,j) -> l
				# X[l] = sum(D[l,] * C[])
				# this is a simple matrix equation:
				# x = Dc; D'x = D'Dc; c = (D'D)^(-1) D' x = D^+ x
				# we need to unfold cube[,,-slice] and A x B,
				# then multiply pseudoinverse of the latter by the former

				X <- cube[,,-slice]
				X <- X * attr(X, 'scale')
				dim(X) <- prod(dim(X))
				mask <- !is.na(X) # skip NA regions in the spectrum
				D <- krprod(fac$B, fac$A)
				# we could have used ginv() from recommended MASS package,
				# but since we already use pracma, let's use its pinv() and
				# save one dependency
				attr(fac, 'Chat') <- t(pinv(D[mask,]) %*% X[mask])
			}
			fac
		}, progress = progress
	)
	structure(
		list(overall = models[[1]], leaveone = models[-1]),
		class = 'feemjackknife'
	)
}

plot.feemjackknife <- function(x, kind = c('estimations', 'RIP', 'IMP'), ...)
	switch(
		match.arg(kind),
		estimations = jkplot(x, ...),
		RIP = jk.RIP(x, ...),
		IMP = jk.IMP(x, ...)
	)

jkplot <- function(
	jk, xlab = quote(lambda*', nm'), ylab = 'Loading values', as.table = T,
	scales = list(x = 'free'), ...
) {
	ovcube <- .pfcube(jk$overall)
	df <- do.call(rbind, lapply(seq_along(jk$leaveone), function(i) rbind(
		data.frame(
			loading = as.vector(jk$leaveone[[i]]$A),
			mode = 'Emission',
			wavelength = attr(ovcube, 'emission'),
			factor = as.factor(col(jk$leaveone[[i]]$A)),
			rep = i
		),
		data.frame(
			loading = as.vector(jk$leaveone[[i]]$B),
			mode = 'Excitation',
			wavelength = attr(ovcube, 'excitation'),
			factor = as.factor(col(jk$leaveone[[i]]$B)),
			rep = i
		)
	)))
	xyplot(
		loading ~ wavelength | mode + factor, df, group = rep,
		type = 'l', as.table = as.table, xlab = xlab, ylab = ylab,
		scales = scales, ...
	)
}

jk.RIP <- function(
	jk, q = .9, xlab = 'Mean squared residuals',
	ylab = 'Mean squared difference in loadings',
	scales = list(alternating = 1), ...
) {
	RIP <- do.call(rbind, lapply(jk$leaveone, function(fac) data.frame(
		msq.resid = mean(resid(fac)^2, na.rm = T),
		Emission = mean((fac$A - jk$overall$A)^2, na.rm = T),
		Excitation = mean((fac$B - jk$overall$B)^2, na.rm = T)
	)))
	cube <- .pfcube(jk$overall)
	if (!is.null(dimnames(cube)[[3]])) rownames(RIP) <- dimnames(cube)[[3]]

	xyplot(
		Emission + Excitation ~ msq.resid, RIP, outer = T,
		xlab = xlab, ylab = ylab, scales = scales,
		panel = function(x, y, ...) {
			panel.xyplot(x, y, ...)
			outl <- x > quantile(x, q) | y > quantile(y, q)
			ltext(x[outl], y[outl], rownames(RIP)[outl])
		},
		...
	)
}

jk.IMP <- function(
	jk, q = .9, xlab = 'Overall model scores',
	ylab = 'Individual model scores', as.table = T,
	scales = list(alternating = 1), ...
) {
	Chat <- do.call(rbind, lapply(jk$leaveone, attr, 'Chat'))
	IMP <- data.frame(
		score.overall = as.vector(jk$overall$C),
		score.predicted = as.vector(Chat),
		factor = as.factor(col(Chat))
	)
	xyplot(
		score.predicted ~ score.overall | factor, IMP,
		xlab = xlab, ylab = ylab, scales = scales, as.table = as.table,
		panel = function(x, y, ...) {
			panel.xyplot(x, y, ...)
			panel.abline(0, 1, lwd = .5)
			outl <- abs(y - x)
			outl <- outl > quantile(outl, q)
			ltext(x[outl], y[outl], rownames(IMP)[outl])
		},
		...
	)
}
