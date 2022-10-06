# The cmocean('haline') palette is somewhat too dark for my tastes, so
# here are some perceptually-uniform (HCL) palette that I like better.

# x from [0,1] to [range]
.scl <- function(range, x) if (length(range) == 1) {
	rep_len(range, length(x))
} else {
	range[1] + diff(range) * x
}

marine.colours <- function(
	n, chroma = .65, luminance = c(.35, 1),
	alpha = 1, gamma = 1, fixup = TRUE
) {
	base <- seq(0, 1, length.out = n) ^ gamma
	hcl(
		h = .scl(c(280, 60), base), # blue-violet to blue, green, yellow
		c = .scl(100 * chroma, base), l = .scl(100 * luminance, base),
		alpha = .scl(alpha, base), fixup = fixup
	)
}

diverging.colours <- function(
	n, chroma = c(.35, .85), luminance = c(.30, .95),
	alpha = 1, gamma = 1, fixup = TRUE
) {
	base <- seq(-1, 1, length.out = n)
	sgn <- sign(base)
	val <- abs(base) ^ gamma

	hcl(
		# bright-blue via violet, reds, yellow to bright-green
		h = 150 * sgn * val,
		# maxima at val == .5, minima at val %in% c(0, 1)
		c = .scl(chroma * 100, (1 - abs(1 - 2 * val))),
		# maximum at endpoints, minimum exactly in the middle
		l = .scl(luminance * 100, val),
		alpha = .scl(alpha, val),
		fixup = fixup
	)
}
