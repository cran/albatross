.seq <- function(range, len) seq(
	range[1],
	ifelse(is.na(range[2]), range[1], range[2]),
	length.out = len
)

marine.colours <- function(
	n, chroma = .65, luminance = c(.35, 1),
	alpha = 1, gamma = 1, fixup = TRUE
) hcl(
	h = .seq(c(280, 60), n), # bluish-violet to yellow, via blue and green
	c = 100 * .seq(chroma, n),
	l = 100 * (.seq(luminance, n) ^ gamma),
	alpha = .seq(alpha, n), fixup = fixup
)
