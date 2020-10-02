.seq <- function(range, len, gamma) if (length(range) == 1) {
	rep_len(range, len)
} else {
	range[1] + diff(range) * seq(0, 1, length.out = len) ^ gamma
}

marine.colours <- function(
	n, chroma = .65, luminance = c(.35, 1),
	alpha = 1, gamma = 1, fixup = TRUE
) hcl(
	h = .seq(c(280, 60), n, gamma), # blue-violet to blue, green, yellow
	c = .seq(100 * chroma, n, gamma), l = .seq(100 * luminance, n, gamma),
	alpha = .seq(alpha, n, gamma), fixup = fixup
)
