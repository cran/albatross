vtolProgressBar <- function(ctol) {
	i <- 0
	maxvtol <- NULL
	up <- function(relSSE, vtol) {
		i <<- i + 1
		if (is.null(maxvtol) || vtol > maxvtol) maxvtol <<- vtol

		msg <- sprintf('#%d SSE/SSX=%e rel.dSSE=%e', i, relSSE, vtol)
		if ((msg.width <- nchar(msg, 'width') + 1) < getOption('width')) {
			pbval <- min(1, max(0, log(vtol/maxvtol) / log(ctol / maxvtol)))
			bar.width <- getOption('width') - msg.width
			bar.fill <- trunc(bar.width * pbval)
			msg <- c(
				strrep(c('=', ' '), c(bar.fill, bar.width - bar.fill)),
				' ', msg
			)
		}
		cat('\r', msg, sep = '')
	}
	close <- function() cat('\n')
	list(up = up, close = close)
}
