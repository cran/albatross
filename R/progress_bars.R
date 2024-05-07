# Currently, this progress bar is only used in feemflame(), and the
# interface is very rough.
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
		flush.console()
	}
	close <- function() cat('\n')
	list(up = up, close = close)
}

makepb <- function(len)
	txtProgressBar(max = len, style = if (interactive()) 3 else 1)

make_pb_callback <- function(pb) {
	force(pb)
	function() setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
}

wrap_pb_callback <- function(fun, pb) {
	force(fun)
	force(pb)
	function(...) {
		on.exit(setTxtProgressBar(pb, getTxtProgressBar(pb) + 1))
		fun(...)
	}
}
