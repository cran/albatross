# "cubeapply" works like lapply or Map on feemcubes or lists of feems.
# If an error occurs, it's augmented with the information about the
# sample that caused it, but the original traceback is preserved (at the
# cost of making it noticeably longer).
# Optionally, there's a progress bar and support for parallel
# processing.
cubeapply <- function(x, ...) UseMethod('cubeapply')
cubeapply.feemcube <- function(x, fun, ...) feemcube(fun(as.list(x), ...), TRUE)

# make function to be called from cubeapply wrapper/error handler
.wraperr <- function(name) function(e) {
	ee <- simpleError(paste0(
		"While processing ", dQuote(name), ": ", conditionMessage(e)
	))
	ee$parent <- e
	class(ee) <- c('feem.wrapped.error', class(ee))
	stop(ee)
}
# wrap user-provided function to unwrap the feem and call handler on error
.wrapfun <- function(fun)
	function(l, ...) withCallingHandlers(fun(l$x, ...), error = .wraperr(l$name))
# wrap user-provided list with names for .wrapfun to unwrap it
.wraplist <- function(x, nm = x)
	Map(function(x, n) list(x = x, name = n), x, nm)

cubeapply.list <- function(x, fun, ..., cl, progress = TRUE, .recycle = FALSE) {
	# prepare to handle errors inside the loop
	nm <- .cubenames(x)
	x <- .wraplist(x, nm)
	wfun <- .wrapfun(fun)
	# run the loop
	if (progress) {
		pb <- txtProgressBar(max = length(x), style = 3)
		on.exit(close(pb))
	}
	if (missing(cl)) { # sequential processing: use library functions
		pfun <- if (progress) function(...) {
			# increment the progress bar every time the function returns
			on.exit(setTxtProgressBar(pb, getTxtProgressBar(pb) + 1))
			wfun(...)
		} else wfun
		if (.recycle) mapply(
			pfun, x, ..., SIMPLIFY = FALSE
		) else lapply(x, pfun, ...)
	} else {
		if (progress) {
			# how much is already done
			i <- 0
			# in chunks of length n
			n <- length(clusterEvalQ(cl, 1))
			# store returned values here
			ret <- list()
			while (i < length(x)) {
				# chunk must not exceed total work size
				todo <- (i+1):min(length(x), i + n)
				ret[todo] <- if (.recycle) do.call(clusterMap, c(
					# .recycle = TRUE is like Map(): need to subset all
					# arguments, not only x
					list(cl, wfun, x[todo]), lapply(list(...), `[`, todo)
				)) else parLapply(cl, x[todo], wfun, ...)
				i <- max(todo)
				if (progress) setTxtProgressBar(pb, i)
			}
			setNames(ret, names(x))
		} else {
			if (.recycle) clusterMap(
				cl, wfun, x, ...
			) else parLapply(
				cl, x, wfun, ...
			)
		}
	}
}

.cubenames <- function(x) UseMethod('.cubenames')

.fixnames <- function(x, n) as.factor(if (is.null(x)) seq_len(n) else make.unique(x))

.cubenames.feemcube <- function(cube) .fixnames(dimnames(cube)[[3]], dim(cube)[3])
.cubenames.list <- function(l) .fixnames(names(l), length(l))
