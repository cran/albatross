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
.wrapfun <- function(fun) {
	force(fun)
	# NOTE: without the force(), the resulting object could serialize to 700
	# megabytes due to `fun` seemingly being a promise. Once the promise
	# resolves, suddenly `fun` is much more compact.
	function(l, ...) withCallingHandlers(fun(l$x, ...), error = .wraperr(l$name))
}
# wrap user-provided list with names for .wrapfun to unwrap it
.wraplist <- function(x, nm = x)
	Map(function(x, n) list(x = x, name = n), x, nm)

cubeapply.list <- function(
	x, fun, ..., cl, progress = TRUE, .recycle = FALSE,
	chunk.size = NULL, .scheduling = c('static', 'dynamic')
) {
	# prepare to handle errors inside the loop
	nm <- .cubenames(x)
	x <- .wraplist(x, nm)
	wfun <- .wrapfun(fun)
	# run the loop
	if (progress) {
		pb <- makepb(length(x))
		on.exit(close(pb))
	}
	if (missing(cl)) { # sequential processing: use library functions
		pfun <- if (progress) wrap_pb_callback(wfun, pb) else wfun
		if (.recycle) mapply(
			pfun, x, ..., SIMPLIFY = FALSE
		) else lapply(x, pfun, ...)
	} else {
		# can't get default cluster on R < 3.5.0, but need length(cl)
		cl <- maybe_default_cluster(cl)
		if (is.null(cl)) progress <- FALSE
		.scheduling <- match.arg(.scheduling)
		if (progress) {
			nchunks <- if (.recycle) length(x) else {
				# logic taken from parallel:::staticNChunks, dynamicNChunks
				if (is.null(chunk.size))
					length(cl) * (1 + (.scheduling == 'dynamic'))
				else if (chunk.size <= 0)
					c(length(cl), length(x))[1 + (.scheduling == 'dynamic')]
				else
					max(1, ceiling(length(x) / chunk.size))
			}
			pb <- makepb(nchunks)
			on.exit(close(pb), add = TRUE)
			cl <- wrap_cluster(cl, make_pb_callback(pb))
			on.exit(unwrap_cluster(cl), add = TRUE)
		}
		if (.recycle) clusterMap(
			cl, wfun, x, ..., .scheduling = .scheduling
		) else switch(.scheduling,
			static = parLapply,
			dynamic = parLapplyLB
		)(cl, x, wfun, ..., chunk.size = chunk.size)
	}
}

.cubenames <- function(x) UseMethod('.cubenames')

.fixnames <- function(x, n) as.factor(if (is.null(x)) seq_len(n) else make.unique(x))

.cubenames.feemcube <- function(x) .fixnames(dimnames(x)[[3]], dim(x)[3])
.cubenames.list     <- function(x) .fixnames(names(x),         length(x))
