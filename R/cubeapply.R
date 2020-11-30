# lapply or Map function fun to each element of x and ...
# optional progress bar and parallel processing
cubeapply <- function(x, ...) UseMethod('cubeapply')
cubeapply.feemcube <- function(x, fun, ...) feemcube(fun(as.list(x), ...), TRUE)

cubeapply.list <- function(x, fun, ..., cl, progress = TRUE, .recycle = FALSE) {
	if (progress) {
		pb <- txtProgressBar(max = length(x), style = 3)
		on.exit(close(pb))
	}
	if (missing(cl)) { # sequential processing: use library functions
		pfun <- if (progress) function(...) {
			# increment the progress bar every time the function returns
			on.exit(setTxtProgressBar(pb, getTxtProgressBar(pb) + 1))
			fun(...)
		} else fun
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
					list(cl, fun, x[todo]), lapply(list(...), `[`, todo)
				)) else parLapply(cl, x[todo], fun, ...)
				i <- max(todo)
				if (progress) setTxtProgressBar(pb, i)
			}
			setNames(ret, names(x))
		} else {
			if (.recycle) clusterMap(
				cl, fun, x, ...
			) else parLapply(
				cl, x, fun, ...
			)
		}
	}
}
