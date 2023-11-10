# "bootstrap" computes parafac decompositions with different subsets and
# parameters, matches and rescales them to be comparable more easily,
# optionally performs some post-processing (see feemjackknife).

match.factors <- function(target, current) {
	# 1. reorder the maximally matching components together
	current <- reorder(current, like = target)

	# 2. rescale the matching components
	# by minimizing L2 norm of reconstruction error per each mode
	current <- rescale(
		current, mode = c('A', 'B'), absorb = 'C', like = target
	)

	current
}

# slices should be a list of indices in cube
# args[[i]], ..., verbose = FALSE are passed to feemparafac()
# postprocess() is called on the result of that
# feemparafac results are reordered&rescaled to fit the first result with
# the same number of components
# returns the list of overall results
bootparafac <- function(
	cube, slices, ..., args = vector('list', length(slices)), postprocess,
	progress = TRUE
) {
	if (progress) {
		pb <- txtProgressBar(
			max = length(slices), style = if (interactive()) 3 else 1
		)
		on.exit(close(pb))
	}

	ret <- list()
	# keep track of nfac because it may differ between runs
	ncomps <- integer()

	# store X in an environment to prevent it from being copied
	env <- new.env(parent = emptyenv())
	env$X <- cube

	for (i in seq_along(slices)) {
		# calculate the model
		ret[[i]] <- do.call(
			feemparafac,
			c(
				list(
					X = 'X', ..., verbose = FALSE,
					subset = slices[[i]], envir = env
				),
				args[[i]]
			)
		)
		# unless we are the first sample to have this nfac, match this run
		# to the first run with the same nfac
		ncomps[i] <- ncol(ret[[i]]$A)
		if (any(ncomps[i] == ncomps[-i])) {
			ret[[i]] <- match.factors(
				ret[[which(ncomps[i] == ncomps[-i])[1]]], ret[[i]]
			)
		}
		# postprocess, if any requested
		if (!missing(postprocess)) ret[[i]] <- postprocess(
			ret[[i]], cube, slices[[i]], args[[i]], ...
		)
		# finally done with this run
		if (progress) setTxtProgressBar(pb, i)
	}

	ret
}
