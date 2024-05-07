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

.bootparafac_taskfun <- function(slice, args, X, ...) {
	ret <- do.call(
		albatross::feemparafac,
		c(
			list(X = X, subset = slice, verbose = FALSE, ...),
			args
		)
	)
	# NOTE: must repair the object once it's received
	structure(ret, cube = NULL)
}

# slices should be a list of indices in cube
# args[[i]], ..., verbose = FALSE are passed to feemparafac()
# postprocess() is called on the result of that
# feemparafac results are reordered&rescaled to fit the first result with
# the same number of components
# returns the list of overall results
bootparafac <- function(
	cube, slices, parallel = FALSE, cl = NULL, ...,
	args = vector('list', length(slices)), postprocess,
	progress = TRUE, .scheduling = 'dynamic'
) {
	fun <- if (progress) {
		pb <- makepb(length(slices))
		on.exit(close(pb), add = TRUE)
		if (parallel) {
			cl <- wrap_cluster(cl, make_pb_callback(pb))
			on.exit(unwrap_cluster(cl), add = TRUE)
			.bootparafac_taskfun
		} else {
			wrap_pb_callback(.bootparafac_taskfun, pb)
		}
	} else .bootparafac_taskfun
	# Let the individual models run concurrently, without the jobs in
	# flight to `nstart` (defaults to 10) and having to wait for the
	# last random start to finish before continuing.
	ret <- if (parallel) clusterMap(
		cl, fun,
		slice = slices, args = args,
		MoreArgs = list(X = cube, ...), RECYCLE = FALSE,
		SIMPLIFY = FALSE, USE.NAMES = FALSE, .scheduling = .scheduling
	) else Map(
		fun, slice = slices, args = args,
		MoreArgs = list(X = cube, ...),
		USE.NAMES = FALSE
	)

	# Reattach the cube back, referenced through an environment
	env <- new.env(parent = emptyenv())
	env$X <- cube
	ret <- lapply(ret, structure, envir = env, cube = 'X')

	# Rearrange models as needed. This is supposedly fast and much
	# easier to perform locally.
	ncomps <- integer()
	for (i in seq_along(ret)) {
		# Unconstrained models may have arbitrary signs. Since we match
		# by maximum cosine similarity of the loadings, make sure that
		# the loadings are maximally positive.
		ret[[i]] <- resign(resign(ret[[i]], 'A', 1, 'C'), 'B', 1, 'C')
		# unless we are the first sample to have this nfac, match this run
		# to the first run with the same nfac
		ncomps[i] <- ncol(ret[[i]]$A)
		if (any(ncomps[i] == ncomps[-i])) {
			ret[[i]] <- match.factors(
				ret[[which(ncomps[i] == ncomps[-i])[1]]], ret[[i]]
			)
		}
	}

	# postprocessing results in large tasks, better done locally
	if (!missing(postprocess)) ret <- Map(
		function(model, cube, slice, args, ...)
			postprocess(model, cube, slice, args, ...),
		model = ret, slice = slices, args = args,
		MoreArgs = list(cube = cube, ...)
	)

	ret
}
