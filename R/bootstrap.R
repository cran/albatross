# "bootstrap" computes parafac decompositions with different subsets and
# parameters, matches and rescales them to be comparable more easily,
# optionally performs some post-processing (see feemjackknife).

# argmin over c[] ||A - diag(c) %*% B||^2
.matscale <- function(A, B) {
	stopifnot(ncol(A) == ncol(B))
	vapply(
	# for every column:
	#  min over c[j] ||A[,j] - c[j] * B[,j]||
	#  c[j] = B[,j]^T A[,j] / B[,j]^T B[,j]
		1:ncol(A), function(j)
			crossprod(B[,j,drop=FALSE], A[,j,drop=FALSE]) /
				crossprod(B[,j,drop=FALSE]),
		numeric(1)
	)
}

match.factors <- function(target, current) {
	nfac <- ncol(target$A)

	# Tucker's congruence coefficient for emission & excitation components
	tcc <- as.matrix(pmin(congru(target$A, current$A), congru(target$B, current$B)))

	# 1. reorder maximally matching components together
	perm <- integer(nfac)
	for (i in seq_along(perm)) {
		# more than one component could fit perfectly
		# (despite it usually shouldn't), so choose the first match
		next.match <- which(tcc == max(tcc), arr.ind=T)[1,]
		# record the match and the distance
		perm[next.match[1]] <- next.match[2]
		# make sure that this pair won't match with anything else
		tcc[next.match[1],] <- -Inf
		tcc[,next.match[2]] <- -Inf
	}
	# NB: multiway <= 1.0-6 doesn't reorder 1-component models correctly
	if (nfac > 1) current <- reorder(current, perm)

	# 2. rescale the matching components
	# by minimizing L2 norm of reconstruction error per each mode
	scaling <- cbind(
		.matscale(target$A, current$A), .matscale(target$B, current$B)
	)
	# scaling matrix such that
	# Acurrent[r] * scaling[r,1] ~ Atarget
	# Bcurrent[r] * scaling[r,2] ~ Btarget
	# Ccurrent[r] * scaling[r,3] ~ Ctarget
	# and product over rows == 1
	scaling <- cbind(scaling, 1 / apply(scaling, 1, prod))

	current$A <- t(t(current$A) * scaling[,1])
	current$B <- t(t(current$B) * scaling[,2])
	current$C <- t(t(current$C) * scaling[,3])

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
		pb <- txtProgressBar(max = length(slices), style = 3)
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
