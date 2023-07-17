(function(e) {
	Raman.shift <- 3400
	feems.param <- albatross:::feems.param
	nfac <- dim(feems.param[[1]])[1]

	# always generate same scores, subject to RNG kind
	if (exists('.Random.seed', .GlobalEnv)) oldseed <- .GlobalEnv$.Random.seed
	set.seed(42)
	on.exit({
		if (exists('oldseed')) .GlobalEnv$Random.seed <- oldseed
		else rm('.Random.seed', envir = .GlobalEnv)
	})
	C <- matrix(runif(12 * nfac), ncol = nfac)

	make.loads <- function(l.em, l.ex) c(
		Map(function(p, x) {
			apply(p, 1, function(p) { # for every factor...
				rowSums(apply(p, 1, function(p) { # sum both peaks...
					p['a'] * exp(-(x - p['m'])^2/(2 * p['s']^2))
				}))
			})
		}, feems.param, list(l.em, l.ex)),
		list(em = l.em, ex = l.ex)
	)

	make.cube <- function(l) {
		dimX <- c(nrow(l$A), nrow(l$B), nrow(C))
		X <- array( # fluorescence
			tcrossprod(l$A, multiway::krprod(C, l$B)), dimX
		)
		X <- X + multiway::nscale( # noise
			array(rnorm(length(X)), dimX),
			0, ssnew = multiway::sumsq(X) * .05
		)
		X <- X + multiway::nscale(
			array(outer(l$em, l$ex, function(em, ex) ex^-4 * (
				dnorm(em, ex, 10) +
				.2 * dnorm(em, 1/(1/ex - Raman.shift/1e7), 9)
			)), dimX), 0, multiway::sumsq(X) * 5
		)
		albatross::feemcube(X, l$em, l$ex, names = letters[1:nrow(C)])
	}

	# a small cube, almost ready for PARAFAC (except for scattering)
	l.em <- seq(240, 435, length.out = 32)
	l.ex <- seq(230, 350, length.out = 10)
	cube <- make.cube(make.loads(l.em, l.ex))

	## a larger dataset with IFE and different sizes
	l.em <- seq(240, 450, 2)
	l.ex <- seq(230, 350, 5)

	## primitive model, assuming that all absorption stems from fluorescence only
	feems <- make.cube(make.loads(l.em, l.ex))
	l.abs <- seq(min(l.em, l.ex), max(l.em, l.ex), 1)
	absorp <- setNames(apply(
		with(make.loads(l.em, l.abs), tcrossprod(B, C)/10),
		2, function(x) cbind(l.abs, x), simplify = FALSE
	), dimnames(feems)[[3]])
	feems <- feems * (feems / albatross::feemife(feems, absorp))
	feems <- as.list(feems)

	e$feems <- feems
	e$cube <- cube
	e$absorp <- absorp
})(environment())
