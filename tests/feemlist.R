library(albatross)

# this should always return paths under tempdir() or at least loop infinitely
mktempd <- function() repeat {
	dir <- tempfile('test_feemlist', tempdir())
	if (dir.create(dir, showWarnings = FALSE)) return(dir)
}

# on Windows we might get either path separator, normalise to hard-coded /
pfix <- if (.Platform$OS.type == 'windows') function(x) {
	gsub('\\', '/', x, fixed = TRUE)
} else identity

mktree <- function(root, tree) for (n in names(tree)) {
	if (is.list(tree[[n]])) {
		stopifnot(dir.create(file.path(root, n)))
		Recall(file.path(root, n), tree[[n]])
	} else {
		stopifnot(file.copy(tree[[n]], file.path(root, n)))
	}
}

assertnames <- function(actual, desired.names)
	stopifnot(setequal(
		names(actual),
		# to allow hardcoding desired.names with '/' in them
		gsub('/', .Platform$file.sep, desired.names, fixed = TRUE)
	))

# format support already tested in import.R; the important part here is
# filesystem walking
f <- system.file(file.path('extdata', 'panorama.dat'), package = 'albatross')

mktree(
	dir <- mktempd(),
	list('a.dat' = f, 'sub' = list('b.dat' = f, 'subsub' = list('c.dat' = f)))
)

# should simplify common file path components by default
assertnames(
	feemlist(dir, 'panorama'),
	c('a.dat', 'sub/b.dat', 'sub/subsub/c.dat')
)

# should use full (well, as full as given) file paths when asked to
assertnames(
	feemlist(dir, 'panorama', simplify.names = FALSE),
	file.path(dir, c('a.dat', 'sub/b.dat', 'sub/subsub/c.dat'))
)

mktree(
	dir2 <- mktempd(),
	list('d.dat' = f, 'sub' = list('e.dat' = f, 'subsub' = list('f.dat' = f)))
)

# less common paths if given 2 different directories
assertnames(
	feemlist(c(dir, dir2), 'panorama'),
	c(
		file.path(basename(dir), c('a.dat', 'sub/b.dat', 'sub/subsub/c.dat')),
		file.path(basename(dir2), c('d.dat', 'sub/e.dat', 'sub/subsub/f.dat'))
	)
)

# use names of the arguments where available
assertnames(
	feemlist(c(foo = dir, bar = dir2), 'panorama'),
	c(
		file.path('foo', c('a.dat', 'sub/b.dat', 'sub/subsub/c.dat')),
		file.path('bar', c('d.dat', 'sub/e.dat', 'sub/subsub/f.dat'))
	)
)

# no common components to shave off in this case
assertnames(
	feemlist(c(foo = dir, dir2), 'panorama'),
	c(
		file.path('foo', c('a.dat', 'sub/b.dat', 'sub/subsub/c.dat')),
		file.path(pfix(dir2), c('d.dat', 'sub/e.dat', 'sub/subsub/f.dat'))
	)
)

mktree(
	dir3 <- mktempd(),
	list('a.dat' = f, 'b.dat' = f)
)

# trees of different depths
assertnames(
	feemlist(c(foo = dir, bar = dir2, baz = dir3), 'panorama'),
	c(
		file.path('foo', c('a.dat', 'sub/b.dat', 'sub/subsub/c.dat')),
		file.path('bar', c('d.dat', 'sub/e.dat', 'sub/subsub/f.dat')),
		file.path('baz', c('a.dat', 'b.dat'))
	)
)

# should also work on single files
assertnames(
	feemlist(c(
		'foo' = file.path(dir3, 'a.dat'),
		'bar' = file.path(dir, 'a.dat')
	), 'panorama'),
	c('foo', 'bar')
)

# ... and on partially-named x
assertnames(
	feemlist(c(
		'foo' = file.path(dir3, 'a.dat'),
		file.path(dir, 'a.dat')
	), 'panorama'),
	c('foo', file.path(pfix(dir), 'a.dat'))
)

# let's mix files and directories, named and un-named
assertnames(
	feemlist(c(
		'foo' = dir,
		dir2,
		'bar' = file.path(dir3, 'b.dat'),
		file.path(dir3, 'a.dat')
	), 'panorama'),
	c(
		file.path('foo', c('a.dat', 'sub/b.dat', 'sub/subsub/c.dat')),
		file.path(pfix(dir2), c('d.dat', 'sub/e.dat', 'sub/subsub/f.dat')),
		'bar',
		file.path(pfix(dir3), 'a.dat')
	)
)

# try with a file of wrong format and filter it out
g <- system.file(file.path('extdata', 'panorama.txt'), package = 'albatross')

# ...using a pattern
mktree(
	dir4 <- mktempd(),
	list('a.DAT' = f, 'b.txt' = g)
)
assertnames(
	feemlist(dir4, 'panorama', pattern = '\\.dat$', ignore.case = TRUE),
	'a.DAT'
)

# ...by ignoring subdirectories and their contents
mktree(
	dir5 <- mktempd(),
	list('a.dat' = f, sub = list('b.txt' = g))
)
assertnames(
	feemlist(dir5, 'panorama', recursive = FALSE),
	'a.dat'
)

unlink(c(dir, dir2, dir3, dir4, dir5), recursive = TRUE)
