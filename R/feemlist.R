# Currently, feemlist is only about importing FEEMs from collections of
# files (or objects from other packages) and returning them in a list.

feemlist <- function(x, ...) UseMethod('feemlist')

# split single path into a vector of components
# cannot just use .Platform$path.sep because it's only a forward slash
# on Windows, but paths can and do contain backslashes too
.splitpath <- function(x) strsplit(
	x,
	if (.Platform$OS.type == 'windows') '/|\\\\' else .Platform$file.sep,
	.Platform$OS.type != 'windows'
)

feemlist.character <- function(
	x, format, pattern = NULL, recursive = TRUE,
	ignore.case = FALSE, simplify.names = TRUE, ...
) {
	# all supplied paths must exist
	stopifnot(file.exists(x))
	# provide names if vector is unnamed...
	if (is.null(names(x))) names(x) <- x
	# ...or of names are empty
	names(x)[names(x) == ''] <- x[names(x) == '']
	# any directories should be replaced by their contents, appropriately named
	isdir <- file.info(x)$isdir
	x <- c(x[!isdir], unlist(lapply(seq_along(x)[isdir], function(i) {
		# iterate over indices in x to preserve names
		dir <- x[i]
		files <- list.files(
			dir, pattern = pattern, full.names = FALSE,
			recursive = recursive, ignore.case = ignore.case
		)
		# directories are always included in non-recursive listings, skip them
		files <- files[!file.info(file.path(dir, files))$isdir]
		# rebuild the names here manually instead of letting unlist()
		# create them for us
		setNames(file.path(dir, files), nm = file.path(names(dir), files))
	})))
	if (simplify.names) {
		# split names into path components
		nc <- .splitpath(names(x))
		while (
			# need at least one remaining component
			all(lengths(nc) > 1) &&
			# otherwise try to drop same-valued components
			length(unique(vapply(nc, `[[`, character(1), 1))) == 1
		) nc <- lapply(nc, `[`, -1)
		names(x) <- vapply(nc, paste, character(1), collapse = .Platform$file.sep)
	}
	if (is.character(format)) {
		fname <- format
		format <- function(f, ...) feem(f, fname, ...)
	}
	lapply(.wraplist(x), .wrapfun(format), ...)
}

# Enhances: eemR
feemlist.eemlist <- function(x, ...) {
	stopifnot(length(list(...)) == 0)
	lapply(
		setNames(x, vapply(x, `[[`, character(1), 'sample')),
		function(eem) feem(eem$x, eem$em, eem$ex)
	)
}

# Enhances: EEM
feemlist.EEM <- function(x, ...) {
	stopifnot(length(list(...)) == 0)
	lapply(x, function(eem)
		feem(eem, as.numeric(rownames(eem)), as.numeric(colnames(eem)))
	)
}
