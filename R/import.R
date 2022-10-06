# The code behind feem.connection / feem.character / feemlist.character.
# All file import routines are here.

# like table as data.frame, but convert factors back to original values
.table <- function(x) {
	x <- as.data.frame(table(x))
	x[,1] <- as.numeric(levels(x[,1]))[x[,1]]
	x
}

# Some formats store unmeasured data as zeroes instead of NA or anything
# like that. Check if (anti-)Stokes area is unmeasured and fill it with
# NAs instead.
.dropStokes <- function(eem, na) {
	wldiff <- outer(
		attr(eem, 'emission'), attr(eem, 'excitation'), `-`
	)
	# count occurences of each lambda.em - lambda.ex
	wldiff.all <- .table(wldiff)
	# ensure the order we need
	wldiff.all <- wldiff.all[order(wldiff.all[,1]),]
	# now count the occurences where potential NAs exist
	wldiff.na <- .table(wldiff[eem %in% na])
	if (nrow(wldiff.na) == 0) return(eem) # maybe there aren't

	# find out where counts match
	matches <- merge(wldiff.all, wldiff.na, by = 1, all.x = TRUE)
	# count of 0 would have meant absence of the row, so we can
	# safely assign it, knowing it wouldn't match
	matches[is.na(matches[,3]),3] <- 0
	matches[,4] <- matches[,2] == matches[,3]
	# only choose unbroken lambda.em-lambda.ex ranges that touch
	# the maximal or minimal value
	runs <- rle(matches[,4])
	runs$idx <- cumsum(runs$lengths)

	# collect the lambda.em-lambda.ex values to empty from matches[,1]
	toremove <- integer()
	if (runs$values[1]) { # empty anti-Stokes zone?
		toremove <- c(toremove, matches[1:runs$idx[1], 1])
	}
	if (runs$values[length(runs$values)]) { # empty Stokes zone?
		toremove <- c(toremove, matches[
			runs$idx[length(runs$idx) - 1]:nrow(matches),
			1
		])
	}

	eem[wldiff %in% toremove] <- NA
	eem
}

read.panorama <- function(filename) {
	lines <- readLines(filename)

	if (any(unlist(header <- regexec(
		# slightly awkward because regexec() returns list of (-1 or
		# integer vectors)
		paste0(
			'^\\s*', # start with optional whitespace
			'([0-9.]+)\\s+([0-9.]+)\\s+([0-9.]+)', # start, end, step
			'\\s*;\\s*', # separated by semicolon and optional whitespace
			'(Emission|Excitation)', # meaning
			'\\((columns|rows)\\)$' # terminated with axis
		),
		lines[1:2]
	)) == -1)) {
		stop('Failed to parse file header')
	}
	# found a header looking like this:
	#   240.0   650.0     1.0	; Emission(columns)
	#   255.0   355.0     5.0	; Excitation(rows)
	header <- regmatches(lines[1:2], header)

	# the fluorescence data
	intensities <- as.matrix(read.table(
		# some files just have broken decimal separators
		textConnection(gsub(',', '.', lines[-(1:2)]))
	))
	stopifnot(is.numeric(intensities))

	wavelengths <- list()
	kind <- character()
	for (desc in header) {
		# 1: full string
		# 2, 3, 4: start, end, step
		# 5: 'Excitation' / 'Emission'
		# 6: 'columns' / 'rows'
		dimnum <- switch(
			desc[6], rows = 2, columns = 1,
			stop("Unknown keyword in header")
		)

		# produce the uniform wavelength grid
		wls <- as.numeric(desc[2:4])
		wls <- seq(wls[1], wls[2], wls[3])
		stopifnot(dim(intensities)[dimnum] == length(wls))
		wavelengths[[dimnum]] <- wls

		# remember the wavelength kind <-> axis mapping
		kind[dimnum] <- desc[5]
	}

	# ensure emission wavelengths along rows, excitation along columns
	if (all(kind == c('Excitation', 'Emission'))) {
		intensities <- t(intensities)
		wavelengths[1:2] <- wavelengths[2:1]
		kind[1:2] <- kind[2:1]
	}

	# last safety check
	stopifnot(kind == c('Emission', 'Excitation'))

	.dropStokes(
		feem(intensities, wavelengths[[1]], wavelengths[[2]]),
		c(0, 100)
	)
}

read.matrix <- function(
	file, transpose = FALSE, na = NULL, fill = TRUE, fileEncoding = '', ...
) {
	# we are asked to decode non-native encoding
	# this is only possible if we are the ones to open the connection
	if (nzchar(fileEncoding)) {
		stopifnot(is.character(file))
		conn <- file(file, encoding = fileEncoding)
		on.exit(close(conn))
		# also, only readLines and scan return UTF-8;
		# read.table uses a function which doesn't
		file <- textConnection(
			readLines(conn),
			encoding = 'UTF-8' # also ask textConnection to return UTF-8
		)
	}

	# read the table itself
	x <- unname(as.matrix(read.table(
		file = file, header = FALSE, colClasses = 'character', fill = fill, ...
	)))
	# warning is rased on text -> numeric conversion resulting in NAs,
	# which we expect to happen
	suppressWarnings(mode(x) <- 'numeric')
	stopifnot(any(!is.na(x)))

	# first row and column should contain wavelengths
	wavelengths <- list(
		x[-1, 1], # the first column describes the rows
		x[1,] # the first row describes the columns
	)
	x <- x[-1, -1]

	# the first row is problematic
	# <text> <sep> <wavelength> <sep> <wavelength> ... <eol> => c(NA, wl... wl)
	# also, when sep is not missing ("whitespace")
	# <sep> <wavelength> <sep> <wavelength> ... <eol> => c(NA, wl... wl)
	# but when sep is missing, leading whitespace is skipped, resulting in
	# a short row padded with NA at the end => c(wl... wl, NA)
	wavelengths[[2]] <- if (is.na(rev(wavelengths[[2]])[1])) {
		wavelengths[[2]][-length(wavelengths[[2]])]
	} else wavelengths[[2]][-1]

	if (transpose) {
		x <- t(x)
		wavelengths[1:2] <- wavelengths[2:1]
	}

	# drop whole NA rows/cols, NA wavelengths
	keep.em <- !apply(is.na(x), 1, all) & !is.na(wavelengths[[1]])
	keep.ex <- !apply(is.na(x), 2, all) & !is.na(wavelengths[[2]])

	# additionally mark unmeasured areas as NA if possible
	.dropStokes(
		feem(
			x[keep.em, keep.ex],
			wavelengths[[1]][keep.em],
			wavelengths[[2]][keep.ex]
		),
		na
	)
}

read.F900txt <- function(f, fileEncoding = 'latin1', sep = ';') {
	# readLines doesn't support arbitrary encodings by itself
	f <- file(f, encoding = fileEncoding)
	on.exit(close(f))
	lines <- readLines(f)

	blocks <- rle(lines == '')
	blocks$end <- cumsum(blocks$lengths)
	# assuming that the file consists of three blocks separated by empty lines
	stopifnot(identical(blocks$values, c(FALSE, TRUE, FALSE, TRUE, FALSE)))

	# first block is typically one line with the sample name
	label <- lines[1:blocks$lengths[1]]

	# Second block is the metadata, a horizontal table which we'll
	# transform into a data.frame. Transposing a data.frame is a bit
	# painful.
	metadata <- t(read.table(
		text = lines[(blocks$end[2]+1):blocks$end[3]], sep = sep,
		colClasses = 'character', header = FALSE, na.strings = ''
	)) # returns a matrix
	colnames(metadata) <- metadata[1,]
	metadata <- as.data.frame(
		metadata[-1,, drop = FALSE], stringsAsFactors = FALSE, check.names = FALSE
	)
	for (i in seq_along(metadata)) metadata[[i]] <- type.convert(
		metadata[[i]], na.strings = character(0), as.is = TRUE
	)

	# last column is typically empty
	if (all(is.na(metadata[nrow(metadata),])))
		metadata <- metadata[-nrow(metadata),]

	# we only support "emission maps"
	stopifnot(metadata[['Type']] == 'Emission Scan')

	# third block is the data
	data <- read.table(text = lines[(blocks$end[4]+1):blocks$end[5]], sep = ';')
	# first column of the data is the scan wavelength
	em <- data[[1]]; data <- data[-1]
	ex <- metadata[['Fixed/Offset']]
	# the last column is typically empty
	if (all(is.na(data[ncol(data)]))) data <- data[-ncol(data)]

	structure(
		albatross::feem(as.matrix(data), em, ex),
		F900.metadata = metadata,
		F900.label = label
	)
}
