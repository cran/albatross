library(tools)

bad <- FALSE
walk <- function(x) {
	# Use \link[pkg]{topic} for external help, \link{topic} for internal help
	if (identical(attr(x, 'Rd_tag'), '\\link')) {
		topic <- as.vector(x[[1]])
		package <- as.vector(attr(x, 'Rd_option'))
		if (is.null(package)) package <- 'albatross'

		h <- help((topic), (package))

		if (length(h) != 1) {
			bad <<- TRUE
			message(
				getSrcFilename(x), ':', getSrcLocation(x), ': ',
				getSrcref(x), '\n',
				'Found ', length(h), ' results for help(',
				topic, ', ', package, ') instead of 1:', '\n',
				deparse(as.vector(h))
			)
		}
	}
	if (is.list(x)) for(tag in x) Recall(tag)
}

walk(Rd_db('albatross'))

if (is.na(packageVersion('albatross')[[1,4]])) {
	# must be release testing
	n <- news(package = 'albatross')
	for (column in c('Date', 'Version', 'Text'))
	if (any(i <- is.na(n[,column]))) {
		i <- which(i)
		message(
			'NEWS: column ', sQuote(column), ' is NA in row[s] ',
			paste(i, collapse = ', ')
		)
		bad <<- TRUE
	}
}

if (bad) stop('Problems found, see above')
