library(tools)
for (p in c('.', '../00_pkg_src/albatross'))
	if (file.exists(file.path(p, 'DESCRIPTION'))) {
		srcdir <- p
		break
	}
if (!exists('srcdir')) stop('Cannot figure out where my source code is')

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

walk(lapply(
	list.files(
		file.path(srcdir, 'man'),
		full.names = TRUE, pattern = '\\.Rd$', recursive = FALSE
	),
	parse_Rd, macros = loadPkgRdMacros(srcdir)
))

if (bad) stop('Problems found, see above')
