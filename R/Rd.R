.eqn3legacy <- function() if (
	getRversion() < '4.2' || nzchar(Sys.getenv('ALBATROSS_LEGACY_EQN3'))
) 'html' else 'FALSE'

.Rdcite.env <- new.env(parent = emptyenv())
# also ensures other invariants (e.g. $cited exists)
.Rdcite.db <- function() {
	if (!exists('db', envir = .Rdcite.env))
		.Rdcite.env$db <- utils::readCitationFile(
			system.file('REFERENCES', package = 'albatross'),
			list(Encoding = 'UTF-8')
		)
	if (!exists('cited', envir = .Rdcite.env))
		.Rdcite.env$cited <- character()
	.Rdcite.env$db
}
.Rdcite <- function(key) {
	db <- .Rdcite.db()
	if (!length(db[[key]])) stop(
		dQuote(key), ' not found in REFERENCES'
	)

	# set of keys already cited
	.Rdcite.env$cited <- union(.Rdcite.env$cited, key)

	utils::cite(key, db)
}
.Rdbibliography <- function() {
	db <- .Rdcite.env$db
	cited <- .Rdcite.env$cited
	stopifnot(!is.null(db), length(cited) > 0)
	.Rdcite.env$cited <- character()

	ret <- tools::toRd(db[cited])
	# work around an Rd2HTML bug (PR18470) where newline-separated
	# paragraphs originating from a \Sexpr aren't recognised -- except
	# there's already a <p> on the first line
	paste(
		c('', rep('\\if{html}{\\out{<p>}}', length(ret)-1)), ret,
		collapse = '\n\n'
	)
}
.Rdreference <- function(key) {
	entry <- .Rdcite.db()[[key]]
	if (!length(entry)) stop(
		dQuote(key), ' not found in REFERENCES'
	)
	tools::toRd(entry)
}
