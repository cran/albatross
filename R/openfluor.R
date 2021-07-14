.maybe.rescale <- function(x, shift, scale) {
	if (shift) x <- x - min(x)
	if (scale) x <- x / max(x)
	x
}

write.openfluor <- function(
	model, filename, name = '?', creator = '?', doi = '?',
	reference = '?', unit = '?', toolbox = NULL, date = NULL,
	fluorometer = '?', constraints = NULL, validation = '?',
	methods = '?', preprocess = '?', sources = '?', ecozones = '?',
	description = '', shift = FALSE, scale = TRUE
) {
	stopifnot(inherits(model, 'feemparafac'))

	# defaults we can provide
	if (missing(toolbox)) toolbox <- paste(
		vapply(
			c('albatross', 'multiway'),
			function(pkg) paste(pkg, packageVersion(pkg)),
			character(1)
		), collapse = ', '
	)
	if (missing(date)) date <- format(Sys.Date(), '%Y-%m-%d')
	if (missing(constraints)) constraints <- model$const
	stopifnot(nchar(description) <= 256)
	cube <- feemcube(model)

	fh <- file(filename, 'wt')
	on.exit(close(fh))

	cat(
		'#\n# Fluorescence Model\n#\n',
		'name',        '\t', name,                                '\n',
		'creator',     '\t', creator,                             '\n',
		'doi',         '\t', doi,                                 '\n',
		'reference',   '\t', reference,                           '\n',
		'unit',        '\t', unit,                                '\n',
		'toolbox',     '\t', toolbox,                             '\n',
		'date',        '\t', date,                                '\n',
		'fluorometer', '\t', fluorometer,                         '\n',
		'nSample',     '\t', dim(cube)[3],                        '\n',
		'constraints', '\t', paste(constraints, collapse = ', '), '\n',
		'validation',  '\t', validation,                          '\n',
		'methods',     '\t', paste(methods,     collapse = ', '), '\n',
		'preprocess',  '\t', paste(preprocess,  collapse = ', '), '\n',
		'sources',     '\t', paste(sources,     collapse = ', '), '\n',
		'ecozones',    '\t', paste(ecozones,    collapse = ', '), '\n',
		'description', '\t', description,                         '\n',
		'#\n# Excitation/Emission (Ex, Em), ',
		'wavelength (nm), component[n] (intensity)\n#\n',
		file = fh, sep = ''
	)
	for (ds in list(
		cbind(
			'Ex', attr(cube, 'excitation'),
			.maybe.rescale(model$B, shift, scale)
		), cbind(
			'Em', attr(cube, 'emission'),
			.maybe.rescale(model$A, shift, scale)
		)
	)) write.table(
		ds, file = fh, sep = '\t', row.names = FALSE,
		col.names = FALSE, quote = FALSE
	)
}
