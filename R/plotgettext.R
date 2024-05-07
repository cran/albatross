pgt <- function(msgid, translate)
	gettext(msgid, domain = if (translate) 'plots-albatross' else NA)
pgtq <- function(msgid, translate)
	parse(text = pgt(msgid, translate))[[1]]

# a few strings are produced dynamically, force their translation here
if (FALSE) {
pgt('Emission')
pgt('Excitation')
pgt('Fluorescence')
pgt('Scattering')
}
