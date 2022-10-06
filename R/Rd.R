.eqn3legacy <- function() if (
	getRversion() < '4.2' || nzchar(Sys.getenv('ALBATROSS_LEGACY_EQN3'))
) 'html' else 'FALSE'
