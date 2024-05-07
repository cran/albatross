# Since R (>= 4.4) officially exports recvData and recvOneData and CRAN
# doesn't seem to mind packages providing methods for it for older R
# versions, we can finally hook into parallel functions and provide a
# progress bar without resorting to:
# - Importing snow internals and calling them manually (doSNOW)
# - Splitting the workload into length(cl)-sized chunks and submitting
#   them as jobs between updating the progress bar(pbapply)
# - Adding local side effects to the function and detecting them while
#   the parallel function is running in a child process (parabar)

.parallel_callback <- function(x) {
	cb <- attr(x, 'albatross_callback')
	if (!is.null(cb)) on.exit(cb())
	NextMethod()
}

maybe_default_cluster <- function(cl) {
	if (
		is.null(cl) &&
		# getDefaultCluster only appeared in 3.5.0
		('getDefaultCluster' %in% getNamespaceExports('parallel'))
	) parallel::getDefaultCluster() else cl
}

wrap_cluster <- function(cl, cb) {
	cl <- maybe_default_cluster(cl)
	if (is.null(cl)) return(cl)
	structure(
		lapply(cl, function(node) structure(node,
			class = c('albatross_progress_node', class(node)),
			albatross_callback = cb
		)),
		class = c('albatross_progress_cluster', class(cl)),
		albatross_callback = cb
	)
}

# There are node objects (but not clusters) in the wild inheriting from
# environment, so we have to undo our changes manually just in case.
# Composition over inheritance would have helped, but recvData is
# only exported in R >= 4.4, and we'd rather not import it manually.
unwrap_cluster <- function(cl) for (node in cl) {
	class(node) <- class(node)[class(node) != 'albatross_progress_node']
	attr(node, 'albatross_callback') <- NULL
}
