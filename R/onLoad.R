.onLoad <- function(libname, pkgname) {
	registerS3method('recvData', 'albatross_progress_node', .parallel_callback, loadNamespace('parallel'))
	registerS3method('recvOneData', 'albatross_progress_cluster', .parallel_callback, loadNamespace('parallel'))
	bindtextdomain(paste0('plots-', pkgname), file.path(libname, pkgname, 'po'))
}
