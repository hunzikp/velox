.onUnload <- function (libpath) {
  library.dynam.unload("velox", libpath)
}

Rcpp::loadModule("BOOSTGEOM", TRUE)
