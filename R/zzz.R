.onLoad <- function(libname, pkgname) {
  # nocov start
  ns <- asNamespace(pkgname)
  assign("GeomAuc", define_GeomAuc(), envir = ns)
  assign("StatAuc", define_StatAuc(), envir = ns)
  # nocov end
}
