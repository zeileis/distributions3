.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)
  assign("GeomAuc", define_GeomAuc(), envir = ns)
  assign("StatAuc", define_StatAuc(), envir = ns)
}
