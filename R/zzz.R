.onAttach <- function(libname, pkgname){
  ## Set verbosity to 1 on package load.
  options(scigenex_verbosity = 1)
}

# .First.lib <- function(lib,pkg) {
#    library.dynam("scigenex",pkg,lib)
#}

  