.onAttach <- function(libname, pkgname){
  ## Set verbosity to 1 on package load.
  options(scigenex_verbosity = 1)
}

# Avoid Notes/Warning about ggplot2
# dataframe variables on R CMD Check.
# ... Ugly...
utils::globalVariables(c("x", 
                         "y",
                         "x1",
                         "x2",
                         "y1",
                         "y2",
                         "Set_1",
                         "Set_2",
                         "jaccard",
                         "cluster",
                         "intensity",
                         "value"))

  