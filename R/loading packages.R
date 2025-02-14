.onAttach <- function(libname, pkgname) {
  # Attach desired packages
  if (!("scales" %in% (.packages()))) {
    suppressPackageStartupMessages(library(scales))
  }
  if (!("bipartite" %in% (.packages()))) {
    suppressPackageStartupMessages(library(bipartite))
  }
  
  if (!("mipfp" %in% (.packages()))) {
    suppressPackageStartupMessages(library(mipfp))
  }
  
  packageStartupMessage("HumPlant: mipfp, scales, and bipartite have been loaded.")
}
