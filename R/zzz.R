#' @useDynLib CircSpaceTime
#' @importFrom Rcpp sourceCpp
NULL

.onAttach <- function(...) {
  tip <- c("CircSpaceTime: \n an R package for spatial and spatio-temporal modeling of Circular data \n
           Please if you are new to the package read the vignette\n
           Bugs and Issues at https://github.com/santoroma/CircSpaceTime/issues")
    packageStartupMessage(paste(strwrap(tip), collapse = "\n"))
}



