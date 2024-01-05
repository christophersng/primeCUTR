#' Get path to primeCUTR example
#'
#' primeCUTR comes bundled with some example VCFs in its `inst/extdata`
#' directory. This function make them easy to access.
#'
#' @param path Name of file. If `NULL`, the example files will be listed.
#' @export
#' @examples
#' primeCUTR_example()
#' primeCUTR_example("vep_HCC1937.vcf")
primeCUTR_example <- function(path = NULL) {
  if (is.null(path)) {
    dir(system.file("extdata", package = "primeCUTR"))
  } else {
    system.file("extdata", path, package = "primeCUTR", mustWork = TRUE)
  }
}
