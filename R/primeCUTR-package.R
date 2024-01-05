#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import dplyr
#' @import stringr
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom Biostrings getSeq
## usethis namespace: end

utils::globalVariables(c("POS","context","orflength","pep","peptides_msg",
                         "start","tail",".","read.delim","write.table"))

NULL
