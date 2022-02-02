getName <- function(path) {
  #' extracts a file name from a path with leading /s and trailing extension
  #'
  #' @param path a file path
  #' @returns a string with the file name between the last / and first `.`
  #' @export
  stringr::str_split(string = path, pattern = "/") |>
    dplyr::last() |>
    dplyr::last() |>
    stringr::str_split(pattern = "\\.") |>
    dplyr::first() |>
    dplyr::first()
}

#' Split each record in a bedfile df into 2 breakpoints
#'
#' @param bed a bed file as a data frame. column names are unimportant
#' but columns 1:3 must follow bed order (chromosome, start, end)
#' @return returns a bed file with
#' @note Each record is split in two records of length 1:
#' Record 1 is start to start+1 and record 2 is end-1 to end.
#' @export
splitBp <- \(bed) {
  bed |>
    dplyr::mutate(!!names(bed)[3] := get(names(bed)[3]) - 1) |>
    tidyr::pivot_longer(colnames(bed[2:3])) |>
    dplyr::select(-name) |>
    dplyr::mutate(end = value + 1) |>
    dplyr::rename(start = value) |>
    dplyr::relocate(!!colnames(bed[1]), start, end)
}
