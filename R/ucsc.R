# USCS + SQL --------------------------------------------------------------
#' deal with the UCSC genome browser mySQL server

queryUcsc <- \(query, params = NULL, genome = "hg38") {
  #' wrapper around dbGetQuery to query the UCSC genome browser mySQL server
  #' @param genome a UCSC genome that represents the server that will connect to
  #' @param query a SELECT SQL query to run
  #' @param params a list of parameters for the query
  #' @returns the result of your query
  #' @export
  con_ucsc <- DBI::dbConnect(RMariaDB::MariaDB(),
    db = genome, user = "genome",
    host = "genome-mysql.soe.ucsc.edu"
  )
  tryCatch({
      if (is.null(params)) {
        query <- DBI::dbGetQuery(con_ucsc, query)
      } else {
        query <- DBI::dbGetQuery(con_ucsc, query, params)
      }
    },
    finally = {
      DBI::dbDisconnect(con_ucsc)
      return(query)
    }
  )
}

# Exporting data is unorthodox - but this is tiny and almost universally needed

#' hg38 centromeres
#' cached human hg38 centromeres, using the median start/end position recorded
#' in
#' @format A data frame with 24 rows and 3 variables:
#' \describe{
#'    \item{chrom}{chromosome name - UCSC naming "chr1"},
#'    \item{chromStart}{starting position - 0-indexed},
#'    \item{chromEnd}{ending position - 0-indexed}
#' }
#' @source USCS Genome Broswer
#' @export
hg38_centromeres <- queryUcsc(
  "SELECT chrom,chromStart,ChromEnd FROM centromeres WHERE chrom in (?)",
  list(c(paste0("chr", 1:22), "chrY", "chrX"))
) |>
  dplyr::group_by(chrom) |>
  dplyr::summarise(chromStart = median(chromStart), chromEnd = median(ChromEnd))


#' hg38 telomeres
#' cached human hg38 telomeric regions using the median start/end position
#' @format A data frame with 48 rows and 3 variables:
#' \describe{
#'    \item{chrom}{chromosome name - UCSC naming "chr1"},
#'    \item{chromStart}{starting position - 0-indexed},
#'    \item{chromEnd}{ending position - 0-indexed}
#' }
#' @source USCS Genome Broswer
#' @export
hg38_telomeres <- queryUcsc(
  "SELECT chrom,chromStart,chromEnd FROM gap WHERE type IN (?)",
  list("telomere")
)

bands <- queryUcsc(
  "SELECT chrom,chromStart,chromEnd, name FROM cytoBand WHERE chrom IN (?)",
  list(c(paste0("chr", 1:22), "chrY", "chrX"))
) |>
  dplyr::mutate(arm = stringr::str_sub(name, 1, 1)) |>
  dplyr::group_by(arm, chrom) |>
  dplyr::summarise(chromStart = median(chromStart), chromEnd = median(chromEnd))
