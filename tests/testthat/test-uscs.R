test_that("We can query UCSC", {
  expect_equal(
    sdamisc::queryUcsc("SELECT chrom,chromStart,chromEnd FROM centromeres LIMIT 1"),
    data.frame("chrom" = "chr1", "chromStart" = 122026459, "chromEnd" = 122224535)
  )
  expect_equal(
    sdamisc::queryUcsc("SELECT chrom,chromStart,chromEnd FROM centromeres WHERE chrom in (?) LIMIT 2", list("chr1")),
    data.frame("chrom" = "chr1", "chromStart" = c(122503247,122026459), "chromEnd" = c(124785432,122224535))
  )
  expect_equal(
    sdamisc::queryUcsc("SELECT DISTINCT chrom FROM centromeres WHERE chrom in (?)", list(c("chr1", "chr2"))),
    data.frame("chrom" = c("chr1", "chr2"))
  )
  # bad query
  expect_error(sdamisc::queryUcsc("SELECT chrom centromeres LIMIT 1"))
})
