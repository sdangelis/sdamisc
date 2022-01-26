test_that("getName works", {
  expect_equal(
    purrr::map_chr(
      c(
        "xyz.txt", "./xyz.txt", "xyz", "~/xyz.txt", "../xyz.txt",
        "c:/abc/def/xyz.txt", "../abc/xyz.txt.gz"
      ),
      \(x) sdamisc::getName(x)
    ), rep("xyz", 7)
  )
  expect_equal(getName("c:\\xyz.txt"), "c:\\xyz")
})

df1 <- tibble(
  "chr" = rep("chr1", 4),
  "start" = c(100, 300, 600, 800),
  "end" = c(200, 400, 700, 900)
)
bp1 <- tibble(
  "chr" = rep("chr1", 8),
  "start" = c(100, 199, 300, 399, 600, 699, 800, 899),
  "end" = c(101, 200, 301, 400, 601, 700, 801, 900)
)

df2 <- tibble(df1, "meta" = c("record1", "record2", "record3", "record4"))
bp2 <- tibble(bp1, "meta" = c(
  "record1", "record1", "record2", "record2",
  "record3", "record3", "record4", "record4"
))

test_that("splitBP works", {
  expect_equal(sdamisc::splitBp(df1), bp1)
  expect_equal(sdamisc::splitBp(df2), bp2)
})
