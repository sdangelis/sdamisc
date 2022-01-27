#' Julia interop is tested manually, due to the need to install Julia,
#' initialise Julia
#' and potentially debug errors on the Julia side

t4 <- data.frame(
  "chromosome" = "chr1", "start" = c(0L, 120L, 220L), "end" = c(20L, 140L, 300L),
  "length" = c(20, 20, 60), "direction" = c("gain", "loss", "normal")
)
t5 <- data.frame(
  "chromosome" = t4[1], "start" = as.double(t4$start),
  "end" = as.double(t4$end)
)

test_that("We", {
  JuliaCall::julia_setup()
  JuliaCall::julia_install_package_if_needed("GenomicFeatures")
  JuliaCall::julia_library("GenomicFeatures")
  expect_invisible(sdamisc::bed2Julia(segs = t4, name = "t4"))
  expect_true(JuliaCall::julia_eval('t4==IntervalCollection(
                                    [Interval("chr1", 0, 20),
                                     Interval("chr1", 120, 140),
                                     Interval("chr1",220,300)])'))
  # can we coherce doubles to int fine. note idk what happens if we
  # pass a rational. But honestly why would you pass a rational there?
  expect_invisible(sdamisc::bed2Julia(segs = t5, name = "t5"))
  expect_true(JuliaCall::julia_eval("t4 ==t5"))
})
