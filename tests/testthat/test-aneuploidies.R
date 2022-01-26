# define some test data


# all gained
t1 <- data.frame("chromosome" = "chr1", "start" = 0, "end" = 1e4, "length" = 1e4, "direction" = "gain")
# all loss
t2 <- data.frame("chromosome" = "chr1", "start" = 0, "end" = 1e4, "length" = 1e4, "direction" = "loss")
# nothing changed
t3 <- data.frame("chromosome" = "chr1", "start" = c(0,1e5), "end" = c(1e4,2e5),
                 "length" = c(1e4,1e5), "direction" = "normal")
# 20/20/60
t4 <- data.frame("chromosome" = "chr1", "start" = c(0, 120, 220), "end" = c(20, 140, 300),
                 "length" = c(20,20,60), "direction" = c("gain", "loss", "normal"))
# 40/60, 60/40
t5 <- data.frame("chromosome" = "chr1", "start" = c(0, 140), "end" = c(40, 200),
                 "length" = c(40, 60), "direction" = c("gain", "loss"))
t6 <- data.frame("chromosome" = "chr1", "start" = c(0, 140), "end" = c(40, 200),
                 "length" = c(40, 60), "direction" = c("loss", "gain"))
# two "chromosome" case
t7 <- data.frame("chromosome" = c("chr1","chr2","chrSILLY"), "start" = c(0, 120, 220), "end" = c(20, 140, 300),
                 "length" = c(20,20,60), "direction" = c("gain", "loss", "normal"))

test_that("Fraction altered works", {


  expect_equal(sdamisc::fractionAltered(t1), list("gain" = 1, "loss" = 0))
  expect_equal(sdamisc::fractionAltered(t2), list("gain" = 0, "loss" = 1))
  expect_equal(sdamisc::fractionAltered(t3), list("gain" = 0, "loss" = 0))
  expect_equal(sdamisc::fractionAltered(t4), list("gain" = .2, "loss" = .2))
  expect_equal(sdamisc::fractionAltered(t5), list("gain" = .4, "loss" = .6))
  expect_equal(sdamisc::fractionAltered(t6), list("gain" = .6, "loss" = .4))
  expect_equal(sdamisc::fractionAltered(t7), list("gain" = .2, "loss" = .2))
  })


test_that("CallWcna works",{
  expect_equal(sdamisc::callWcna(t1, "chr1"), "gain")
  expect_warning(expect_equal(sdamisc::callWcna(t1, "chrSILLY"), NA_character_))
  expect_equal(sdamisc::callWcna(t2, "chr1"), "loss")
  expect_equal(sdamisc::callWcna(t3, "chr1"), "normal")
  expect_equal(sdamisc::callWcna(t7, "chr1"), "gain")
  expect_equal(sdamisc::callWcna(t7, "chr2"), "loss")
  expect_equal(sdamisc::callWcna(t7, "chrSILLY"), "normal")
  # with a custom threshold these change
  expect_equal(sdamisc::callWcna(t5, "chr1"), sdamisc::callWcna(t6, "chr1"))
  expect_equal(sdamisc::callWcna(t5, "chr1", threshold = 0.6),"loss")
  expect_equal(sdamisc::callWcna(t6, "chr1", threshold = 0.6),"gain")
  # this is the biased thing
  expect_equal(sdamisc::callWcna(t4, "chr1"), "normal")
  expect_warning(expect_equal(sdamisc::callWcna(t4, "chr1", threshold = 0.2), "gain"))
})

test_that("CallArmcna works", {
  # not ready
})
