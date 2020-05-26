# Test for plot
context("Test Plot")
load("../testdata/decidua.small.rda")

test_that("Test spot plot", {
  test.plot <- PlotSpot(decidua.small, which.cell = 1:4, which.gene = 1:8, order = TRUE, do.return = TRUE)
  expect_is(test.plot$ggplot, "ggplot")
  expect_equal(test.plot$strength["EVT|dM1","ACVR2A_BMP2"],decidua.small@strength["EVT|dM1","ACVR2A_BMP2"])
  expect_equal(test.plot$pvalues["dM1|EVT","ACVR2A_INHBC"],decidua.small@pvalue["dM1|EVT","ACVR2A_INHBC"])
  test.plot <- PlotSpot(decidua.small, which.cell = 1:4, which.gene = 1:8, aggregate = "cell", counted = TRUE, do.return = TRUE)
  expect_equal(dim(test.plot$ggplot$data), c(4,3))
  expect_equal(max(test.plot$ggplot$data$Counts), 6)
  test.plot <- PlotSpot(decidua.small, which.cell = 1:4, which.gene = 1:8, aggregate = "gene", counted = TRUE, do.return = TRUE)
  expect_equal(dim(test.plot$ggplot$data), c(64,3))
})