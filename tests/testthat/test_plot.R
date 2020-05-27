# Test for plot
context("Test Plot")
load("../testdata/decidua.small.rda")
library(vdiffr)

test_that("Test spot plot", {
  test.plot <- PlotSpot(decidua.small, which.cell = 1:4, which.gene = 1:8, order = TRUE, do.return = TRUE)
  expect_is(test.plot$ggplot, "ggplot")
  expect_equal(test.plot$strength["EVT|dM1","ACVR2A_BMP2"],decidua.small@strength["EVT|dM1","ACVR2A_BMP2"])
  expect_equal(test.plot$pvalues["dM1|EVT","ACVR2A_INHBC"],decidua.small@pvalue["dM1|EVT","ACVR2A_INHBC"])
  expect_doppelganger("PlotSpot Strength", test.plot)
  test.plot <- PlotSpot(decidua.small, which.cell = 1:4, which.gene = 1:8, aggregate = "cell", counted = TRUE, do.return = TRUE)
  expect_equal(dim(test.plot$ggplot$data), c(4,3))
  expect_equal(max(test.plot$ggplot$data$Counts), 6)
  expect_doppelganger("PlotSpot By Cell", test.plot)
  test.plot <- PlotSpot(decidua.small, which.cell = 1:4, which.gene = 1:8, aggregate = "gene", counted = TRUE, do.return = TRUE)
  expect_equal(dim(test.plot$ggplot$data), c(64,3))
  expect_doppelganger("PlotSpot By Gene", test.plot)
  test.plot <- PlotSpot(decidua.small, which.cell = 1:4, which.gene = 1:8, order = TRUE, do.return = TRUE, type = "bar")
  expect_doppelganger("PlotSpot Strength Barplot", test.plot)
})

test_that("Test heatmap plot", {
  test.plot <- function() PlotHeatmap(decidua.small, which.cell = 1:4, which.gene = 1:8, order = TRUE, mar = c(12,12))
  expect_doppelganger("PlotHeatmap Strength", test.plot)
  test.plot <- function() PlotHeatmap(decidua.small, aggregate = "cell", order = TRUE, mar = c(12,12))
  expect_doppelganger("PlotHeatmap By Cell", test.plot)
  test.plot <- function() PlotHeatmap(decidua.small, aggregate = "gene", order = TRUE, mar = c(12,12))
  expect_doppelganger("PlotHeatmap By Gene", test.plot)
})

test_that("Test network plot", {
  test.plot <- function() PlotNetwork(decidua.small, aggregate = "cell")
  expect_doppelganger("PlotNetwork By Cell", test.plot)
  test.plot <- function() PlotNetwork(decidua.small, aggregate = "gene")
  expect_doppelganger("PlotNetwork By Gene", test.plot)
})

test_that("Test histogram plot", {
  test.plot <- function() PlotHistogram(decidua.small, ligand = "ACVR2A", receptor = "BMP2", ligand.ident = "EVT", receptor.ident = "dM1")
  expect_doppelganger("PlotHistogram", test.plot)
})

test_that("Test scatter plot", {
  test.plot <- function() PlotScatter(decidua.small, ident1 = "EVT", ident2 = "dM1", ligands = "CD74", receptors = "MIF", 
                                      use_raw = FALSE, add.lines = TRUE, label.offset = -3)
  expect_doppelganger("PlotScatter", test.plot)
})

