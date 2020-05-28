# Test for plot
skip_if(getRversion() > 4.0)
context("Test Plot")
load("../testdata/decidua.small.rda")
library(vdiffr)

test_that("Test spot plot", {
  test.plot <- PlotSpot(decidua.small, which.cell = 1:4, which.gene = 1:8, order = FALSE, do.return = FALSE)
  expect_doppelganger("PlotSpot Strength", test.plot)
  test.plot <- PlotSpot(decidua.small, aggregate = "gene", order = FALSE, counted = TRUE, do.return = FALSE)
  expect_doppelganger("PlotSpot By Gene", test.plot)
  test.plot <- PlotSpot(decidua.small, order = FALSE, do.return = FALSE, type = "bar")
  expect_doppelganger("PlotSpot Strength Barplot", test.plot)
})

test_that("Test heatmap plot", {
  test.plot <- function() PlotHeatmap(decidua.small, which.cell = 1:4, which.gene = 1:8, order = FALSE, mar = c(12,12))
  expect_doppelganger("PlotHeatmap Strength", test.plot)
  test.plot <- function() PlotHeatmap(decidua.small, aggregate = "gene", order = FALSE, mar = c(12,12))
  expect_doppelganger("PlotHeatmap By Gene", test.plot)
})

test_that("Test network plot", {
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

