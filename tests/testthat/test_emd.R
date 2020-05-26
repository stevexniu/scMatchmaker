# Test for base calculation
context("Base Model Calculation")
library(stats)

test.data <- scMatchmaker::decidua[,c(1:50,551:600)]
test.meta <- scMatchmaker::decidua.annotation$annotation[c(1:50,551:600)]
test.interactions <- cellphonedb[1:10,]
object <- Screening(data = test.data, annotation = test.meta, interaction = test.interactions, project_name = "test")
test.data <- object@data

test_that("Stats calculated correctly", {
  test.stats <- statsCells(test.data, test.meta, stats = "mean")
  test.r <- cbind(apply(test.data[,1:50],1,mean), apply(test.data[,51:100],1,mean))
  expect_equal(test.stats, test.r, check.attributes = FALSE)
  test.stats <- statsCells(test.data, test.meta, stats = "median")
  test.r <- cbind(apply(test.data[,1:50],1,median), apply(test.data[,51:100],1,median))
  expect_equal(test.stats, test.r, check.attributes = FALSE)
})

test_that("Strengths calculated correctly", {
  # f = "+"
  test.stats <- statsCells(test.data, test.meta, stats = "mean")
  test.strengths <- calcStrength(test.stats[test.interactions[,1],], test.stats[test.interactions[,2],], 
                                                f = "+", threshold = -Inf)
  test.strengths.single <- calcStrength.single(test.stats[test.interactions[3,1],1,drop = FALSE], test.stats[test.interactions[3,2],1,drop = FALSE], 
                                     f = "+", threshold = -Inf)
  test.strengths.single.r <- mean(c(test.stats[test.interactions[3,1],1], test.stats[test.interactions[3,2],1]))
  expect_equal(test.strengths[1,3], test.strengths.single, check.attributes = FALSE)
  expect_equal(test.strengths.single, test.strengths.single.r, check.attributes = FALSE)
  # findInteractions
  test.interaction <- findInteractions.single(test.data,  test.meta, test.interactions[,1], test.interactions[,2], threshold = -Inf)
  expect_equal(test.interaction[1,3], test.strengths.single.r, check.attributes = FALSE)
  # f = "*"
  test.strengths <- calcStrength(test.stats[test.interactions[,1],], test.stats[test.interactions[,2],], 
                                                f = "*", threshold = -Inf)
  test.strengths.single <- calcStrength.single(test.stats[test.interactions[3,1],1,drop = FALSE], test.stats[test.interactions[3,2],1,drop = FALSE], 
                                                              f = "*", threshold = -Inf)
  test.strengths.single.r <- prod(c(test.stats[test.interactions[3,1],1], test.stats[test.interactions[3,2],1]))
  expect_equal(test.strengths[1,3], test.strengths.single, check.attributes = FALSE)
  expect_equal(test.strengths.single, test.strengths.single.r, check.attributes = FALSE)
  test.interaction <- findInteractions.single(test.data,  test.meta, test.interactions[,1], test.interactions[,2], pair.fxn = "*", threshold = -Inf)
  expect_equal(test.interaction[1,3], test.strengths.single.r, check.attributes = FALSE)
})

# Test for EMD calculation
context("EMD Model Calculation")
library(emdist)

test_that("Density calculated correctly", {
  test.dens <- calcDensity(1:10, breaks = 0:10)
  expect_equal(test.dens, rep(0.1,10))
})

test_that("EMD calculated correctly", {
  # unweighted EMD
  nbreaks <- seq(from = min(test.data, na.rm = TRUE), to = max(test.data, na.rm = TRUE), length.out = 10 + 1)
  test.dens1 <- calcDensity(test.data[test.interactions[3,1],1:50], breaks = nbreaks)
  test.dens2 <- calcDensity(test.data[test.interactions[3,2],1:50], breaks = nbreaks)
  emd.r <- emd2d(as.matrix(test.dens1), as.matrix(test.dens2))
  test.dens1 <- calcDensity(test.data[test.interactions[1,1],1:50], breaks = nbreaks)
  test.dens2 <- calcDensity(test.data[test.interactions[1,2],51:100], breaks = nbreaks)
  emd.max <- emd2d(as.matrix(test.dens1), as.matrix(test.dens2))
  test.dens1 <- calcDensity(test.data[test.interactions[4,1],1:50], breaks = nbreaks)
  test.dens2 <- calcDensity(test.data[test.interactions[4,2],1:50], breaks = nbreaks)
  emd.min <- emd2d(as.matrix(test.dens1), as.matrix(test.dens2))
  emd.r <- 1 -  (emd.r - emd.max) / (emd.min - emd.max)
  test.emd <- calcEMD(test.data, test.meta, test.interactions[,1], test.interactions[,2], nbins = 10, weighted = FALSE)
  expect_equal(test.emd[1,3], emd.r)
  # weighted EMD
  weight_breaks <- matrix(data = rev(x = nbreaks[-1]))
  test.dens1 <- calcDensity(test.data[test.interactions[3,1],1:50], breaks = nbreaks)
  test.dens2 <- calcDensity(test.data[test.interactions[3,2],1:50], breaks = nbreaks)
  emd.r <- emdw(as.matrix(test.dens1), weight_breaks, as.matrix(test.dens2), weight_breaks)
  test.dens1 <- calcDensity(test.data[test.interactions[1,1],1:50], breaks = nbreaks)
  test.dens2 <- calcDensity(test.data[test.interactions[1,2],51:100], breaks = nbreaks)
  emd.max <- emdw(as.matrix(test.dens1), weight_breaks, as.matrix(test.dens2), weight_breaks)
  test.dens1 <- calcDensity(test.data[test.interactions[4,1],1:50], breaks = nbreaks)
  test.dens2 <- calcDensity(test.data[test.interactions[4,2],51:100], breaks = nbreaks)
  emd.min <- emdw(as.matrix(test.dens1), weight_breaks, as.matrix(test.dens2), weight_breaks)
  emd.r <- 1 -  (emd.r - emd.max) / (emd.min - emd.max)
  test.emd <- calcEMD(test.data, test.meta, test.interactions[,1], test.interactions[,2], nbins = 10, weighted = TRUE)
  expect_equal(test.emd[1,3], emd.r)
})

