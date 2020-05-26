# Test for CellPhoneDB loading
context("Test Object Operations")
load("../testdata/decidua.small.rda")

test_that("Test Matchmaking", {
  obj.res <- Matchmaking(decidua.small, n_perm = 2, emd = FALSE, n_cores = 2)
  expect_equal(max(obj.res@strength), obj.res@strength["EVT|dM1","ACVR2B_INHBA"])
  expect_equal(obj.res@pvalue["EVT|dM1","ACVR2B_INHBA"], 0)
  obj.res <- Matchmaking(decidua.small, n_perm = 2, emd = TRUE,  nbins = 10, n_cores = 2)
  expect_equal(max(obj.res@strength), obj.res@strength["EVT|EVT","ACVR2B_INHBA"])
  expect_equal(obj.res@pvalue["EVT|EVT","ACVR2B_INHBA"], 0)
  obj.res <- Matchmaking(decidua.small, n_perm = 2, emd = TRUE, weighted = TRUE, nbins = 10, n_cores = 2)
  expect_equal(max(obj.res@strength), obj.res@strength["EVT|EVT","ACVR2B_INHBA"])
  expect_equal(obj.res@pvalue["EVT|EVT","ACVR2B_INHBA"], 0)
})

test_that("Test Selecting", {
  obj.res <- Selecting(decidua.small, filter.cells = TRUE)
  expect_equal(nrow(obj.res@selected$strength), 2)
  obj.res <- Selecting(obj.res, filter.cells = FALSE)
  expect_equal(nrow(obj.res@selected$strength), 4)
})

test_that("Test Converting", {
  convert.res <- Converting(decidua.small, selected = FALSE)
  expect_equal(nrow(convert.res), 16)
  obj.res <- Selecting(decidua.small, filter.cells = TRUE)
  convert.res <- Converting(obj.res, selected = TRUE)
  expect_equal(nrow(convert.res), 2)
})

test_that("Test Subsetting", {
  obj.res <- Subsetting(decidua.small, ident1 = "dM1", ident2 = "EVT")
  expect_equal(nrow(obj.res@strength), 2)
  obj.res <- Subsetting(decidua.small, ident1 = "dM1", ident2 = "EVT", directed = TRUE)
  expect_equal(nrow(obj.res@strength), 1)
  obj.res <- Subsetting(decidua.small, ident1 = "dM1", ident2 = "EVT", 
                        partner_a = "ACVR2A", partner_b = "INHBA", directed = TRUE)
  expect_equal(length(obj.res@strength), 1)
})

test_that("Test Complexing", {
  obj.res <- Complexing(decidua.small)
  expect_equal(ncol(obj.res@strength), 9)
  expect_equal(ncol(obj.res@strength), nrow(obj.res@interaction))
  expect_equal(obj.res@strength["EVT|dM1","ACVR1B:ACVR2B_INHBA:INHBB"],0.04915128)
})

test_that("Test Merging", {
  obj.res <- Merging(decidua.small, strength_merge = "max")
  expect_equal(nrow(obj.res@strength), 3)
  expect_equal(rownames(obj.res@strength), c("dM1|EVT","dM1|dM1","EVT|EVT"))
  expect_equal(obj.res@strength["dM1|EVT","ACVR2A_INHBA"],
               max(obj.res@misc$unmerged_strength["dM1|EVT","ACVR2A_INHBA"],obj.res@misc$unmerged_strength["EVT|dM1","ACVR2A_INHBA"]))  
})

test_that("Test Resetting", {
  obj.res <- Merging(decidua.small)
  obj.res <- Resetting(obj.res,by = "merge")
  expect_equal(obj.res@strength, decidua.small@strength)
  expect_equal(obj.res@pvalue, decidua.small@pvalue)
  obj.res <- Complexing(decidua.small)
  obj.res <- Resetting(obj.res,by = "complex")
  expect_equal(obj.res@strength, decidua.small@strength)
  expect_equal(obj.res@pvalue, decidua.small@pvalue) 
})