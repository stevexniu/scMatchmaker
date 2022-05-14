# Test for utility functions
context("Test for Matchmaker class")
load("../testdata/decidua.small.rda")
dir.name <- "../testdata/cellphonedb/"
cellphonedb.data <- LoadCellPhoneDB(interaction_input = paste(dir.name, "interaction_input.txt", sep = ""),
                                    gene_input = paste(dir.name ,"gene_input.txt", sep = ""),
                                    complex_input = paste(dir.name, "complex_input.txt", sep = ""),
                                    path = TRUE)

test_that("Test Matchmaker class", {
  decidua.interaction <- Screening(data = decidua.small@data, annotation = decidua.small@annotation, 
                                   interaction = cellphonedb.data, 
                                   project_name = "Decidua")
  expect_is(decidua.interaction@data, "matrices")
  expect_is(decidua.interaction@annotation, "data.frame")
  expect_is(decidua.interaction@interaction, "data.frame")
  expect_true(is(decidua.interaction@strength, "matrices"))
  expect_true(is(decidua.interaction@pvalue, "matrices"))
  expect_is(decidua.interaction@selected, "list")
  expect_is(decidua.interaction@project_name, "character")
  expect_is(decidua.interaction@command, "list")
  expect_is(decidua.interaction@misc, "list")
})



