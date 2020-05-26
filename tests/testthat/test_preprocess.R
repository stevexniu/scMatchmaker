# Test for CellPhoneDB loading
context("Load CellPhoneDB")

dir.name <- "../testdata/cellphonedb/"
cellphonedb.data <- LoadCellPhoneDB(interaction_input = paste(dir.name, "interaction_input.txt", sep = ""),
                             gene_input = paste(dir.name ,"gene_input.txt", sep = ""),
                             complex_input = paste(dir.name, "complex_input.txt", sep = ""),
                             path = TRUE)

test_that("CellPhoneDB Data Parsing", {
  expect_is(cellphonedb.data, "data.frame")
  expect_equal(ncol(cellphonedb.data), 15)
  expect_equal(nrow(cellphonedb.data), 1574)
  expect_equal(cellphonedb.data[1,1], "ACE2")
  expect_equal(cellphonedb.data[16,11], "ACVR2A")
  expect_equal(cellphonedb.data[3,13], NA_character_)
})

# Test for data preprocessing
context("Data Preprocessing")
library(Matrix)
library(stats)

mat <- matrix(1:9, nrow = 3)
mat.sp <- Matrix(1:9, nrow = 3, sparse = TRUE)

test_that("LogTPM Normalization", {
  mat.norm <- Normalization(mat, normalization = "logTPM")
  mat.norm.sp <- Normalization(mat.sp, normalization = "logTPM")
  mat.logTPM <- logTPM(mat)
  expect_equal(mat.norm, mat.logTPM)
  expect_equal(mat.norm, as.matrix(mat.norm.sp), check.attributes = FALSE)
  expect_equal(mat.norm[1,1], log(1/sum(1:3)*1e4+1,2))
})

test_that("Cosine Normalization", {
  mat.norm <- Normalization(mat, normalization = "cosine")
  mat.norm.sp <- Normalization(mat.sp, normalization = "cosine")
  mat.cosine <- CosineNorm(mat, display_progress = FALSE)
  mat.cosine.sp <- CosineNormSparse(mat.sp, display_progress = FALSE)
  expect_equal(mat.norm, mat.cosine)
  expect_equal(mat.norm, as.matrix(mat.norm.sp), check.attributes = FALSE)
  expect_equal(mat.norm, as.matrix(mat.cosine.sp), check.attributes = FALSE)
  expect_equal(mat.norm[1,1], 1/sum(c(1:3)^2)^0.5)
})

test_that("Scaling", {
  mat.scaled <- Scaling(mat)
  mat.scaled.sp <- Scaling(mat.sp)
  mat.scaled.r <- t(scale(t(mat)))
  expect_equal(mat.scaled, mat.scaled.sp, check.attributes = FALSE)
  expect_equal(mat.scaled, mat.scaled.r, check.attributes = FALSE)
  expect_equal(rowMeans(mat.scaled), rep(0,nrow(mat.scaled)))
  expect_equal(apply(mat.scaled, 1, sd), rep(1,nrow(mat.scaled)))
})

# Test for object creation
context("Object Creation")

test.data <- decidua[,c(1:50,551:600)]
test.meta <- decidua.annotation$annotation[c(1:50,551:600)]
test.interactions <- cellphonedb[1:10,]
object <- Screening(data = test.data, annotation = test.meta, interaction = test.interactions, project_name = "test")

test_that("Object created correctly", {
  expect_is(object, "Matchmaker")
  expect_equal(ncol(object@data), nrow(object@annotation))
  expect_equal(sort(rownames(object@data)), sort(unique(unlist(test.interactions[,1:2]))))
})

test_that("Sketching correctly", {
  size <- 10
  object.new <- Sketching(object, size = size, seed = 123)
  set.seed(123)
  sample.id <- c(sample(1:50,5), sample(51:100,5))
  expect_equal(object.new@misc$sketch_id, sample.id)
})
