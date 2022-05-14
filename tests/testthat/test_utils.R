# Test for utility functions
context("Test for utility functions")

test_that("Zero one normalization", {
  test.data <- 1:10
  test.res <- zeroOne(test.data)
  expect_equal(max(test.res), 1)
  expect_equal(min(test.res), 0)
})

test_that("Sparse matrix conversion", {
  test.res <- asSparse(matrix(1:9,3), 0.7)
  expect_is(test.res, "matrix")
  test.res <- asSparse(matrix(c(rep(0,5),1:4),3), 0.5)
  expect_is(test.res, c("dgCMatrix", "dsCMatrix"))
})

test_that("Filter matrix", {
  test.res <- filterMatrix(matrix(1:9,3),low.cutoff = 3,high.cutoff = 6)
  expect_equal(test.res[,1,drop = TRUE], 4:6)
})

test_that("Aggregate by name", {
  mat <- matrix(1:9,3,dimnames = list(paste(letters[1:3],1:3,sep = "|"), paste(LETTERS[1:3],1:3,sep = "-")))   
  test.res <- aggregateName(mat, aggregate = "row")
  expect_equal(rownames(test.res), c(1:3,letters[1:3]))
  expect_equal(colnames(test.res), c(1:3,letters[1:3]))
  test.res <- aggregateName(mat, aggregate = "col", split.by = "-")
  expect_equal(rownames(test.res), c(1:3,LETTERS[1:3]))
  expect_equal(colnames(test.res), c(1:3,LETTERS[1:3]))
  test.res <- aggregateName(mat, aggregate = "row", cutoff = 10)
  expect_true(all(test.res == 0))
  test.res <- aggregateName(mat, aggregate = "col", split.by = "-", counted = TRUE)
  expect_equal(test.res["A", 1], 3)
  expect_equal(test.res["B", 2], 3)
  expect_equal(test.res["C", 3], 3)
})

test_that("Random identity generation", {
  test.res <- randomIdents(rep(1:2, each = 10), n = 1, seed = 123)
  set.seed(123)
  sample.id <- sample(x = rep(1:2, each = 10))
  expect_equal(test.res[,1,drop = TRUE], sample.id)
})

test_that("Order matrix", {
  test.res <- orderMatrix(matrix(1:9,3))
  expect_equal(test.res[1,1], 9)
  expect_equal(test.res[3,3], 1)
})
            
test_that("Merge lists", {
  test.res <- mergeLists(list(a = 1), list(b = 2))
  expect_equal(test.res$a, 1)
  expect_equal(test.res$b, 2)
})

test_that("Command Log", {
  foo <- function(object = "a", b = 2, test = TRUE){
    logCommand()
  }
  test.log <- foo()
  expect_equal(test.log$b, 2)
  expect_equal(test.log$test, TRUE)
  expect_equal(format(test.log$Time,"%Y-%m-%d"), format(Sys.time(),"%Y-%m-%d"))   
})

m <- matrix(rnorm(10), 5, 2, dimnames = list(1:5, 1:2))
m
orderMatrix(m)
## permute rows and columns
seriation::permute(m, seriation::ser_permutation(5:1, 2:1))
