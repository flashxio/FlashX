library("FlashR")
context("Reliability Distribution")

test_that("rdf errors for incorrect dimensions",{


})

test_that("rdf errors for incorrect data types", {

  expect_error(rdf('a string', ids))

})

test_that("rdf ranks properly in simple case",{
  A <- diag(4)-diag(4)
  A[upper.tri(A)] <- c(1, 4, 5, 5, 4, 2)
  A <- A+t(A)

  expect_true(all(rdf(A, c(1, 1, 2, 2)) == c(1.0, 1.0, 1.0, 1.0)))
  expect_true(all(rdf(A, c(1, 2, 1, 2)) == c(0.5, 0.5, 0.5, 0.5)))
})

test_that("tied ranks resolve in fractions", {
  A <- diag(4)-diag(4)
  A[upper.tri(A)] <- c(2, 2, 4, 5, 4, 2)
  A <- A+t(A)

  expect_true(all(rdf(A, c(1, 1, 2, 2)) == c(0.75, 1.0, 0.75, 1.0)))
})

test_that("mnr is mean of rdf",{


})

test_that("infinity and NaN are considered valid entries", {


})
