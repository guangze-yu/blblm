library(testthat)
test_that("blblm()", {
  envir = environment
  fit <- blblm(mpg ~ wt + hp, mtcars, m=5, B=100)
  expect_s3_class(fit, "blblm")
  expect_equal(length(coef(fit)), 3)
  splitdata <- split_data(mtcars, m=5)
  expect_equal(length(splitdata), 5)
  fit2 <- blbgm(am ~ wt + hp, mtcars, m=5, B=100, cl = 2,family = gaussian)
  expect_s3_class(fit2, "blbgm")
  expect_equal(length(fit2$estimates), 5)

})
