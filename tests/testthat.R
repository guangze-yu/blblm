library(testthat)
library(blblm)

test_check("blblm")
envir = environment()
fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
expect_s3_class(fit, "blblm")

expect_equal(co(fit),n(nrow(mtcars,3)))
