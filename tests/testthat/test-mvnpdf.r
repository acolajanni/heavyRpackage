test_that("correct result for univariate gaussian", {
  expect_equal(mvnpdf(x=matrix(1.96), Log=FALSE)$y, dnorm(1.96))
  expect_equal(mvnpdf(x=matrix(c(1.96, -0.5), ncol = 2), Log=FALSE)$y,
               dnorm(c(1.96, -0.5)))
})

test_that("correct results for bivariate gaussian", {
  expect_equal(mvnpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1), Log=FALSE)$y,
               mvtnorm::dmvnorm(rep(1.96, 2)))
})
