test_that("ABT fits a gaussian model", {
  set.seed(1)
  d <- data.frame(y = rnorm(100), a = rnorm(100), b = rnorm(100))
  fit <- ABT(y ~ a + b, data = d, n.trees = 20, cv.folds = 3,
             monitor = FALSE, seed = 1)
  expect_s3_class(fit, "ABT")
  expect_true(fit$best.n.trees >= 1)
  p <- predict(fit, newdata = d)
  expect_length(p, 100)
})