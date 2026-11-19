test_that("micro_data validates inputs", {
  expect_error(micro_data(otu_data = matrix(1:6, 2)),
               "env_data must be provided")
  expect_error(micro_data(env_data = data.frame(a = 1:2)),
               "otu_data must be provided")
})

test_that("micro_data produces pair table", {
  set.seed(1)
  env <- data.frame(x = rnorm(5), y = rnorm(5), row.names = paste0("S", 1:5))
  otu <- matrix(rpois(40, 5), nrow = 8, ncol = 5,
                dimnames = list(paste0("OTU", 1:8), paste0("S", 1:5)))
  md <- micro_data(env_data = env, otu_data = otu,
                   datatype = "otu/asv", method = "bray")
  expect_equal(nrow(md), choose(5, 2))
  expect_true("otu_value" %in% names(md))
  expect_equal(attr(md, "n_samples"), 5)
})