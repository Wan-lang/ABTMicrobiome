test_that("bundled example data are readable", {
  otu_path  <- system.file("extdata", "otutab.txt", package = "ABTMicrobiome")
  env_path  <- system.file("extdata", "env.txt", package = "ABTMicrobiome")
  tree_path <- system.file("extdata", "otus.tre", package = "ABTMicrobiome")
  bar_path  <- system.file("extdata", "barramundi.csv", package = "ABTMicrobiome")

  skip_if(otu_path == "", "otutab.txt not installed")
  skip_if(env_path == "", "env.txt not installed")
  skip_if(bar_path == "", "barramundi.csv not installed")

  otu <- read.delim(otu_path, row.names = 1, check.names = FALSE)
  env <- read.delim(env_path, row.names = 1, check.names = FALSE)
  bar <- read.csv(bar_path, check.names = FALSE)

  expect_true(ncol(otu) == 18)
  expect_true(nrow(env) == 18)
  expect_setequal(rownames(env), colnames(otu))
  expect_true(all(c("Sr","Ba","Mn","Mg","Fe","Zn","Na","K","P","S","Fresh") %in% names(bar)))
})

test_that("micro_data works on bundled env/otu (bray)", {
  otu_path <- system.file("extdata", "otutab.txt", package = "ABTMicrobiome")
  env_path <- system.file("extdata", "env.txt", package = "ABTMicrobiome")
  skip_if(otu_path == "" || env_path == "")

  otu <- read.delim(otu_path, row.names = 1, check.names = FALSE)
  env <- read.delim(env_path, row.names = 1, check.names = FALSE)

  md <- micro_data(env_data = env, otu_data = otu,
                   datatype = "otu/asv", method = "bray")
  expect_equal(nrow(md), choose(18, 2))
  expect_true("otu_value" %in% names(md))
  expect_true(all(c("pH","EC","Ca","Mg","K","Na","CO","HCO","Cl","SO") %in% names(md)))
})

test_that("micro_data works on bundled env/otu (weighted unifrac) if tree exists", {
  otu_path  <- system.file("extdata", "otutab.txt", package = "ABTMicrobiome")
  env_path  <- system.file("extdata", "env.txt", package = "ABTMicrobiome")
  tree_path <- system.file("extdata", "otus.tre", package = "ABTMicrobiome")
  skip_if(otu_path == "" || env_path == "" || tree_path == "")

  otu <- read.delim(otu_path, row.names = 1, check.names = FALSE)
  env <- read.delim(env_path, row.names = 1, check.names = FALSE)
  md <- micro_data(env_data = env, otu_data = otu, tree = tree_path,
                   datatype = "tre/nwk", method = "weighted unifrac")
  expect_equal(nrow(md), choose(18, 2))
  expect_true(all(md$otu_value >= 0))
})

test_that("ABT fits gaussian model on pairwise distance data", {
  otu_path <- system.file("extdata", "otutab.txt", package = "ABTMicrobiome")
  env_path <- system.file("extdata", "env.txt", package = "ABTMicrobiome")
  skip_if(otu_path == "" || env_path == "")

  otu <- read.delim(otu_path, row.names = 1, check.names = FALSE)
  env <- read.delim(env_path, row.names = 1, check.names = FALSE)
  md  <- micro_data(env_data = env, otu_data = otu,
                    datatype = "otu/asv", method = "bray")

  fit <- ABT(otu_value ~ pH + EC + Ca + Mg + K + Na + CO + HCO + Cl + SO,
             data = md, distribution = "gaussian",
             n.trees = 30, cv.folds = 3, n.minobsinnode = 2,
             seed = 42, monitor = FALSE)
  expect_s3_class(fit, "ABT")
  expect_true(fit$best.n.trees >= 1)
  expect_true(all(fit$var.importance$rel.inf >= 0))
})

test_that("ABT fits bernoulli model on barramundi data", {
  bar_path <- system.file("extdata", "barramundi.csv", package = "ABTMicrobiome")
  skip_if(bar_path == "")
  bar <- read.csv(bar_path, check.names = FALSE)
  bar$Fresh <- as.factor(bar$Fresh)
  fit <- ABT(Fresh ~ Sr + Ba + Mn + Mg + Fe + Zn + Na + K + P + S,
             data = bar, distribution = "bernoulli",
             n.trees = 30, cv.folds = 3, seed = 1, monitor = FALSE)
  expect_s3_class(fit, "ABT")
  pr <- predict(fit, newdata = bar, type = "response")
  expect_true(all(pr >= 0 & pr <= 1))
})