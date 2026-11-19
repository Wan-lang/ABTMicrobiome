#' Aggregated Boosted Trees
#'
#' Fit \code{cv.folds} gradient boosting models with cross-validation and
#' aggregate their predictions and variable importances. This is a
#' lightweight wrapper around \pkg{gbm} designed for the microbiome
#' distance analysis workflow.
#'
#' @param formula A model formula.
#' @param data A data frame (typically output of \code{\link{micro_data}}).
#' @param distribution Character; one of \code{"gaussian"}, \code{"bernoulli"},
#'   \code{"poisson"}, \code{"adaboost"}, \code{"laplace"}, \code{"coxph"}.
#' @param weights Optional observation weights.
#' @param offset Optional offset.
#' @param var.monotone Optional monotonicity constraints.
#' @param n.trees Maximum number of trees per fold.
#' @param interaction.depth Tree depth.
#' @param n.minobsinnode Minimum observations per terminal node.
#' @param shrinkage Learning rate.
#' @param bag.fraction Bagging fraction.
#' @param keep.data Logical; keep data in component models.
#' @param verbose Logical; print progress.
#' @param stratify Logical; stratify CV folds by the outcome.
#' @param monitor Logical; print a model summary.
#' @param use Either \code{"best"} or \code{"all"} (reserved).
#' @param seed Random seed.
#' @param na.omit.y Drop rows with missing outcome.
#' @param cv.folds Number of CV folds (\eqn{\geq 2}).
#' @param train.fraction Reserved.
#' @param class.stratify.cv Reserved.
#'
#' @return An object of class \code{"ABT"} (a list) with components
#'   \code{component.trees}, \code{cv.predictions}, \code{cv.folds},
#'   \code{cv.best.trees}, \code{best.n.trees}, \code{best.cv.error},
#'   \code{cv.error}, \code{cv.error.matrix}, \code{var.importance}, and
#'   model hyper-parameters.
#'
#' @export

ABT <- function(formula = formula(data), data = list(),
                distribution = "gaussian",
                weights, offset = NULL, var.monotone = NULL, n.trees = 200,
                interaction.depth = 3, n.minobsinnode = NULL,
                shrinkage = 0.05, bag.fraction = 0.5, keep.data = FALSE,
                verbose = FALSE, stratify = TRUE, monitor = TRUE,
                use = c("best", "all")[1], seed = 0, na.omit.y = TRUE,
                cv.folds = 5, train.fraction = 1, class.stratify.cv = NULL)
{
  if (!is.numeric(cv.folds) || length(cv.folds) != 1 ||
      cv.folds <= 1 || cv.folds != as.integer(cv.folds)) {
    stop("cv.folds must be a single integer greater than 1.")
  }
  use <- match.arg(use, c("best", "all"))
  if (is.null(n.minobsinnode)) n.minobsinnode <- ceiling(log10(nrow(data)) + 1)

  m <- match.call(expand.dots = FALSE)
  m <- m[c(1, match(c("formula", "data", "weights", "offset"),
                    names(m), 0))]
  m[[1]] <- as.name("model.frame")
  m$na.action <- na.pass
  m <- eval(m, parent.frame())

  Terms <- attr(m, "terms")
  a <- attributes(Terms)
  y <- model.response(m)
  x <- model.frame(delete.response(Terms), data, na.action = na.pass)
  w <- model.extract(m, weights)
  if (is.null(w)) w <- rep(1, length(y))

  if (!is.null(model.extract(m, "offset"))) {
    stop("In ABT the offset term needs to be specified using the offset parameter.")
  }

  if (na.omit.y && any(nays <- is.na(y))) {
    y <- y[!nays]; x <- x[!nays, , drop = FALSE]; w <- w[!nays]
    if (!is.null(offset)) offset <- offset[!nays]
  }
  if (any(is.na(x))) {
    warning("Missing values found in predictors. Consider imputation.")
  }
  nrows <- nrow(x)
  if (length(y) != nrows) stop("Internal error: x/y length mismatch after NA handling.")

  var.names <- a$term.labels
  response.name <- as.character(attr(Terms, "variables"))[[2]]

  distributions <- c("bernoulli", "gaussian", "poisson", "adaboost", "laplace", "coxph")
  distribution <- distributions[pmatch(distribution, distributions)]
  if (is.na(distribution)) {
    stop("Unsupported distribution. Choose from: ",
         paste(distributions, collapse = ", "))
  }
  invisible(gbm::getStratify(class.stratify.cv, d = list(name = distribution)))

  if (seed > 0) set.seed(seed)

  if (!stratify) {
    nsamps <- sample(rep(1:cv.folds, length = nrows))
  } else {
    nsamps <- integer(nrows)
    if (is.factor(y) || is.character(y)) {
      grp <- as.factor(y)
    } else {
      nbin <- min(cv.folds, max(2, floor(nrows / 10)))
      brks <- unique(quantile(y, probs = seq(0, 1, length.out = nbin + 1)))
      grp  <- cut(y, breaks = brks, include.lowest = TRUE, labels = FALSE)
    }
    for (g in unique(grp)) {
      idx <- which(grp == g)
      nsamps[idx] <- sample(rep(1:cv.folds, length.out = length(idx)))
    }
  }

  valid.err.matrix <- matrix(NA_real_, nrow = n.trees, ncol = cv.folds)
  cv.best.trees   <- integer(cv.folds)
  component.trees <- vector("list", cv.folds)
  preds <- rep(NA_real_, nrows)
  fold.n <- integer(cv.folds)

  for (i in seq_len(cv.folds)) {
    trainees     <- which(nsamps != i)
    non.trainees <- which(nsamps == i)
    nTrain <- length(trainees)
    fold.n[i] <- length(non.trainees)

    ord <- c(trainees, non.trainees)
    temp.x <- x[ord, , drop = FALSE]
    temp.y <- y[ord]
    temp.w <- w[ord]
    temp.offset <- if (!is.null(offset)) offset[ord] else NULL

    gbm.obj <- gbm::gbm.fit(
      x = temp.x, y = temp.y, offset = temp.offset,
      distribution = distribution, w = temp.w,
      var.monotone = var.monotone,
      n.trees = n.trees, interaction.depth = interaction.depth,
      n.minobsinnode = n.minobsinnode, shrinkage = shrinkage,
      bag.fraction = bag.fraction, nTrain = nTrain,
      keep.data = keep.data, verbose = verbose,
      var.names = var.names, response.name = response.name
    )

    ve <- gbm.obj$valid.error
    if (length(ve) < n.trees) ve <- c(ve, rep(NA_real_, n.trees - length(ve)))
    valid.err.matrix[, i] <- ve[seq_len(n.trees)]

    bi <- which.min(ve)
    cv.best.trees[i] <- bi
    preds[non.trainees] <- gbm::predict.gbm(
      gbm.obj, newdata = x[non.trainees, , drop = FALSE], n.trees = bi
    )
    component.trees[[i]] <- gbm.obj
  }

  valid.err <- as.numeric(valid.err.matrix %*% fold.n) / sum(fold.n)
  if (all(is.na(valid.err))) stop("All CV errors are NA; cannot select tree count.")

  best <- which.min(valid.err)
  best.error <- valid.err[best]
  if (best == n.trees && monitor)
    warning("Optimal trees == n.trees; consider increasing n.trees.")
  if (best == 1 && monitor)
    warning("Optimal trees == 1; model may be under-fitted.")

  for (i in seq_len(cv.folds)) component.trees[[i]]$n.trees <- best

  var.importance <- .aggregate_importance(component.trees, best)

  if (monitor) {
    cat("\n========================================\n")
    cat("ABT Model Summary\n========================================\n")
    cat("N trees (max)      =", n.trees, "\n")
    cat("Best N trees       =", best, "\n")
    cat("Depth              =", interaction.depth, "\n")
    cat("Minimum node size  =", n.minobsinnode, "\n")
    cat("Distribution       =", distribution, "\n")
    cat("Shrinkage          =", signif(shrinkage, 3), "\n")
    cat("Bag Fraction       =", bag.fraction, "\n")
    cat("CV Folds           =", cv.folds, "\n")
    cat("Best CV error      =", signif(best.error, 5), "\n")
    cat("Per-fold best trees:", paste(cv.best.trees, collapse = ", "), "\n")
    cat("========================================\n\n")
  }

  result <- list(
    component.trees = component.trees,
    cv.predictions  = preds,
    cv.folds        = nsamps,
    cv.best.trees   = cv.best.trees,
    best.n.trees    = best,
    best.cv.error   = best.error,
    cv.error        = valid.err,
    cv.error.matrix = valid.err.matrix,
    shrinkage       = shrinkage,
    interaction.depth = interaction.depth,
    n.minobsinnode  = n.minobsinnode,
    bag.fraction    = bag.fraction,
    distribution    = distribution,
    var.importance  = var.importance,
    use             = use,
    formula         = formula,
    data            = data,
    var.names       = var.names,
    response.name   = response.name
  )
  class(result) <- "ABT"
  result
}