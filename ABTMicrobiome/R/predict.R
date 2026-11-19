#' Predict Method for ABT Objects
#'
#' @param object An \code{"ABT"} object.
#' @param newdata A data frame of new observations.
#' @param n.trees Number of trees (defaults to \code{object$best.n.trees}).
#' @param type Either \code{"link"} or \code{"response"}.
#' @param se Logical; return standard errors across component trees.
#' @param ... Unused.
#'
#' @return A vector (or matrix for multiclass) of predictions, or a list
#'   with \code{preds} and \code{preds.se} if \code{se = TRUE}.
#' @export
predict.ABT <- function(object, newdata, n.trees = NULL,
                        type = c("link", "response"), se = FALSE, ...) {
  type <- match.arg(type)

  if (is.null(n.trees)) n.trees <- object$best.n.trees
  if (length(n.trees) != 1) stop("n.trees must be a single integer.")

  vars <- object$var.names
  missing <- setdiff(vars, names(newdata))
  if (length(missing) > 0) {
    stop("newdata is missing required variables: ",
         paste(missing, collapse = ", "))
  }
  newdata <- newdata[, vars, drop = FALSE]

  max_trees <- length(object$component.trees[[1]]$trees)
  if (n.trees > max_trees) {
    warning(sprintf("n.trees = %d exceeds component trees (%d); using %d.",
                    n.trees, max_trees, max_trees))
    n.trees <- max_trees
  }
  if (n.trees <= 0) stop("n.trees must be a positive integer.")

  dist <- object$component.trees[[1]]$distribution$name
  is_multinomial <- dist == "multinomial"

  per_tree <- lapply(object$component.trees, function(obj) {
    gbm::predict.gbm(obj, newdata = newdata, n.trees = n.trees,
                     type = "link")
  })

  if (!is_multinomial) {
    mat <- do.call(cbind, per_tree)
    fit <- rowMeans(mat, na.rm = TRUE)

    if (se) {
      se_val <- apply(mat, 1, stats::sd, na.rm = TRUE) / sqrt(ncol(mat))
      return(.convert_response(list(preds = fit, preds.se = se_val),
                               dist, type))
    }
    return(.convert_response(fit, dist, type))
  }

  cls <- object$component.trees[[1]]$classes
  k <- length(cls)
  arr <- lapply(per_tree, function(v) {
    matrix(v, nrow = nrow(newdata), ncol = k, dimnames = list(NULL, cls))
  })
  fit <- Reduce(`+`, arr) / length(arr)
  colnames(fit) <- cls

  if (se) {
    se_arr <- array(NA_real_, dim = c(nrow(newdata), k),
                    dimnames = list(NULL, cls))
    for (j in seq_len(k)) {
      mj <- sapply(arr, function(m) m[, j])
      if (!is.matrix(mj)) mj <- matrix(mj, nrow = nrow(newdata))
      se_arr[, j] <- apply(mj, 1, stats::sd, na.rm = TRUE) / sqrt(ncol(mj))
    }
    return(.convert_response(list(preds = fit, preds.se = se_arr), dist, type))
  }
  .convert_response(fit, dist, type)
}