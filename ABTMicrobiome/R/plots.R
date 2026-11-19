#' Plot Variable Importance
#'
#' @param object An \code{"ABT"} object.
#' @param top Number of top variables to show.
#' @param gradient Logical; use a gradient fill.
#' @param fill Solid fill colour when \code{gradient = FALSE}.
#' @return A \code{ggplot} object.
#' @export
plot_importance <- function(object, top = NULL, gradient = TRUE,
                            fill = "#2c7fb8") {
  df <- object$var.importance
  if (!is.null(top)) df <- utils::head(df, top)
  df$var <- factor(df$var, levels = rev(df$var))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = rel.inf, y = var))
  if (gradient) {
    p <- p + ggplot2::geom_col(ggplot2::aes(fill = rel.inf), width = 0.85) +
      ggplot2::scale_fill_gradientn(
        colours = c("#FF6347", "#FFA500", "#32CD32", "#20B2AA",
                    "#00BFFF", "#9370DB", "#DA70D6", "#FF69B4"),
        guide = "none"
      )
  } else {
    p <- p + ggplot2::geom_col(fill = fill, width = 0.85)
  }
  p +
    ggplot2::labs(x = "Relative Influence", y = "Variable",
                  title = "Relative Influence of Variables") +
    .abt_theme() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

#' Plot CV Error Curve
#' @param object An \code{"ABT"} object.
#' @param log_x Logical; log-scale the x-axis.
#' @return A \code{ggplot} object.
#' @export
plot_cv_error <- function(object, log_x = FALSE) {
  df <- data.frame(trees = seq_along(object$cv.error),
                   cv_error = object$cv.error)
  p <- ggplot2::ggplot(df, ggplot2::aes(trees, cv_error)) +
    ggplot2::geom_line(color = "#d95f02", linewidth = 1) +
    ggplot2::geom_vline(xintercept = object$best.n.trees,
                        linetype = "dashed", color = "grey40") +
    ggplot2::labs(x = "Number of trees", y = "CV error",
                  title = sprintf("CV Error (best = %d)", object$best.n.trees)) +
    .abt_theme()
  if (log_x) p <- p + ggplot2::scale_x_log10()
  p
}

#' Plot Partial Dependence
#' @param object An \code{"ABT"} object.
#' @param vars Character vector of variable names.
#' @param n.trees Number of trees.
#' @param ncol Number of facet columns.
#' @return A \code{ggplot} object.
#' @export
plot_pdp <- function(object, vars, n.trees = NULL, ncol = NULL) {
  if (is.null(n.trees)) n.trees <- object$best.n.trees
  dfs <- lapply(vars, function(v) {
    per_tree <- lapply(object$component.trees, function(obj) {
      g <- gbm::plot.gbm(obj, i.var = v, n.trees = n.trees,
                         return.grid = TRUE)
      names(g)[2] <- "yhat"
      g
    })
    xcol <- names(per_tree[[1]])[1]
    grid <- per_tree[[1]][[xcol]]
    yhat <- rowMeans(sapply(per_tree, function(g) g$yhat))
    data.frame(var = v, x = grid, yhat = yhat)
  })
  df <- do.call(rbind, dfs)
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = yhat)) +
    ggplot2::geom_line(color = "#1b9e77", linewidth = 1) +
    ggplot2::facet_wrap(~ var, scales = "free_x", ncol = ncol) +
    ggplot2::labs(x = NULL, y = "Predicted response",
                  title = "Partial Dependence (ABT)") +
    .abt_theme() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 12))
}

#' Plot 2D Partial Dependence
#' @param object An \code{"ABT"} object.
#' @param v1 First variable.
#' @param v2 Second variable.
#' @param n.trees Number of trees.
#' @return A \code{ggplot} object.
#' @export
plot_pdp2 <- function(object, v1, v2, n.trees = NULL) {
  if (is.null(n.trees)) n.trees <- object$best.n.trees
  per_tree <- lapply(object$component.trees, function(obj) {
    gbm::plot.gbm(obj, i.var = c(v1, v2), n.trees = n.trees,
                  return.grid = TRUE)
  })
  base <- per_tree[[1]]
  base[[ncol(base)]] <- rowMeans(sapply(per_tree, function(g) g[[ncol(g)]]))
  names(base)[1:2] <- c("x", "y")
  names(base)[3] <- "pred"
  ggplot2::ggplot(base, ggplot2::aes(x, y, fill = pred)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colours = c("#FF6347", "#FFA500", "#32CD32", "#20B2AA",
                  "#00BFFF", "#9370DB", "#DA70D6", "#FF69B4"),
      name = "Pred"
    ) +
    ggplot2::labs(x = v1, y = v2,
                  title = sprintf("2D Partial Dependence: %s x %s", v1, v2)) +
    .abt_theme() +
    ggplot2::theme(legend.position = "right")
}

#' Plot Pairwise Interaction Strength
#' @param object An \code{"ABT"} object.
#' @param vars Variables to include (defaults to top 8).
#' @param n.trees Number of trees.
#' @return A \code{ggplot} object.
#' @export
plot_interaction <- function(object, vars = NULL, n.trees = NULL) {
  if (is.null(n.trees)) n.trees <- object$best.n.trees
  if (is.null(vars)) vars <- utils::head(object$var.importance$var, 8)
  k <- length(vars)
  H <- matrix(NA_real_, k, k, dimnames = list(vars, vars))

  pdp2 <- function(v1, v2) {
    per_tree <- lapply(object$component.trees, function(obj) {
      gbm::plot.gbm(obj, i.var = c(v1, v2), n.trees = n.trees,
                    return.grid = TRUE)
    })
    base <- per_tree[[1]]
    base[[ncol(base)]] <- rowMeans(sapply(per_tree, function(g) g[[ncol(g)]]))
    base
  }

  total_var <- stats::var(object$cv.predictions, na.rm = TRUE)
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      if (i == j) { H[i, j] <- 1; next }
      if (j < i)  { H[i, j] <- H[j, i]; next }
      g2 <- pdp2(vars[i], vars[j])
      H[i, j] <- stats::var(g2[[ncol(g2)]]) / total_var
    }
  }
  df <- as.data.frame(as.table(H))
  ggplot2::ggplot(df, ggplot2::aes(Var1, Var2, fill = Freq)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "white", high = "#b2182b",
                                 name = "Interaction") +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Pairwise Interaction Strength") +
    .abt_theme() +
    ggplot2::theme(
      legend.position = "right",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.major = ggplot2::element_blank()
    )
}