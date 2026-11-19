.aggregate_importance <- function(component.trees, best) {
  vi_list <- lapply(component.trees, function(obj) {
    s <- gbm::summary.gbm(obj, plotit = FALSE, n.trees = best)
    stats::setNames(s$rel.inf, s$var)
  })
  all_vars <- unique(unlist(lapply(vi_list, names)))
  vi_mat <- sapply(vi_list, function(v) {
    out <- stats::setNames(rep(0, length(all_vars)), all_vars)
    out[names(v)] <- v
    out
  })
  df <- data.frame(var = all_vars,
                   rel.inf = rowMeans(vi_mat, na.rm = TRUE),
                   stringsAsFactors = FALSE)
  df <- df[order(df$rel.inf, decreasing = TRUE), ]
  rownames(df) <- NULL
  df
}

.convert_response <- function(x, dist, type) {
  if (type == "link") return(x)

  tf <- function(p) {
    if (dist %in% c("bernoulli", "pairwise")) return(1 / (1 + exp(-p)))
    if (dist == "poisson")  return(exp(p))
    if (dist == "adaboost") return(1 / (1 + exp(-2 * p)))
    if (dist == "coxph")    return(exp(p))
    if (dist == "laplace")  return(p)
    if (dist == "tdist")    return(p)
    p
  }

  if (dist == "multinomial") {
    softmax <- function(m) {
      e <- exp(m)
      e / rowSums(e)
    }
    if (is.list(x)) {
      x$preds <- softmax(x$preds)
      return(x)
    }
    return(softmax(x))
  }

  if (is.list(x)) {
    x$preds <- tf(x$preds)
    return(x)
  }
  tf(x)
}

.abt_theme <- function(base_size = 14) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold",
                                            size = base_size + 4),
      axis.title.x  = ggplot2::element_text(face = "bold", size = base_size + 1),
      axis.title.y  = ggplot2::element_text(face = "bold", size = base_size + 1),
      axis.text     = ggplot2::element_text(size = base_size - 3),
      panel.border  = ggplot2::element_rect(color = "black", fill = NA,
                                            linewidth = 0.6),
      panel.grid.major = ggplot2::element_line(color = "grey90",
                                               linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none"
    )
}

utils::globalVariables(c(
  "x", "y", "yhat", "pred",
  "trees", "cv_error",
  "rel.inf", "var",
  "Var1", "Var2", "Freq"
))