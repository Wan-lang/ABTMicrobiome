#' @export
print.ABT <- function(x, ...) {
  cat("Aggregated Boosted Trees (ABT)\n")
  cat("  component trees :", length(x$component.trees), "\n")
  cat("  trees per comp. :", x$best.n.trees, "\n")
  cat("  distribution    :", x$component.trees[[1]]$distribution$name, "\n")
  cat("  CV error        :", signif(x$best.cv.error, 5), "\n")
  invisible(x)
}

#' @export
summary.ABT <- function(object, ...) {
  print(object$var.importance)
  invisible(object$var.importance)
}

#' @export
plot.ABT <- function(x, type = c("importance", "cv"), ...) {
  type <- match.arg(type)
  if (type == "importance") plot_importance(x, ...)
  else plot_cv_error(x, ...)
}