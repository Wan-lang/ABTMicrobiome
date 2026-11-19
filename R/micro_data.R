#' Build Pairwise Microbiome-Environment Data
#'
#' Compute pairwise sample dissimilarities (Bray-Curtis, Jaccard, weighted
#' or unweighted UniFrac) and corresponding pairwise differences in
#' environmental variables. The resulting data frame is intended to be
#' used as `data` in \code{\link{ABT}}.
#'
#' @param env_data A data frame of environmental variables with samples
#'   in rows and variables in columns. Row names must be sample IDs.
#' @param otu_data A data frame of OTU/ASV abundances. Samples in columns
#'   when \code{margin = 2} (default) or in rows when \code{margin = 1}.
#' @param tree Optional path to a Newick tree file. Required when
#'   \code{datatype = "tre/nwk"}.
#' @param datatype Either \code{"otu/asv"} or \code{"tre/nwk"}.
#' @param method Distance method. For \code{"otu/asv"} use
#'   \code{"bray"} or \code{"jaccard"}; for \code{"tre/nwk"} use
#'   \code{"weighted unifrac"} or \code{"unweighted unifrac"}.
#' @param include_pairs Logical; include sample pair columns in the output.
#' @param transform Either \code{"abs"} (absolute differences) or
#'   \code{"signed"} (signed differences).
#' @param margin Margin of \code{otu_data} holding samples: 2 (columns,
#'   default) or 1 (rows).
#'
#' @return A data frame with one row per sample pair and columns
#'   \code{env_diff} (one per environmental variable) and
#'   \code{otu_value} (the dissimilarity). Attributes store metadata:
#'   \code{datatype}, \code{method}, \code{n_samples}, \code{pair_idx},
#'   \code{sample_names}.
#'
#' @examples
#' \dontrun{
#' otu  <- read.delim(system.file("extdata", "otutab.txt", package = "ABTMicrobiome"),
#'                    row.names = 1, check.names = FALSE)
#' env  <- read.delim(system.file("extdata", "env.txt", package = "ABTMicrobiome"),
#'                    row.names = 1, check.names = FALSE)
#' tree <- system.file("extdata", "otus.tre", package = "ABTMicrobiome")
#' md <- micro_data(env_data = env, otu_data = otu, tree = tree,
#'                  datatype = "tre/nwk", method = "weighted unifrac")
#' }
#' @export
micro_data <- function(env_data = NULL, otu_data = NULL, tree = NULL,
                       datatype = c("otu/asv", "tre/nwk"), method,
                       include_pairs = TRUE, transform = c("abs", "signed"),
                       margin = 2) {
  datatype  <- match.arg(datatype)
  transform <- match.arg(transform)

  if (is.null(env_data)) stop("env_data must be provided")
  if (is.null(otu_data)) stop("otu_data must be provided")
  env_data <- as.data.frame(env_data)
  otu_data <- as.data.frame(otu_data)

  env_data[] <- lapply(env_data, function(x) {
    if (is.factor(x)) x <- as.character(x)
    if (is.character(x)) as.numeric(trimws(x)) else x
  })
  if (anyNA(env_data)) stop("env_data has NA after numeric conversion")

  if (!is.character(method) || length(method) != 1)
    stop("method must be a single character string")
  method <- match.arg(method, c("bray", "jaccard",
                                "weighted unifrac", "unweighted unifrac"))

  if (nrow(env_data) != if (margin == 2) ncol(otu_data) else nrow(otu_data))
    stop("Sample count mismatch between env_data and otu_data.")

  if (!is.null(rownames(env_data)) && !is.null(colnames(otu_data))) {
    sample_names_otu <- if (margin == 2) colnames(otu_data) else rownames(otu_data)
    if (!setequal(rownames(env_data), sample_names_otu))
      stop("Sample names in env_data and otu_data do not match")
    if (margin == 2) otu_data <- otu_data[, rownames(env_data), drop = FALSE]
    else             otu_data <- otu_data[rownames(env_data), , drop = FALSE]
  }

  otu_t <- if (margin == 2) t(otu_data) else as.matrix(otu_data)
  otu_t <- otu_t / rowSums(otu_t)

  if (datatype == "otu/asv") {
    if (!method %in% c("bray", "jaccard"))
      stop("For 'otu/asv', method must be 'bray' or 'jaccard'")
    if (!is.null(tree)) warning("Tree not used in 'otu/asv' mode")
    if (!requireNamespace("vegan", quietly = TRUE))
      stop("Package 'vegan' is required for 'otu/asv' distance.")
    D <- as.matrix(vegan::vegdist(otu_t, method = method))
  } else {
    if (!method %in% c("weighted unifrac", "unweighted unifrac"))
      stop("For 'tre/nwk', method must be 'weighted unifrac' or 'unweighted unifrac'")
    if (is.null(tree) || !file.exists(tree))
      stop("tree file path must be provided and exist")
    if (!requireNamespace("ape", quietly = TRUE) ||
        !requireNamespace("GUniFrac", quietly = TRUE))
      stop("Packages 'ape' and 'GUniFrac' are required for 'tre/nwk'.")
    tree_data <- ape::read.tree(tree)
    common <- intersect(colnames(otu_t), tree_data$tip.label)
    if (length(common) == 0) stop("No matching OTU names between otu_data and tree")
    tree_data <- ape::keep.tip(tree_data, common)
    otu_t <- otu_t[, tree_data$tip.label, drop = FALSE]
    uf <- GUniFrac::GUniFrac(otu_t, tree_data)$unifracs
    key <- if (method == "weighted unifrac") "d_1" else "d_UW"
    if (!key %in% dimnames(uf)[[3]])
      stop(sprintf("GUniFrac output missing '%s'; available: %s",
                   key, paste(dimnames(uf)[[3]], collapse = ", ")))
    D <- as.matrix(as.dist(uf[, , key]))
  }

  n <- nrow(D)
  if (n < 2) stop("At least 2 samples required")
  if (is.null(rownames(D))) rownames(D) <- colnames(D) <- paste0("S", seq_len(n))

  pair_idx <- as.data.frame(t(utils::combn(n, 2)))
  names(pair_idx) <- c("i", "j")

  env_mat <- as.matrix(env_data)
  diff_raw <- env_mat[pair_idx$i, , drop = FALSE] -
    env_mat[pair_idx$j, , drop = FALSE]
  if (transform == "abs") diff_raw <- abs(diff_raw)
  env_diff <- as.data.frame(diff_raw)
  names(env_diff) <- colnames(env_data)

  out <- data.frame(env_diff,
                    otu_value = D[cbind(pair_idx$i, pair_idx$j)],
                    stringsAsFactors = FALSE)
  if (include_pairs) {
    out <- cbind(
      Sample_i = rownames(D)[pair_idx$i],
      Sample_j = rownames(D)[pair_idx$j],
      out,
      stringsAsFactors = FALSE
    )
  }

  attr(out, "datatype")    <- datatype
  attr(out, "method")      <- method
  attr(out, "n_samples")   <- n
  attr(out, "pair_idx")    <- pair_idx
  attr(out, "sample_names")<- rownames(D)
  out
}