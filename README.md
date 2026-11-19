# ABTMicrobiome

**Aggregated Boosted Trees for Microbiome–Environment Association**

`ABTMicrobiome` combines **Aggregated Boosted Trees (ABT)** with **microbial community distance analysis** to answer a central question:

> **To what extent, and in what nonlinear ways, do environmental factors drive microbial community dissimilarity?**

It models pairwise "environmental difference × community distance" data with cross-validated gradient boosted trees and provides diagnostics for variable importance, CV error, partial dependence, and interaction strength.

---

## Installation

```r
# From GitHub
remotes::install_github("Wan-lang/ABTMicrobiome")
```

**Dependencies:** `ggplot2`, `gbm`, `vegan`, `ape`, `GUniFrac`, `rlang`.

---

## Quick Start

### 1. Load example data

```r
library(ABTMicrobiome)

otu  <- read.delim(system.file("extdata", "otutab.txt", package = "ABTMicrobiome"),
                   row.names = 1, check.names = FALSE)
env  <- read.delim(system.file("extdata", "env.txt", package = "ABTMicrobiome"),
                   row.names = 1, check.names = FALSE)
tree <- system.file("extdata", "otus.tre", package = "ABTMicrobiome")
```

### 2. Build pairwise data

`micro_data()` computes pairwise community distances and environmental differences, returning a pairwise data frame ready for ABT:

```r
# Bray-Curtis (no tree required)
md_bray <- micro_data(env_data = env, otu_data = otu,
                      datatype = "otu/asv", method = "bray")

# Weighted UniFrac (requires a Newick tree)
md_uf <- micro_data(env_data = env, otu_data = otu, tree = tree,
                    datatype = "tre/nwk", method = "weighted unifrac")
```

Output columns: `Sample_i`, `Sample_j`, environmental differences (absolute values when `transform = "abs"`), and `otu_value` (community distance).

### 3. Fit ABT

```r
fit <- ABT(
  otu_value ~ pH + EC + Ca + Mg + K + Na + CO + HCO + Cl + SO,
  data           = md_uf,
  distribution   = "gaussian",
  n.trees        = 500,
  cv.folds       = 5,
  n.minobsinnode = 2,
  seed           = 42
)
```

### 4. Inspect and visualize

```r
summary(fit)                                           # variable importance table
plot_importance(fit)                                   # importance bar chart
plot_cv_error(fit)                                     # CV error curve (check best.n.trees)
plot_pdp(fit, vars = head(fit$var.importance$var, 3)) # partial dependence
plot_interaction(fit)                                  # pairwise interaction heatmap
```

### 5. Binary classification example (barramundi)

```r
bar <- read.csv(system.file("extdata", "barramundi.csv", package = "ABTMicrobiome"),
                check.names = FALSE)
bar$Fresh <- as.factor(bar$Fresh)

fit_bin <- ABT(Fresh ~ Sr + Ba + Mn + Mg + Fe + Zn + Na + K + P + S,
               data = bar, distribution = "bernoulli",
               n.trees = 100, cv.folds = 5)

pr <- predict(fit_bin, newdata = bar, type = "response")
```

---

## Main Functions

| Function | Purpose |
| :--- | :--- |
| `micro_data()` | Compute pairwise community distances and environmental differences |
| `ABT()` | Cross-validated aggregated boosted trees (supports `gaussian`, `bernoulli`, etc.) |
| `predict(fit, newdata)` | Predict from an ABT object; supports `type = "link"/"response"` and `se = TRUE` |
| `plot_importance()` | Variable relative-importance bar chart |
| `plot_cv_error()` | CV error vs. number of trees, with best-tree marker |
| `plot_pdp()` / `plot_pdp2()` | 1D / 2D partial dependence plots |
| `plot_interaction()` | Pairwise interaction-strength heatmap |

---

## Data Format

* **`env_data` (environmental variables)**
  * Row names = sample IDs
  * Columns = numeric environmental variables (character columns are coerced to numeric)
  * Must be convertible to numeric, otherwise an error is raised
* **`otu_data` (OTU/ASV table)**
  * Column names = sample IDs (must match `rownames(env_data)`)
  * Rows = OTU/ASV
  * Automatically converted to relative abundance (each sample divided by its total)
* **`tree` (phylogenetic tree)**
  * Path to a Newick file
  * Tip labels must overlap with `rownames(otu_data)`; `micro_data()` keeps the intersection and prunes the tree accordingly

---

## Method Overview

1. **Pairwise construction:** For $n$ samples, build $\binom{n}{2}$ sample pairs.
2. **Response:** Community distance $d_{ij}$ per pair (Bray-Curtis or UniFrac).
3. **Predictors:** Environmental differences $|\Delta \text{env}_{ij}|$ per pair.
4. **Modeling:** Fit $d_{ij} = f(|\Delta \text{env}_{ij}|)$ with `gbm`, selecting the optimal number of trees via cross-validation.
5. **Aggregation:** Average across fold models for more robust variable importance and predictions.

---

## Citation

If you use this package in your research, please cite:

```bibtex
@software{ABTMicrobiome,
  author  = {Wang, Xiaoyu},
  title   = {ABTMicrobiome: Aggregated Boosted Trees for Microbiome-Environment Association},
  year    = {2025},
  url     = {https://github.com/your-username/ABTMicrobiome}
}
```

Please also cite the underlying methods:

* **`gbm`:** Ridgeway G. (2007). Generalized Boosted Models: A guide to the gbm package.
