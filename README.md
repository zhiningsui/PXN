# PXN <img src="https://img.shields.io/badge/R-Package-276DC3.svg?logo=R&logoColor=white" alt="R"/>

PXN provides cross-platform normalization and prediction of gene expression using a probabilistic PCA (PPCA) model with platform-specific affine transforms. Train on paired multi-platform data, optionally refine with gradient descent, make cross-platform predictions, and (optionally) integrate several pairwise models into a grand model for indirect paths (e.g., A to C via B).

  - Paired multi-platform training (`InitEst`)
  - Optional refinement (`GDfun`)
  - Cross-platform prediction (`Prediction`)
  - Multi-model integration (`ModIntegrate`)
  - Built-in toy data for quick starts: `sim_expres`, `sim_covars`

## Installation Guide

This package is currently installed from source, typically from GitHub or from
a local clone of this repository. The guide below lists all required and
optional software dependencies.

### 1. Install R

PXN requires:

- R version 4.1.0 or newer
- A standard R package library that you can write to

Check your R version with:

```r
R.version.string
```

If your R version is older than 4.1.0, install a newer R release from
<https://cran.r-project.org/>.

### 2. Install System Build Tools

Some dependencies may be installed from source, so your computer may need
standard compilation tools.

- macOS: install Xcode Command Line Tools with `xcode-select --install`
- Windows: install the Rtools version that matches your R version from
  <https://cran.r-project.org/bin/windows/Rtools/>
- Linux: install your distribution's R development and build packages, such as
  `r-base-dev`, `build-essential`, `gfortran`, and the development libraries
  required by your local R package setup

These tools are not PXN runtime dependencies, but they are often needed during
installation.

### 3. Install Required R Package Dependencies

PXN depends on packages from both CRAN and Bioconductor. Run this once before
installing PXN:

```r
install.packages(c("Rfast", "remotes", "BiocManager"))
BiocManager::install("limma")
```

Required runtime dependencies:

| Dependency | Source | Why it is needed |
| --- | --- | --- |
| R >= 4.1.0 | CRAN/R Project | Base R runtime required by PXN |
| Rfast | CRAN | Fast row, column, matrix, and SVD helper operations |
| limma | Bioconductor | Empirical Bayes linear modeling utilities |
| stats | Base R | Statistical functions such as p-values and adjustment |
| graphics | Base R | Base plotting support |
| grDevices | Base R | Graphics device support |
| parallel | Recommended R package | Optional multicore cross-validation support |

PXN does not require Python, Java, a database, or external web services at
runtime.

### 4. Install PXN from GitHub

For a regular installation:

```r
remotes::install_github("zhiningsui/PXN", dependencies = TRUE)
```

To install without optional suggested packages:

```r
remotes::install_github("zhiningsui/PXN", dependencies = c("Depends", "Imports"))
```

To build vignettes during installation, install the optional vignette packages
first, make sure Pandoc is available, and then enable vignette building:

```r
install.packages(c("knitr", "rmarkdown", "ggplot2", "kableExtra"))
remotes::install_github(
  "zhiningsui/PXN",
  dependencies = TRUE,
  build_vignettes = TRUE
)
```

### 5. Install PXN from a Local Clone

If you already cloned this repository, open R in the repository root and run:

```r
remotes::install_local(".", dependencies = TRUE)
```

For a minimal local install with only runtime dependencies:

```r
remotes::install_local(".", dependencies = c("Depends", "Imports"))
```

If you prefer command-line R tooling, this also works from the repository root:

```sh
R CMD INSTALL .
```

### 6. Verify the Installation

After installation, start a fresh R session and run:

```r
library(PXN)
packageVersion("PXN")

data("sim_expres")
data("sim_covars")

stopifnot(
  is.list(sim_expres),
  is.list(sim_covars),
  ncol(sim_expres$AB) == ncol(sim_covars$AB)
)
```

You can also run a tiny model initialization check:

```r
genes <- 1:10
m <- nrow(sim_expres$AB) / 2
Y <- sim_expres$AB[c(genes, m + genes), 1:8]
X <- sim_covars$AB[, 1:8]
fit <- InitEst(X = X, Y = Y, K = 2, L = 1)
names(fit)
```

## Test Dataset and Example Run

The repository includes a small simulated test dataset for trying PXN end to
end:

- Expression data: [`data/sim_expres.rda`](data/sim_expres.rda)
- Covariate data: [`data/sim_covars.rda`](data/sim_covars.rda)

The same data are installed with the package, so you usually do not need to
download the files manually. Load them in R with `data("sim_expres")` and
`data("sim_covars")`.

Input:

- `sim_expres$AB`: stacked expression matrix for platforms A and B. Rows are
  genes/features, columns are samples, with A rows followed by B rows.
- `sim_covars$AB`: covariate matrix for the same AB samples. Columns must match
  the expression sample columns.
- `K = 2`: number of platforms in `sim_expres$AB`.
- `L = 1`: latent dimension used in this small test run.

Output:

- `fit`: trained PXN model parameters returned by `InitEst()`.
- `YB_pred`: predicted platform-B expression for held-out samples measured on
  platform A.
- `mse`: mean squared error comparing `YB_pred` with the held-out true platform-B
  expression.

Run the method on the test dataset:

```r
library(PXN)

data("sim_expres")
data("sim_covars")

# Use 100 genes and 20 samples so the example runs quickly.
m_total <- nrow(sim_expres$AB) / 2
genes <- 1:100
sample_ids <- 1:20

Y_AB <- sim_expres$AB[c(genes, m_total + genes), sample_ids, drop = FALSE]
X_AB <- sim_covars$AB[, sample_ids, drop = FALSE]

# Hold out the last 4 samples for A-to-B prediction.
m <- length(genes)
test_idx <- 17:20
train_idx <- setdiff(seq_along(sample_ids), test_idx)

Y_train <- Y_AB[, train_idx, drop = FALSE]
X_train <- X_AB[, train_idx, drop = FALSE]

YA_test <- Y_AB[1:m, test_idx, drop = FALSE]
YB_true <- Y_AB[(m + 1):(2 * m), test_idx, drop = FALSE]
X_test <- X_AB[, test_idx, drop = FALSE]

fit <- InitEst(X = X_train, Y = Y_train, K = 2, L = 1)

# Align held-out data to genes retained by InitEst().
retained_genes <- rownames(fit$b0)
YA_test <- YA_test[retained_genes, , drop = FALSE]
YB_true <- YB_true[retained_genes, , drop = FALSE]

YB_pred <- Prediction(
  Ysource = YA_test,
  X = X_test,
  trained.model = fit,
  k.source = 1,
  k.target = 2
)

mse <- mean((YB_pred - YB_true)^2)

dim(YB_pred)
mse
```

For this subset, `YB_pred` is a numeric matrix with retained genes as rows and
the 4 held-out samples as columns. The scalar `mse` summarizes prediction error
on the same held-out samples.

### Optional Dependencies

These packages are not required to use PXN, but they support documentation,
examples, testing, and development workflows.

| Dependency | Source | Used for |
| --- | --- | --- |
| knitr | CRAN | Building vignettes |
| rmarkdown | CRAN | Rendering vignette documents |
| ggplot2 | CRAN | Vignette plots and examples |
| testthat >= 3.0.0 | CRAN | Unit tests |
| covr | CRAN | Test coverage reports |
| kableExtra | CRAN | Rich vignette tables |
| remotes | CRAN | Installing PXN from GitHub or local source |
| devtools | CRAN | Alternative development workflow |
| BiocManager | CRAN | Installing Bioconductor packages such as `limma` |
| Pandoc | System tool | Rendering R Markdown vignettes |

Install all optional development dependencies with:

```r
install.packages(c(
  "knitr",
  "rmarkdown",
  "ggplot2",
  "testthat",
  "covr",
  "kableExtra",
  "remotes",
  "devtools",
  "BiocManager"
))
BiocManager::install("limma")
```

Pandoc is required only when rendering R Markdown vignettes. RStudio includes
Pandoc, or you can install it separately from <https://pandoc.org/installing.html>.

### Troubleshooting

If `limma` is not found, install it with Bioconductor:

```r
install.packages("BiocManager")
BiocManager::install("limma")
```

If `Rfast` fails to compile, confirm that your system build tools are installed
and that R can compile packages from source. On macOS and Windows, this usually
means installing Xcode Command Line Tools or Rtools.

If vignette building fails with a Pandoc error, install Pandoc or install PXN
without vignette rendering:

```r
remotes::install_github("zhiningsui/PXN", dependencies = TRUE, build_vignettes = FALSE)
```

If R reports that your package library is not writable, choose a user library:

```r
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
.libPaths(Sys.getenv("R_LIBS_USER"))
```

If an installation is interrupted by a stale lock directory, restart R and
remove the lock directory shown in the error message before trying again.

### Development Checks

Developers can run:

```r
devtools::document()
devtools::check()
testthat::test_dir("tests")
```

The repository currently ships package documentation and a vignette under
`vignettes/introduction.Rmd`.
