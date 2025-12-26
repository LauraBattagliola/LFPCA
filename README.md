# LFPCA: Localized Functional Principal Component Analysis

**LFPCA** performs Localized Functional Principal Component Analysis by detecting block-diagonal structure in the empirical covaraince function and computing FPCA within each block. 


For methodological details, see our preprint:

Battagliola, M.L. & Bauer, J.O. (2025). *Localized Functional Principal Component Analysis Based on Covariance Structure*. arXiv:2506.02836. [https://arxiv.org/abs/2506.02836](https://arxiv.org/abs/2506.02836)


## Usage

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("LauraBattagliola/LFPCA")
```

LFPCA depends on:
- [bdsvd](https://CRAN.R-project.org/package=bdsvd) - for block-diagonal structure detection
- [refund](https://CRAN.R-project.org/package=refund) - for functional PCA via FACE

Example using the pinch force dataset from the `fda` package:

```r
library(LFPCA)
library(fda)

# Load example data
data(pinch)
X <- t(pinch)  # 301 time points × 24 curves

# Run LFPCA
result <- LFPCA(
  data = X,
  Tmax = ncol(X),
  range = c(1, ncol(X)),
  Nsel = 3,
  smooth.data = TRUE,
  window_type = "cos"
)

# Access results
phi_df <- result$phi.df          # Eigenfunctions data frame
blocks <- result$blocks          # Detected block indices
eigenvalues <- result$eigenvalues # Top eigenvalues
subprocesses <- result$subprocess_numbers  # Block assignments
```

See the `inst/examples/pinch_example.R` file for a complete worked example.


Key parameters for `LFPCA()`:

- `data`: Matrix of functional observations (rows) × grid points (columns)
- `Tmax`: Number of grid points
- `range`: Time domain range, e.g., `c(0, 1)`
- `Nsel`: Number of principal components to select
- `window_type`: `"cos"` or `"epan"` for tapering style
- `delta`: Width of transition region for window functions
- `epsilon`: Minimum block size threshold (proportion of total points)
- `anp`: Which regularization function should be used for the HBIC (from bdsvd)

## License
MIT

