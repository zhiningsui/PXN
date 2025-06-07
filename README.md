# PPCAXNORM: Probabilistic Principal Component Analysis for Cross-Platform Normalization

## Overview
`PPCAXNORM` is an R package that provides a powerful framework for normalizing and integrating high-dimensional data from multiple platforms, such as different microarray or RNA-sequencing experiments. 
It implements a probabilistic principal component analysis (PPCA) based model that can account for platform-specific effects, adjust for covariates, and model shared latent structures in the data.

The key features of `PPCAXNORM` include:
- Cross-platform normalization: Normalize data from a source platform to a target platform.
- Covariate adjustment: Incorporate and adjust for known covariates, such as age or treatment group.
- Model training and parameter estimation: Functions for initializing model parameters and refining them using gradient descent.
- Model integration: Combine models trained on different datasets into a single, integrated model.
- Cross-validation: Tools for selecting the optimal number of latent factors.

## Installation
You can install the development version of PPCAXNORM from GitHub with:
```R
# install.packages("devtools")
devtools::install_github("yourusername/PPCAXNORM")
```

## Quick Example
Here is a simple example of how to use `PPCAXNORM` to normalize data between two platforms.
```R
library(PPCAXNORM)

# 1. Generate some example data
set.seed(123)
m <- 100 # number of genes
n <- 50  # number of samples
K <- 2   # number of platforms

# Create a stacked data matrix Y (platform 1 and platform 2)
Y <- rbind(
  matrix(rnorm(m * n, mean = 2), nrow = m),
  matrix(rnorm(m * n, mean = 5), nrow = m)
)

# Create a covariate matrix X
X <- matrix(rnorm(n * 2), nrow = 2)
rownames(X) <- c("age", "treatment")

# 2. Get initial parameter estimates
# We'll use L=3 latent factors
init_params <- InitEst(X = X, Y = Y, K = K, L = 3)

# 3. Refine the model using Gradient Descent
trained_model <- GDfun(X = X, Y = Y, K = K, params = init_params)

# 4. Predict data from platform 1 to platform 2
# Let's use the first 10 samples as our source data
Y_source <- Y[1:m, 1:10]
X_source <- X[, 1:10]

Y_predicted <- Prediction(
  Ysource = Y_source,
  X = X_source,
  trained.model = trained_model,
  k.source = 1,
  k.target = 2
)

# View the dimensions of the predicted data
dim(Y_predicted)
```
## Citation
If you use `PPCAXNORM` in your research, please cite it as follows:
