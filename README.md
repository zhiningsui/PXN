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

