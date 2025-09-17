# PXN <img src="https://img.shields.io/badge/R-Package-276DC3.svg?logo=R&logoColor=white" alt="R"/>

PXN provides cross-platform normalization and prediction of gene expression using a probabilistic PCA (PPCA) model with platform-specific affine transforms. Train on paired multi-platform data, optionally refine with gradient descent, make cross-platform predictions, and (optionally) integrate several pairwise models into a grand model for indirect paths (e.g., A→C via B).

  - Paired multi-platform training (`InitEst`)
  - Optional refinement (`GDfun`)
  - Cross-platform prediction (`Prediction`)
  - Multi-model integration (`ModIntegrate`)
  - Built-in toy data for quick starts: `sim_expres`, `sim_covars`
