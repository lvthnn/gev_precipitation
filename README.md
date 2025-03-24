# Bayesian Modelling of Extreme Precipitation Events

This repository contains the implementation of Bayesian modelling methods for extreme precipitation events using Generalized Extreme Value (GEV) distributions. The analysis focuses on 132 annual maximum precipitation measurements from a weather station in southern England.

## Overview

The repository implements several models (excluding the trend effects $\beta_j$):

- Stationary GEV model with fixed parameters
- Non-stationary models with evolving location parameters:
  - Linear time trend model
  - Piecewise linear trend models with various prior structures (IID, AR(1), random walk)

Model evaluation using WAIC and LOOIC information criteria demonstrates that non-stationary models better characterize the data, with the random walk piecewise model yielding the most favourable results.

## Repository Structure

- `data/`: Precipitation datasets
  - `dt_maxprecip.csv`: Maximum precipitation data
  - `precip_data.csv`: General precipitation data
- `pdf/`: Documentation
  - `report/report_full.pdf`: Comprehensive report with methodology and results
- `src/`: Implementation
  - `stan/`: Stan model implementations
    - `gev.stan`: Basic GEV distribution model
    - `hierarchical_ar1.stan`: Hierarchical model with AR(1) process
    - `hierarchical_iid.stan`: Hierarchical model with IID structure
    - `hierarchical_rw.stan`: Hierarchical model with random walk
    - `nullmodel.stan`: Baseline model for comparison
    - `timetrend.stan`: Linear time trend model
  - `beta_function_basis.R`: Beta function basis implementation
  - `model_fitting.R`: Model fitting functionality
  - `post_analysis.R`: Post-processing and analysis
  - `utils.R`: Utility functions

## Key Results

The random walk piecewise model with m=18 coefficients produced the most favourable results. The optimal model shows an upward trend in annual maximum precipitation for approximately 100 years, followed by a potential downward trend in the remaining years.

Parameter estimates for the optimal model:

| Parameter | Mean    | Median  | 95% Credible Interval |
|-----------|---------|---------|------------------------|
| $\mu_0$        | 20.7    | 20.7    | (18.0, 23.3)           |
| $\sigma$        | 6.56    | 6.54    | (5.72, 7.50)           |
| $\xi$   | 0.0593  | 0.0574  | (-0.0587, 0.187)       |

## Documentation

For a complete description of the methodology, theoretical background, and detailed results,
please refer to the [full report](pdf/report/report_full.pdf) in the repository.
