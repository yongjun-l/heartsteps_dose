# mrt.dose

## Description

The mHealthDose package implements methods to analyze micro-randomized trials (MRTs) with a focus on estimating the optimal dose of just-in-time adaptive interventions (JITAI). The methods are inspired by the HeartSteps study by Klasnja et al. (2019), a mobile health intervention aimed at increasing physical activity through notifications sent at predefined decision points throughout the day. This package includes tools to simulate data, fit models using causal excursion effects, and analyze the treatment effects of various intervention doses.

## Dataset

The original HeartSteps dataset is not included in this package. However, the package provides functions to simulate datasets under different settings that closely resemble the HeartSteps study. Users can customize parameters such as baseline covariates, time-varying covariates, and treatment effects to simulate realistic data for analysis. The original dataset is available [here](https://github.com/klasnja/HeartStepsV1/tree/main?tab=readme-ov-file)

## Installation

You can install the development version of mrt.dose from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("yongjun-l/heartsteps_dose")
```

## Files and Functions

1.  `boot_samples.R`

This file contains functions for generating bootstrap samples for longitudinal data. These functions are integral to performing variance estimation and calculating confidence intervals for model parameters.

2.  `boruvka_functions.R`

This file includes functions for implementing comparison methods based on approaches introduced by Boruvka et al. (2018) These methods can be used to benchmark the causal excursion effects model against alternative approaches.

3.  `estimating_equations.R`

This file houses all the core functionality for fitting the causal excursion effects model. It includes functions to:

-   Define and solve estimating equations for treatment effects.

-   Calculate and print model coefficients, confidence intervals, and p-values.

-   Summarize results for single or multiple doses.

4.  `simulate_mHealth.R`

This file contains functions to simulate data for micro-randomized trials. Key features include:

-   Specification of baseline and time-varying covariates.

-   Customization of treatment effects and correlation structures.

-   Output of simulated datasets suitable for use with the model fitting functions.

Please refer to the vignette for a detailed description of the usage of the functions.
