# Tesi-KI

# README: Simulation Study to Compare Different Methods for Multi-State Models with Panel Data

## Project Description

This project implements a pipeline to simulate data and analyze approaches for fitting multi-state models with panel data under both Markovian and Semi-Markovian assumptions. The chosen multi-state model is an irreversible illness-death model, where the intermediate state represents dementia. Each transition is influenced by a set of covariates with different magnitudes.

The simulation process generates different scenarios to test various models by varying the dataset sample size and the observation scheme under which individuals are observed.

The workflow consists of four main phases:

0.  **Data Simulation**: Generating data under different assumptions based on a pre-existing ground truth.
1.  **Model Fitting**: Estimating various multi-state models for each simulated scenario.
2.  **Performance Evaluation**: Assessing model performance using statistical metrics.
3.  **Results Visualization**: Creating plots to analyze performance measures.

------------------------------------------------------------------------

## File Structure

The code is organized into the following main scripts:

### **Data Simulation**

-   `0.MM_Simulation_exe.R`: Simulates data under the Markovian assumption.
-   `0.SM_Simulation_exe.R`: Simulates data under the Semi-Markovian assumption.
-   Both scripts load simulation functions, ground truth parameters, and generate datasets of different sizes (500, 2K, 5K, 10K patients), saving them in the `Simulated_data_MM` or `Simulated_data_SM` folder.

### **Model Fitting**

-   `1.MM_fit_exe.R`: Fits various multi-state models to Markovian simulated data.
-   `1.SM_fit_exe.R`: Fits various multi-state models to Semi-Markovian simulated data.
-   Uses specific functions for data preparation and model fitting, saving the results in the `wrapper_MM` or `wrapper_SM` folder according to the Markovian or Semi-Markovian simulation assumption. Inside this folder, the results are saved in different subfolders categorized by sample size (`results_"size"`), which further contain results nested by observation scheme and seed.

### **Performance Evaluation**

-   `2.MM_fit_exe.R`: Computes performance metrics for multi-state models fitted to Markovian simulated data.
-   `2.SM_fit_exe.R`: Computes performance metrics for multi-state models fitted to Semi-Markovian simulated data.
-   Uses specific functions for performance evaluation, saving the results in the `wrapper_MM` or `wrapper_SM` folder. The results are stored in different folders categorized by sample size (`saved_performance_"size"`), further divided by observation scheme and seed.

The following files store different performance measures for all models analyzed and the benchmark model:

-   `convergence.RData`: Convergence level for each model for each seed.
-   `all_estimates.RData`: Estimates of baseline parameters, covariate effects, and HR for each seed.
-   `mean_estimates.RData`: Average estimates of baseline parameters, covariate effects, and HR over all seeds.
-   `bias_all.RData`: Absolute bias of baseline parameters, covariate effects, and HR for each seed.
-   `res_bias.RData`: Mean and confidence intervals of absolute bias for baseline parameters and covariate effects.
-   `baseline_bias.RData`: Mean and confidence intervals of relative bias for baseline parameters and covariate effects.
-   `mean_estimates_baseline.RData`: Average estimates of baseline parameters over all seeds.
-   `baseline_estimates_all.RData`: Estimates of baseline parameters for each seed.
-   `all_coverage.RData`: Coverage of baseline parameters and covariate effects for each seed.
-   `95%coverage.RData`: Mean coverage and confidence intervals for baseline parameters and covariate effects.
-   `se_all.RData`: Standard errors of covariate effects for each seed.
-   `mean_se.RData`: Mean standard error of covariate effects.
-   `comp_time.RData`: Mean computational time.
-   `significancy.RData`: Mean probability that covariates are predicted to be significant.
-   `significancy_all.RData`: Whether each covariate is predicted to be significant for each seed.
-   `significancy_covs.RData`: Presence/absence of covariate effects in the simulation process.
-   `all_estimates_lfe.RData`: Average time spent in each state for each seed.
-   `mean_estimates_lfe.RData`: Average time spent in each state averaged over seeds.
-   `bias_lfe.RData`: Bias of the average time spent in each state, averaged over seeds.
-   `power.RData`: Mean statistical power for each transition.
-   `typeIerr.RData`: Mean Type I error for each feasible transition.

### **Results Visualization**

-   `3.MM_plots.R`: Markdown report to visualize key results for different models compared across observation schemes and sample sizes from data generaed under Markov assumption.
-   `3.SM_plots.R`: Markdown report to visualize key results for different models compared across observation schemes and sample sizes from data generaed under Semi-Markov assumption.

------------------------------------------------------------------------

## Running the Code

To execute the pipeline:

1.  **Data Simulation**:

    ``` r
    source("0.MM_Simulation_exe.R") # For Markovian data
    source("0.SM_Simulation_exe.R") # For Semi-Markovian data
    ```

2.  **Model Fitting**: The process must be repeated for each sample size (500, 2K, 5K, 10K) and set manually by the user:

    ``` r
    source("1.MM_fit_exe.R") # Fit Markovian models
    source("1.SM_fit_exe.R") # Fit Semi-Markovian models
    ```

3.  **Performance Evaluation and Visualization**: The process must be repeated for each sample size (500, 2K, 5K, 10K) and set manually by the user:

    ``` r
    source("2.MM_performance_exe.R")
    source("2.SM_performance_exe.R")
    ```

------------------------------------------------------------------------

## Requirements

-   **Required R Libraries**: data.table, deSolve, dplyr, elect, flexsurv, fs, future, future.apply, ggplot2, haven, here, hesim, kableExtra, lubridate, mstate, msm, nhm, parallel, plotly, survival, tidyr, webshot

## Contact

For questions or issues, contact: [alessandra.pescina\@gmail.com](mailto:alessandra.pescina@gmail.com){.email}
