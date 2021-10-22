# rrr_kappa

rrr_kappa is an R library that implements Auto-DML method to estimate complier parameters.  (Note that Auto-DML was previously named DML-RRR, which is the reference in the script.)


## Usage
To apply the Auto-DML function, load `primitives.R`, `stage1_kappa.R` and `stage2_kappa.R`.  The `psi` function in `primitives.R` specifies whether we truncate the propensity scores for extreme values.

To replicate the simulation studies that compare estimates by DML and kappa weighting, execute and follow in-line instructions in `main_sim_kappa.R`, which calls various function in this library.  In particular, `specifications_sim.R` varies the simulation design; `stage2_kappa_dml.R` implements both Auto-DML and DML; `naive_kappa.R` implements kappa weighting.

To replicate the analysis of the simulation results, execute `plot_dml_rrr.R`, which calculates the bias, MSE, and coverage probability based on simulated estimates.

To replicate the empirical illustration of Auto-DML in the context of effect of 401(k) on assets, weobtain data and construct transformation of the covariates based the replication files for 
Belloni, A., Chernozhukov, V., Fernández-Val, I., and Hansen, C. (2017). Program Evaluation and Causal Inference With High-Dimensional Data. Econometrica 85, 233-298.

We then execute `main_sipp_kappa.R`, which implements the Auto-DML estimator to estimate compliers' counterfactual outcome distribution.  To replicate the plot of the distribution estimates, execute `main_sipp_plot.R`.

## Authors and acknowledgment
Rahul Singh, Liyang Sun

Paper available on [arXiv](https://arxiv.org/abs/1909.05244).