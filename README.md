# rrr_kappa

rrr_kappa is an R library that implements Auto-DML method to estimate complier parameters.  (Note that Auto-DML was previously named DML-RRR, which is the reference in the script.)


## Usage
To apply the Auto-DML function, load `primitives.R`, `stage1_kappa.R` and `stage2_kappa.R`.  The `psi` function in `primitives.R` specifies whether we truncate the propensity scores for extreme values.

To replicate the simulation studies that compare estimates by DML and kappa weighting, execute and follow in-line instructions in `main_sim_kappa.R`, which calls various function in this library.  In particular, `specifications_sim.R` varies the simulation design; `stage2_kappa_dml.R` implements both Auto-DML and DML; `naive_kappa.R` implements kappa weighting.

To replicate the analysis of the simulation results, execute `plot_dml_rrr.R`, which calculates the bias, MSE, and coverage probability based on simulated estimates.

To replicate the empirical illustration of Auto-DML in the context of effect of 401(k) on assets, we obtain data and construct transformation of the covariates based the replication files for 
Belloni, A., Chernozhukov, V., Fernández-Val, I., and Hansen, C. (2017). Program Evaluation and Causal Inference With High-Dimensional Data. Econometrica 85, 233-298.

We then execute `main_sipp_kappa.R`, which implements the Auto-DML estimator to estimate compliers' counterfactual outcome distribution.  `specifications_kappa.R' constructs the dictionary.  To replicate the plot of the distribution estimates, execute `main_sipp_plot.R`.

To replicate the empirical illustration of Auto-DML in the context of effect of childbearing on female labor supply, we obtain [data](https://economics.mit.edu/faculty/angrist/data1/data/angev98) and construct transformation of the covariates based the [replication files](https://economics.mit.edu/faculty/angrist/data1/data/angrist3) for (see `m_d_806_extrapolate.do`)
Angrist, J. D. and Fernández-Val, I.,(2013). ExtrapoLATE-ing: External validity and overidentification in the LATE framework. In Advances in Economics and Econometrics

We then execute `main_AE98_kappa.R`, which implements the Auto-DML estimator to estimate compliers' covariate means and conducts the semi-parametric test on complier means.  `specifications_kappa_AE98.R` constructs the dictionary.    

## Authors and acknowledgment
Rahul Singh, Liyang Sun

Paper available on [arXiv](https://arxiv.org/abs/1909.05244).