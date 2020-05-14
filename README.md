# rrr_kappa

rrr_kappa is a MATLAB / R library that implements DML-RRR algorithm to learn complier parameters.


## Usage
For the simulations related to DML and kappa weighting, execute and follow in-line instructions in `master.m`, that calls various function in this library.
For the simulations related to DML-RRR, execute and follow in-line instructions in `main_sim_kappa.m`, that calls various function in this library.

For the empirical illustration of DML-RRR in the context of effect of 401(k) on assets, we obtain data and construct dictionaries from the replication files for 
Belloni, A., Chernozhukov, V., Fernández-Val, I., and Hansen, C. (2017). Program Evaluation and Causal Inference With High-Dimensional Data. Econometrica 85, 233-298.

We then apply DML-RRR algorithm to learn compliers' counterfactual outcome distribution using the relevant section in `main_sipp_kappa.R`.

## Authors and acknowledgment
