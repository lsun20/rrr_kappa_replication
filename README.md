# rrr_kappa

rrr_kappa is a MATLAB / R library that implements Auto-DML algorithm to learn complier parameters.  (Note that Auto-DML was previously named DML-RRR, which is the reference in the script.)


## Usage
For the simulations related to DML and kappa weighting, execute and follow in-line instructions in `master.m`, that calls various function in this library.
For the simulations related to Auto-DML, execute and follow in-line instructions in `main_sim_kappa.m`, that calls various function in this library. 

For the empirical illustration of Auto-DML in the context of effect of 401(k) on assets, we obtain data and construct dictionaries from the replication files for 
Belloni, A., Chernozhukov, V., Fernández-Val, I., and Hansen, C. (2017). Program Evaluation and Causal Inference With High-Dimensional Data. Econometrica 85, 233-298.

We then apply Auto-DML algorithm to learn compliers' counterfactual outcome distribution using the relevant section in `main_sipp_kappa.R`.

## Authors and acknowledgment
