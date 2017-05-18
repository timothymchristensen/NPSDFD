# NPSDFD
Replication data and files for Christensen (2017): "Nonparametric Stochastic Discount Factor Decomposition".

### Replicating the empirical application
The file "app_boot.R" reproduces the empirical application. Running the file will import the data, perform estimation, and generate all graphs. Estimates and confidence intervals will be stored in the file "npsdfd_app_output.txt".

### Replicating the simulation exercise
The file "sim_boot.R" will run the simulation exercise and generate all graphs. Change basis from 1 (hermite polynomials) to 2 (B-splines) in line 28 and run again to generate both sets of results. Bias and RMSE will be stored in the file “npsdfd_sim_hp.txt” and "npsdfd_sim_bs.txt".
