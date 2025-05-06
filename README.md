# nonlinearEngelcurves
Code (currently replication codes are available for both Python and Stata, R coming soon) to replicate the nonlinear price index and welfare estimation procedure in Atkin, Faber, Fally and Gonzalez-Navarro (2023).

### Inputs
Requires household consumption data that contains expenditures for households on goods $g$ organized into groups $G$. Household consumption data can either be a repeated cross section or panel, but in either case must only contain two time periods. Households are ranked according to expenditures within a given market, which also must be provided. 

### Outputs
Smoothed expenditure and Engel curves, as well as a dataset containing welfare indices.

### Usage 
Python: Run engel_curve_estimation.ipynb to produce engel curves, compute welfare metrics, and plot them.

Stata: Run engel_curve_estimation.do first to produce engel curve, then  welfare_estimation.do for welfare metrics, and finally plot_welfare_indices.do for plots. 

### Example
This code is already set up to estimate Engel curves and compute welfare metrics using LSMS data available here:
https://github.com/lsms-worldbank/LSMS-ISA-harmonised-dataset-on-agricultural-productivity-and-welfare
We reprovide @TBentze 's "Household_dataset.dta" here to generate example output.

### Contact
Please contact jsayre@ucdavis.edu for any bugs discovered, or make a merge request. Suggestions welcome!
