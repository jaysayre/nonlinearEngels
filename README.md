# nonlinearEngelcurves
Code (available for Python, Stata, R, and Julia) to apply the price index and welfare estimation procedure in "Measuring Welfare and Inequality with Incomplete Price Information" (with David Atkin, Ben Faber and Marco Gonzalez-Navarro), *Quarterly Journal of Economics*, February 2024.

### Inputs
Requires household consumption data that contain expenditures for households on goods $g$ organized into groups $G$. Household consumption data can either be a repeated cross section or panel. Households are ranked according to total expenditure per capita within a given market identifier. 

### Outputs
Smoothed relative Engel curves, as well as a dataset containing price indices and welfare changes.

### Usage
Python: Run engel_curve_estimation.ipynb to produce engel curves, compute welfare metrics, and plot them.

R: Run engel_curve_estimation.R from the R_package directory. It sources engel_curves.R (function library) and produces engel curves, welfare metrics, and plots. Supports exact price correction (AFFG Proposition 1) and first-order price correction (AFFG Equation 8).

Stata: Run engel_curve_estimation.do first to produce engel curve, then  welfare_estimation.do for welfare metrics, and finally plot_welfare_indices.do for plots. Alternately, simply run do run_all.do.

Julia: Include engel_curves.jl and import the EngelCurves module. Requires DataFrames.jl, Statistics, and LinearAlgebra (all in the standard library or commonly installed). Provides the same functions as the Python/R packages: lpoly, monotonicity_tails, monotonicity_check, identify_horizontal_shifts, gen_welfare_df, identify_non_crossings, and both exact and first-order price corrections.

### Example
This code is set up to estimate Engel curves and compute welfare metrics for an artificial data set I provide here. This artifical example
is based on Global Consumption Data available below, which I use the script "prep_lsms_data.ipynb" to clean and create a 2 period example from.
https://microdata.worldbank.org/index.php/catalog/4424/data-dictionary/f01?file_name=WB_GCD_2010_v2014-03_survey_data.xlsx

### Contact
Please contact jsayre@ucdavis.edu for any bugs discovered, or make a merge request. Suggestions welcome!
