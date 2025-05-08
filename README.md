# nonlinearEngelcurves
Code (currently codes are available for both Python and Stata, R coming soon) to apply the price index and welfare estimation procedure in Atkin, Faber, Fally and Gonzalez-Navarro (2023).

### Inputs
Requires household consumption data that contain expenditures for households on goods $g$ organized into groups $G$. Household consumption data can either be a repeated cross section or panel. Households are ranked according to total expenditure per capita within a given market identifier. 

### Outputs
Smoothed relative Engel curves, as well as a dataset containing price indices and welfare changes.

### Usage 
Python: Run engel_curve_estimation.ipynb to produce engel curves, compute welfare metrics, and plot them.

Stata: Run engel_curve_estimation.do first to produce engel curve, then  welfare_estimation.do for welfare metrics, and finally plot_welfare_indices.do for plots. Alternately, simply run do run_all.do. 

### Example
This code is set up to estimate Engel curves and compute welfare metrics for an artificial data set I provide here. This artifical example
is based on Global Consumption Data available below, which I use the script "prep_lsms_data.ipynb" to clean and create a 2 period example from.
https://microdata.worldbank.org/index.php/catalog/4424/data-dictionary/f01?file_name=WB_GCD_2010_v2014-03_survey_data.xlsx

### Contact
Please contact jsayre@ucdavis.edu for any bugs discovered, or make a merge request. Suggestions welcome!
