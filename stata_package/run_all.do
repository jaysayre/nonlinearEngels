// run_all.do
clear all
set more off

do engel_curve_estimation.do
do welfare_estimation.do
do plot_welfare_indices.do
