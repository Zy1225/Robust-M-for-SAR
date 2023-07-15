# On the Robust Estimation of Spatial Autoregressive Models

This repository code contain template code associated with the manuscript "On the Robust Estimation of Spatial Autoregressive Models" published in Econometrics and Statistics by [Tho](https://rsfas.anu.edu.au/about/staff-directory/zhi-yang-tho), [Ding](https://cbe.anu.edu.au/about/staff-directory/dr-ding-ding), [Hui](https://francishui.netlify.app/), [Welsh](https://cbe.anu.edu.au/about/staff-directory/professor-alan-welsh), and [Zou](https://cbe.anu.edu.au/about/staff-directory/dr-tao-zou).

# Getting started

There are currently two directories in this repository:

-   `Code`, which contains `Matlab` scripts implementing the proposed robust M-estimator for the spatial autoregressive (SAR) model along with other scripts for estimating standard errors and performing simulation studies. Many of the functions in these scripts contain pseudo-help files. Additionally, `etoolbox` (Matlab Spatial Econometrics Toolbox by James P. LeSage) is included due to the dependencies of the scripts in using its functions to compute maximum likelihood estimates (MLE) of the SAR model for comparison purpose;

-   `Simulation`, which contains template script to implement the simulation study in the manuscript. Users can run `Simulation_Script_for_W_d_no_contamination.m` to replicate the results from Table 1 in the manuscript.

# If you find any bugs and issues...

If you find something that looks like a bug/issue, please use Github issues and post it up there. As much as possible, please include in the issue:

1.  A description of the bug/issue;
2.  Paste-able code along with some comments that reproduces the problem. If you also have an idea of how to fix the problem, then that is also much appreciated;
3.  Required data files etc...

Alternatively, please contact the corresponding author at [ZhiYang.Tho\@anu.edu.au](mailto:ZhiYang.Tho@anu.edu.au).
