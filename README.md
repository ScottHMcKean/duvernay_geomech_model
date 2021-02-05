# A Geostatistical Model of Geomechanical Properties in the Duvernay

This repository provides the code and data for reproducing the machine learning 
and geostatistical workflow from our study. 

The code is roughly broken into two parts. The first generates a
machine learning model to predict static geomechanical properties (Young's 
Modulus, Poisson's Ratio, and Brittleness). The second uses that model to predict
geomechanical properties across the Duvernay basin using geostatistics. 

The whole code base is written in R, with a tidyverse philosophy.

## Machine Learning Model

The data in the /data folder is a sanitized version of interpreted geomechanical
tests from the basin. Locations and well identification have been preserved in a 
spirit of fairness to operators and data providers (these tests cost a lot). The
tests include:
- Brazilian tests (brazilian_tests.csv)
- Unconfined compressive strength / Triaxial testing (static_tests.csv)
- Ultrasonic pulse transmission tests (upt_tests.csv)
- Merged static / dynamic test results with or without residuals (merged_tests_*.csv)

## Geostatistical Analysis

The /data/dvisopach folder contains a copy of the data from Switzer et al. (2008)'s
"Duvernay Interval Isopach and Lithofacies", found here:http://www.ags.gov.ab.ca/publications/DIG/ZIP/DIG_2008_0124.zip.

The /data/shen2018 folder contains a copy of the data from 
Shen et al. (2018)'s "In-Situ Stress Measurements for the Duvernay Formation, Alberta", found here: http://ags.aer.ca/document/DIG/DIG_2018_0013.zip. Vertical stress, 
minimum horizontal stress, and pore-pressure measurements are used in the study.

The /data/ags3dmodel folder contains layers from the Alberta Geological Survey's
3D Geological Model: https://ags.aer.ca/publication/3d-pgf-model-v2 that are loaded
by the study.



