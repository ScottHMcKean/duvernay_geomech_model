# A Geostatistical Model of Geomechanical Properties in the Duvernay

This repository provides the code and data for reproducing the machine learning 
and geostatistical workflow from our study. 

The code is roughly broken into two parts. The first generates a
machine learning model to predict static geomechanical properties (Young's 
Modulus, Poisson's Ratio, and Brittleness). The second uses that model to predict
geomechanical properties across the Duvernay basin using geostatistics. 

The whole code base is written in R, with a tidyverse philosophy.