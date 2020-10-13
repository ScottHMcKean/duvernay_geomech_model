# A Duvernay Geomechanical Model
A machine learning repository provides the code and data for reproducing the 
geomechanical model and geostatistical analysis of geomechanical properties 
in the Duvernay. This code is used for the paper entitled "A Predictive 
Model of Static Elastic Properties in the Duvernay and Waterways Formations".

The code can roughly be broken into two parts. The first trains a
machine learning model to predict static geomechanical properties (Young's 
Modulus, Poisson's Ratio, and Brittleness). The second uses that model to predict
geomechanical properties across the Duvernay basin using a geostatistical approach.

The whole code base is written in R, with a tidyverse philosophy.