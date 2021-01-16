# Replication Material

This repository hosts the data and R codes for the paper "A Comparison of Direct and Iterated Multistep Smooth Transition Autoregressive Methods for Forecasting Commodity Prices." 

## Data

- dataset.RData: this file contains the price series that were collected from online portals of the World Bank and the International Monetary Fund. 



## R Codes

### Main

- subset.r -- subsets the commodity price series that exhibit evidence of STAR-type nonlinear dynamics based on Terasvirta's linearity test
- select.r -- using the subset of commodity price series, estimates the STAR models for each commodity series; records results such as the functional form of the STAR models (i.e., logistic or exponential) as well as the probability values of the linearity tests.
- forecast.r
- evaluate.r

### Auxillary

These codes store functions that facilitate testing and estimating of STAR models

- startest.r
- lstar.r
- estar.r