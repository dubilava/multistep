# Replication Material

This repository hosts the data and R codes for the paper "A Comparison of Multistep Commodity Price Forecasts Using Direct and Iterated Smooth Transition Autoregressive Methods," accepted in Agricultural Economics. 

## Data

- dataset.RData: this file contains the price series that were accessed on 27 November 2021 from online portals of the World Bank and the International Monetary Fund. Specifically: 
  * The World Bank data are collected from https://www.worldbank.org/en/research/commodity-markets
  * The International Monetary Fund data are collected from https://www.imf.org/en/Research/commodity-prices
- Note: presently (as of February 2022), the IMF online portal is undergoing some maintenance, so the price series don't appear to be available, hopefully the issue will get resolved soon.


## R Codes

### Main

- subset.r -- subsets the commodity price series that exhibit evidence of one-step STAR-type nonlinear dynamics, based on Terasvirta's linearity test, in each expanding window.
- forecast.r -- for the selected commodity price series, runs the expanding window routine to generate one--to--twelve--step--ahead point forecasts for direct and iterated STAR methods and their linear counterparts. The lag order and the STAR function is selected in each expanding window and, in the case of direct STAR, for every horizon.
- evaluate.r -- compares forecasts using root mean square forecast errors and the Diebold--Mariano type test statistics.
- tables.r -- generates Table 1 and Appendix Table 1
- figures.r -- generates Figure 1 and Figure 3

### Auxillary

These codes store functions that facilitate testing and estimating of STAR models

- startest.r
- lstar.r
- estar.r