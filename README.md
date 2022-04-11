# Replication Material

This repository hosts the data and R codes for the paper "[A Comparison of Multistep Commodity Price Forecasts Using Direct and Iterated Smooth Transition Autoregressive Methods](https://onlinelibrary.wiley.com/doi/full/10.1111/agec.12707)." 

The codes should run with no issue and replicate the results of the study as long as all the supplied material (i.e., the data file as well as the R codes) are stored in the same folder.

## Data

- `dataset.RData` file contains the price series that were accessed on 27 November 2021 from online portals of the World Bank and the International Monetary Fund. Specifically: 
  * The World Bank data are collected from https://www.worldbank.org/en/research/commodity-markets
  * The International Monetary Fund data are collected from https://www.imf.org/en/Research/commodity-prices
- Note: The IMF Primary Commodity Price System online portal has been undergoing some maintenance, so the price series don't appear to be available. Hopefully the issue will get resolved soon.


## R Codes

### Main

- `01-subset.r` -- subsets the commodity price series that exhibit evidence of one-step STAR-type nonlinear dynamics, based on Terasvirta's linearity test, in each expanding window.
- `02-forecast.r` -- for the selected commodity price series, runs the expanding window routine to generate one--to--twelve--step--ahead point forecasts for direct and iterated STAR methods and their linear counterparts. The lag order and the STAR function is selected in each expanding window and, in the case of direct STAR, for every horizon.
- `03-evaluate.r` -- compares forecasts using root mean square forecast errors and the Diebold--Mariano type test statistics.
- `04-tables.r` -- generates Table 1 and Appendix Table 1
- `05-figures.r` -- generates Figure 1 and Figure 3

### Auxillary

These codes store functions that facilitate testing and estimating of STAR models

- `startest.r`
- `lstar.r`
- `estar.r`


## License

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)