# CoronavirusGraphs_Fits

R code to make country specific plots of the covid19 data from Johns Hopkins datasets. The goal is to be able to plot data from different locations to visualize the similarities and differences in the spread.

I've also included a very simple Bayesian model fitting approach where we simply assume that new cases are Poisson distributed with a rate that is the product of the number of current cases and the per case R0. Because we are not using and mostly do not have data on the probability of testing given infected, this approach basically assumes that a constant fraction of infected cases become identified.  

The disease incidence data are from Johns Hopkins CSSE available at https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series


NOTE: Johns Hopkins changed their data format.

Two other sources of data are now being used here. One is the NY Times dataset. The other is from covidtracker.com

Added is an Rmd file that looks at state and county level data.
