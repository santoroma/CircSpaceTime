# CircSpaceTime
Spatial and Spatio-Temporal Bayesian Model for Circular Data

Implementation of Bayesian models for spatial and spatio-temporal interpolation of circular data using Gaussian Wrapped and Gaussian Projected distributions.

Currently the following models are implemented:  
Spatial Wrapped Normal   
Spatial Projected Normal   
Spatio-Temporal Wrapped Normal   
Spatio-Temporal Projected Normal   

## Installation

### From source
If you are linux/linux-like users or simply you want to compile from source the best way is to use "devtools"

``` r
devtools_installed <- require(devtools)
 if (!devtools_installed){
   install.packages("devtools", dep = TRUE)
    library(devtools)
    }
  install_github("santoroma/CircSpaceTime")  
 ``` 
 
 Dependencies: Rcpp, RcppArmadillo, circular, ggplot2, coda   
 Suggested: foreach, parallel, iterators, doParallel, knitr, rmarkdown, gridExtra   
 
 ### From CRAN
 The package is in submission on CRAN.
 
 ``` r
   install.packages("CircSpaceTime", dep = TRUE)
 ``` 
 
 ## Using the package
 
 ```r
 library(CircSpaceTime)
 ```
 
 For further information on the package you can read the help or take a look at the [vignette](https://github.com/santoroma/CircSpaceTime/tree/master/inst/doc)

## Issues

Please help us to improve the package!  
For any issue/error/"what is this?" report the best way is to visit the [issues page](https://github.com/santoroma/CircSpaceTime/issues) and:
 1. Find if already exist a similar issue, read it and if the case write a precise comment with reproducible example.
 2. If not, open a new one writing a precise comment with reproducible example.
 
 ## Thanks

 Mario, Gianluca and Giovanna
