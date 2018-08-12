# CircSpaceTime
An R package for spatial, spatio-temporal and temporal model for circular, cylindrical and spherical data

This package implements spatial, spatio-temporal and temporal models for circular data.  

Currently the following models are implemented:  
Spatial Wrapped Normal   
Spatial Projected Normal   

Yet to come:      
Spatio-Temporal Wrapped Normal   
Spatio-Temporal Projected Normal   

We are going to update constantly the library with Bayesian and classical models dealing with complex dependence structures for circular, cylindrical and spherical variable.

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
 Suggested: foreach, doParallel, knitr, rmarkdown, gridExtra   
 
 ### From binary
 The package will be released on CRAN before the end of the 2018.
 For now you can download (win or mac) package [here](https://github.com/santoroma/CircSpaceTime/binary)
 
 Dependencies: circular, ggplot2, coda   
 Suggested: foreach, doParallel, knitr, rmarkdown, gridExtra   
 
 ## Using the package
 
 ```r
 library(CircSpaceTime)
 ```
 
 For further information on the package you read the help or take a first look at the [vignette](https://github.com/santoroma/CircSpaceTime/blob/master/inst/doc/CircSpaceTime.html)

## Issues

Please help us to improve the package!!!  
For any issue/error/"what is this?" report the best way is to visit the [issues page](https://github.com/santoroma/CircSpaceTime/issues) and:
 1. Find if already exist a similar issue, read it and if the case write a precise comment with reproducible example
 2. If not, open a new one writing a precise comment with reproducible example.
 
 ## Thanks

 Mario, Gianluca and Giovanna
