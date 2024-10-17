# fast.EOT
Fast version (parallelized) to calculate EOTs (Empirical orthogonal functions) in _R_.
It is based on the 'eof' function from [`remote`](https://cran.r-project.org/package=remote), 
which already utilises efficient _C++_ code, but implements an additional parallelization step,
making the calculation of EOTs approximately 7x faster. Currently 'fast.EOT' only calculates 
the EOTs of a single spatio-temporal field, contrary to the 'eot' function in 'remote', which 
also allows to calculate EOTs based on two seperate fields (predictor & predictant). 

**To load** (using `devtools`):
```
library(devtools)
install_github("BernhardKuehn/fast.EOT")
```

**to cite**:
```
citation("fast.EOT")
```

## Non-CRAN required packages

Utilising _C++_ code for speed in a parallelized way, requires the installation 
of an additional _R_ package (`Rcpp2doParallel`) that is not (yet) available on CRAN. 
The `Rcpp2doParallel` _R_ package allows the implementation of _C++_
functions and parallelization using the [`doParallel`](https://cran.r-project.org/package=doParallel)
and [`foreach`](https://cran.r-project.org/package=foreach) backend. Therefore the overhead of compiling _C++_ 
code on each parallel worker is reduced, allowing a faster code execution. 

```{r}
# Non-CRAN
# an Rcpp-extension that handles C++ code in a parallel backend
devtools::install_github("coatless-rd-rcpp/rcpp-and-doparallel")
library("Rcpp2doParallel")
```
