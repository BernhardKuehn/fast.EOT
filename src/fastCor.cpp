#include <Rcpp.h>
using namespace Rcpp;

///////////////// linear model functions //////////////////////////////////
// [[Rcpp::export]]
double corC(NumericVector x, NumericVector y) {
  int nx = x.size(), ny = y.size();
  
  if (nx != ny) stop("Input vectors must have equal length!");
  
  double sum_x = sum(x), sum_y = sum(y);
  
  NumericVector xy = x * y;
  NumericVector x_squ = x * x, y_squ = y * y;
  
  double sum_xy = sum(xy);
  double sum_x_squ = sum(x_squ), sum_y_squ = sum(y_squ);
  
  double out = ((nx * sum_xy) - (sum_x * sum_y)) / sqrt((nx * sum_x_squ - pow(sum_x, 2.0)) * (nx * sum_y_squ - pow(sum_y, 2.0)));
  
  return out;
}


// [[Rcpp::export]]
NumericVector fastCor(NumericMatrix x_i, NumericVector y_i, 
                      bool standardised) {
  // Number of rows of input matrices
  int nrow_x = x_i.nrow(), nrow_y = y_i.size();
  
  // For the current predictor cell, loop through all response cells
  // and calculate the corresponding R-squared value
  NumericVector lm_rsq(nrow_x);
  for (int i = 0; i < nrow_x; i++) {
    
    lm_rsq[i] = pow(corC(x_i(i,_), y_i), 2.0);
      
      // Perform standardisation (optional)
      if (!standardised) {
        lm_rsq[i] = lm_rsq[i] * var(y_i);
      }
      
      // Replace possible NaN with 0
      if (lm_rsq[i] != lm_rsq[i]) {
        lm_rsq[i] = 0;
      }
    }
    
    // return R-squared values of current predictor cell
  return lm_rsq;  
}