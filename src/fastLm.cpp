#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List lmC(NumericVector x, NumericVector y) {
  // Check if input vectors have equal length
  int nx = x.size(); // C++ size() is equal to R length()
  int ny = y.size();
  // Throw error message in case of unequal length
  if (nx != ny) stop("Input vectors must have equal length!");

  // Calculate r value based on linear regression formula
  double sum_x = sum(x); // C++ sum() is equal to R sum()
  double sum_y = sum(y);
  
  NumericVector xy = x * y; // Arithmetic operations work the same way as in R
  NumericVector x_squ = x * x;
  NumericVector y_squ = y * y;
  
  double sum_xy = sum(xy);
  double sum_x_squ = sum(x_squ);
  double sum_y_squ = sum(y_squ);
  
  // C++ pow() is equal to R ^
  double lm_r = ((nx * sum_xy) - (sum_x * sum_y)) / sqrt((nx * sum_x_squ - pow(sum_x, 2.0)) * (nx * sum_y_squ - pow(sum_y, 2.0)));

  // Calculate intercept and slope
  double lm_intercept = ((sum_y * sum_x_squ) - (sum_x * sum_xy)) / ((nx * sum_x_squ) - pow(sum_x, 2.0));
  double lm_slope = ((nx * sum_xy) - (sum_x * sum_y)) / ((nx * sum_x_squ) - pow(sum_x, 2.0));

  // Calculate residuals
  NumericVector lm_predicted_values = lm_slope * x + lm_intercept;
  NumericVector lm_residuals = y - lm_predicted_values;
  
  // Calculate test statistics
  int lm_df = ny - 2; // Degrees of freedom 
  
  double lm_stderr = sqrt(sum(pow(lm_residuals, 2.0)) / lm_df) / 
                     sqrt(sum(pow(x - mean(x), 2.0))); // Standard error
                     
  double lm_tscore = lm_slope / lm_stderr; // t score

  
  // List and return single parameters
  return List::create(lm_r, lm_intercept, lm_slope, lm_residuals, lm_tscore, lm_df);
}


// [[Rcpp::export]]
List respLmParam(NumericMatrix x, NumericMatrix y, int cell) {
  int nrow_y = y.nrow();
  
  List lmC_out(nrow_y);
  for (int i = 0; i < nrow_y; i++) {
    lmC_out[i] = lmC(x(cell, _), y(i, _));
  }
  
  return lmC_out;
}
