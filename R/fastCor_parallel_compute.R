#' Function to compute Sum of correlations of the whole field in parallel
#' @param X A matrix with rows as space and columns as time.
#' @param Y Same format as X.
#' @param n_cores Specify the number of CPU-cores to use.
#' @importFrom foreach %dopar%

#' @noRd
fastCor_parallel_compute = function(X,Y,n_cores = parallel::detectCores()) {
  # Construct cluster
  cl = parallel::makeCluster(n_cores)

  # After the function is run, shutdown the cluster.
  on.exit(parallel::stopCluster(cl))

  # Register parallel backend
  doParallel::registerDoParallel(cl)
  # Compute estimates
  estimates = foreach::foreach(i = 1:nrow(X), # Perform n simulations
                               .combine = "rbind",           # Combine results
                               # Self-load
                               .packages = c("Rcpp","Rcpp2doParallel")) %dopar% {
                                 #Rcpp::sourceCpp("./functions/fastCor.cpp")
                                 result = fastCor(x_i = as.matrix(X),
                                                  y_i = as.vector(Y[i,]),standardised = F)
                                 result
                               }
  # aggregate results
  out = apply(estimates,2,sum)
  return(out)
}
