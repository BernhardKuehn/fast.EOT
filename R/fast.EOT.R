#' A parallelized version to calculate  EOTs (Empirical orthogonal teleconnections)
#'
#' This function builds on the work done for the 'remote' package, but parallelizes the code, in order to reduce computation times.
#' Given a spatio-temporal field (as a raster/brick file) it calculates a predefined number of EOTs.
#'
#' @param x The spatio-temporal data, organised in a 'RasterBrick'-format.
#' @param n The number of EOTs to calculate.
#' @param standardised Boolean. Indicates whether the output should be scaled or not.
#' @param n_cores Integer. The number of CPU-cores for the function to use.
#' @param verbose Boolean. Should intermediate messages be printed to the console (TRUE) or not (FALSE). Default to TRUE.
#' @return A List of the EOTs.
#' @references 
#' \bold{Empirical Orthogonal Teleconnections}\cr
#' H. M. van den Dool, S. Saha, A. Johansson (2000)\cr
#' Journal of Climate, Volume 13, Issue 8, pp. 1421-1435\cr
#' \url{http://journals.ametsoc.org/doi/abs/10.1175/1520-0442%282000%29013%3C1421%3AEOT%3E2.0.CO%3B2
#' }
#'  
#' \bold{Empirical methods in short-term climate prediction}\cr
#' H. M. van den Dool (2007)\cr
#' Oxford University Press, Oxford, New York\cr
#' \url{https://global.oup.com/academic/product/empirical-methods-in-short-term-climate-prediction-9780199202782?cc=de&lang=en&}
#' 
#' @exportClass EOTmode
#' @export
fast.eot = function (x, n = 1, standardised,n_cores = parallel::detectCores()-1, verbose = T)
{
  x.vals <- raster::getValues(x)
  orig.var <- calcVar(x, standardised = standardised)
  # preprocessing
  noNA = stats::complete.cases(x.vals)
  coords = raster::coordinates(x)
  coords_noNA = coords[noNA,]
  x.vals_noNA = x.vals[noNA,]
  y.vals_noNA = x.vals_noNA
  # iterate over predefined number of EOT-modes
  out_list = list()
  for(it in 1:n){
  if (verbose) {
    cat("\nCalculating linear model ...", "\n")
  }
  # vectorize & parallelize
  a = fastCor_parallel_compute(X = x.vals_noNA,
                               Y = y.vals_noNA,
                               n_cores = n_cores)

  if (verbose) {
    cat("Locating ", it, ". EOT ...", "\n", sep = "")
  }

  maxxy.all <- which(a == max(a, na.rm = TRUE))
  maxxy <- maxxy.all[1]
  if (length(maxxy.all) != 1) {
    if (verbose) {
      cat("WARNING:", "\n", "LOCATION OF EOT AMBIGUOUS!",
          "\n", "MULTIPLE POSSIBLE LOCATIONS DETECTED, USING ONLY THE FIRST!\n\n")
    }
  }
  if (verbose) {
    cat("Location:", coords_noNA[maxxy,], "\n",
        sep = " ")
  }
  # regress x on x
  x.lm.param.t <- respLmParam(x.vals_noNA, x.vals_noNA, maxxy - 1)
  x.lm.param.p <- lapply(x.lm.param.t, function(i) {
    tmp <- i
    tmp[[5]] <- 2 * stats::pt(-abs(tmp[[5]]), df = tmp[[6]])
    return(tmp)
  })
  rst.x.r <- raster::rasterFromXYZ(data.frame(coords_noNA,sapply(x.lm.param.p, "[[", 1)))
  rst.x.rsq <- raster::rasterFromXYZ(data.frame(coords_noNA,sapply(x.lm.param.p, "[[", 1)^2))
  rst.x.intercept <- raster::rasterFromXYZ(data.frame(coords_noNA,sapply(x.lm.param.p, "[[", 2)))
  rst.x.slp <- raster::rasterFromXYZ(data.frame(coords_noNA,sapply(x.lm.param.p, "[[", 3)))
  rst.x.p <- raster::rasterFromXYZ(data.frame(coords_noNA,sapply(x.lm.param.p, "[[", 5)))
  rst.x.rsq.sums <- raster::rasterFromXYZ(data.frame(coords_noNA,t(a)))
  # get residuals from the linear model
  x.resids <- sapply(x.lm.param.p, "[[", 4)
  rst.x.resids = raster::rasterFromXYZ(data.frame(coords_noNA,t(x.resids)))

  # get the eot-ts
  eot.ts <- as.numeric(x.vals_noNA[maxxy,])
  # regress x on y
  y.lm.param.t <- respLmParam(x.vals_noNA, y.vals_noNA, maxxy - 1)
  y.lm.param.p <- lapply(y.lm.param.t, function(i) {
    tmp <- i
    tmp[[5]] <- 2 * stats::pt(-abs(tmp[[5]]), df = tmp[[6]])
    return(tmp)
  })
  rst.y.r <- raster::rasterFromXYZ(data.frame(coords_noNA,sapply(y.lm.param.p, "[[", 1)))
  rst.y.rsq <- raster::rasterFromXYZ(data.frame(coords_noNA,sapply(y.lm.param.p, "[[", 1)^2))
  rst.y.intercept <- raster::rasterFromXYZ(data.frame(coords_noNA,sapply(y.lm.param.p, "[[", 2)))
  rst.y.slp <- raster::rasterFromXYZ(data.frame(coords_noNA,sapply(y.lm.param.p, "[[", 3)))
  rst.y.p <- raster::rasterFromXYZ(data.frame(coords_noNA,sapply(y.lm.param.p, "[[", 5)))
  rst.y.rsq.sums <- raster::rasterFromXYZ(data.frame(coords_noNA,t(a)))
  # get residuals from the linear model
  y.resids <- sapply(y.lm.param.p, "[[", 4)
  rst.y.resids = raster::rasterFromXYZ(data.frame(coords_noNA,t(y.resids)))

  # calc. variance explained
  resid.var <- calcVar(rst.y.resids, standardised = standardised)
  cum.expl.var <- (orig.var - resid.var)/orig.var
  if (verbose) {
    cat("Cum. expl. variance (%):", cum.expl.var * 100, "\n",
        sep = " ")
  }
  xy <- t(coords_noNA[maxxy,])
  location.df <- as.data.frame(cbind(xy, paste("mode", sprintf("%02.f",
                                                               n), sep = "_"), cum.expl.var, if (length(maxxy.all) !=
                                                                                                 1)
                                                                 "ambiguous"
                                     else "ok"), stringsAsFactors = FALSE)
  names(location.df) <- c("x", "y", "mode", "cum_expl_var",
                          "comment")
  mode(location.df$x) <- "numeric"
  mode(location.df$y) <- "numeric"
  mode(location.df$cum_expl_var) <- "numeric"
  out <- methods::new("EotMode", 
                      mode = as.integer(it),
             name = paste("mode", sprintf("%02.f",it), sep = "_"),
             eot = eot.ts, 
             coords_bp = xy, 
             cell_bp = maxxy,
             cum_exp_var = cum.expl.var,
             r_predictor = rst.x.r,
             rsq_predictor = rst.x.rsq,
             rsq_sums_predictor = rst.x.rsq.sums,
             int_predictor = rst.x.intercept,
             slp_predictor = rst.x.slp,
             p_predictor = rst.x.p,
             resid_predictor = rst.x.resids,
             r_response = rst.y.r,
             rsq_response = rst.y.rsq,
             int_response = rst.y.intercept, 
             slp_response = rst.y.slp,
             p_response = rst.y.p,
             resid_response = rst.y.resids
             )
  if(n > 1){
    # calculate more EOTs based on the remaining residuals
    y.vals_noNA = t(y.resids)
  }
    # free memory
    rm(list = c("eot.ts", "maxxy", "location.df",
                "rst.x.r", "rst.x.rsq", "rst.x.rsq.sums", "rst.x.intercept",
                "rst.x.slp", "rst.x.p", "rst.x.resids"))
    gc()
  out_list[[it]] = out
  }
  names(out_list) = paste0("EOT",1:n)
  return(out_list)
}
