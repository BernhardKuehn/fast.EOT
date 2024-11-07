#' Function to plot the EOT pattern (spatial/temporal)
#'
#' This functions is a rather easy plotting function to see the spatial & temporal 
#' representation of the first n calculated EOT modes. 
#' Typically a maximum of 15 modes should not be exceeded. 
#' 
#' 
#' @param eot.out An output object from the \code{fast.EOT} function
#' @param n The number of EOTs to show.
#' @param type Indicates which spatial representation should be shown. Either the regression weights as 
#' in the original work of van Dool 2000 with \code{"regression weights"} 
#' or squared correlation coefficient per grid cell with \code{"R2"}. 

#' @export
plot_eot = function(eot.out,n,type = "regression weights"){
  opar = graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
  # get eot time series
  EOTs_temporal = sapply(eot.out,function(x) x@eot)
  # get spatial pattern (= slope or regression weights/R2)
  if(type == "regression weights"){
  EOTs_spatial = raster::brick(lapply(eot.out,function(x) x@slp_predictor))
  } else if(type == "R2"){
    EOTs_spatial = raster::brick(lapply(eot.out,function(x) x@rsq_predictor))
  } else{
    stop("'type' must be specified as either 'regression weights' or 'R2'!")
  }
  
  graphics::par(mfrow = c(1,1),mar = c(2,4,1,3))
  # create plot layout 
  tmp = matrix(1:(2*n),nrow = 2)
  if(n > 3){
  frst.row = tmp[,1:ceiling(n/2)]
  scnd.row = tmp[,(1+ceiling(n/2)):n]
  if(ncol(frst.row) > ncol(scnd.row)){
    scnd.row = cbind(scnd.row,c(max(scnd.row)+1,max(scnd.row)+2))
  }
  mat = rbind(frst.row,scnd.row)
  } else if(n >= 9){
    indx.split = split(1:n,f = rep(1:3,each = ceiling(n/3))[1:n])
    frst.row = tmp[,indx.split[[1]]]
    scnd.row = tmp[,indx.split[[2]]]
    thrd.row = tmp[,indx.split[[3]]]
    if(ncol(scnd.row) > ncol(thrd.row)){
      thrd.row = cbind(thrd.row,c(max(thrd.row)+1,max(thrd.row)+2))
    }
    mat = rbind(frst.row,scnd.row,thrd.row)
  } else {
  mat = tmp
  }
  graphics::layout(mat,
         height = rep(c(1.5,1),nrow(mat)/2))
  
  for (i in 1:n) {
      terra::image(EOTs_spatial[[i]],
                   col = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(10, 
                                                                                                   "RdYlBu")))(100),
                   las = 1,xlab = "",ylab = "lat")
      terra::contour(x = EOTs_spatial[[i]], add = TRUE, 
                     col = "gray30",labcex = 0.8,vfont = c("sans serif","bold"))
      # add point of eot ts
      graphics::points(eot.out[[i]]@coords_bp,pch = 24,bg = "deeppink4",cex = 1.1,
             col = "wheat")
      maps::map(add = TRUE, fill = TRUE, col = "gray90")
    graphics::box()
    graphics::title(paste0("EOT", i),adj = 0)
    graphics::mtext("lon",side = 1,line = 1.5,cex = 0.8)
    # add time series
    plot(EOTs_temporal[,i],type = "l",xlab = "year",
         ylab = "",las = 1)    
   }
}
  