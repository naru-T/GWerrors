gwerrors <- function (x, vars, fp, adapt = NULL, bw, kernel, longlat = NULL, distMatrix = NULL, random=FALSE)
{
  `%+%` <- function(x,y){ paste0(x,y)}

  p4s <- as.character(NA)
  Polys <- NULL
  fp.missing <- missing(fp)
  if (is(x, "SpatialPolygonsDataFrame")) {
    Polys <- as(x, "SpatialPolygons")
    gridded <- gridded(x)
    dp <- coordinates(x)
    p4s <- proj4string(x)
    data <- as(x, "data.frame")
    if (is.null(longlat) || !is.logical(longlat)) {
      if (!is.na(is.projected(x)) && !is.projected(x)) {
        longlat <- TRUE
      }
      else {
        longlat <- FALSE
      }
    }
  }
  else if (is(x, "SpatialPointsDataFrame")) {
    gridded <- gridded(x)
    dp <- coordinates(x)
    p4s <- proj4string(x)
    data <- as(x, "data.frame")
    if (is.null(longlat) || !is.logical(longlat)) {
      if (!is.na(is.projected(data)) && !is.projected(data)) {
        longlat <- TRUE
      }
      else {
        longlat <- FALSE
      }
    }
  }
  else stop("x must be a Spatial Polygons or Points DataFrame")
  if (is.null(longlat) || !is.logical(longlat))
    longlat <- FALSE
  if (is.integer(vars))
    vars <- names(x)[vars]
  stopifnot(is.character(vars))
  x <- as.matrix(cbind(data[, vars], as.vector(data[, vars[1]] - data[, vars[2]])))
  
  if (any(is.na(x)))
    stop("x contains NAs")
  nc <- ncol(x)
  if (is.null(colnames(x)))
    colnames(x) <- "V" %+% 1:nc
  cn <- colnames(x)
  n1 <- nrow(dp)
  if (n1 != nrow(x))
    stop("Differing row numbers between x and dp")

  if (missing(fp)) {
    fp <- dp
    colnames(fp) <- colnames(dp)
  }
  else {
    stopifnot(is(fp, "Spatial"))
    gridded <- gridded(fp)
  }

  ##NA of pixels in background grid  are excluded
  fp.mask <- !is.na(fp@data[,1])
  fp.na <- is.na(fp@data[,1])
  fp.coord <- coordinates(fp)
  fp.mask.coord <- fp.coord[fp.mask, ]
  fp.na.coord <- fp.coord[fp.na, ]
  
  #coordiniations in all pixels (not masked yet)
  fp <- coordinates(fp)

  n2 <- nrow(fp.mask.coord)
  if (is.null(adapt)) {
    if (!missing(bw))
      bw <- rep(bw, nrow(fp.mask.coord))
    else stop("Bandwidth must be given for non-adaptive weights")
  }
  else {bw <- gw.adapt(dp = dp, fp = fp, quant = adapt, longlat = longlat)

  ###mapping bw###
  ###grid <- SpatialPointsDataFrame(coords=fp, data=data.frame(bw))
  ###spplot(grid)

  bw <- bw[fp.mask]
  }

  ###change here from for loop to apply func.
  if(random==FALSE){


    if (!missing(distMatrix)) {
      dxs.all <- distMatrix
    } else {
      dxs.all <- spDists(dp, fp.mask.coord, longlat = longlat)
      dxs.all[!is.finite(dxs.all)] <- 0
    }


    ###mapping distance (from point 1)###
    ###grid <- SpatialPointsDataFrame(coords=fp.mask.coord, data=data.frame(dxs.all[1, ]))
    ###spplot(grid)

    #weight
    if(kernel == "gaussian")
      wts.all <- apply(dxs.all, 1, function(x){ifelse(x>bw,0, exp(-0.5*(x/bw)^2))})
    else if(kernel == "exponential")
      wts.all <- apply(dxs.all, 1, function(x){ifelse(x>bw,0, exp(-x/bw))})
    else if(kernel == "bisquare")
      wts.all <- apply(dxs.all, 1, function(x){ifelse(x>bw,0, (1-(x/bw)^2)^2)})
    else if(kernel == "tricube")
      wts.all <- apply(dxs.all, 1, function(x){ifelse(x>bw,0, (1-(x/bw)^3)^3)})
    else if(kernel == "boxcar")
      wts.all <- apply(dxs.all, 1, function(x){ifelse(x>bw,0, 1)})
    else stop("ivalid kernel type")
    ###mapping distance (from point 1)###
    ###grid <- SpatialPointsDataFrame(coords=fp.mask.coord, data=data.frame(wts.all[,1]))
    ###spplot(grid)

  } else {
    mcs <- sample(n1)
    dxs.all <- spDists(dp[mcs, ], fp.mask.coord, longlat = longlat)
    dxs.all[!is.finite(dxs.all)] <- 0

    #weight
    if(kernel == "gaussian")
      wts.all <- apply(dxs.all[mcs, ], 1, function(x){ifelse(x>bw,0, exp(-0.5*(x/bw)^2))})
    else if(kernel == "exponential")
      wts.all <- apply(dxs.all[mcs, ], 1, function(x){ifelse(x>bw,0, exp(-x/bw))})
    else if(kernel == "bisquare")
      wts.all <- apply(dxs.all[mcs, ], 1, function(x){ifelse(x>bw,0, (1-(x/bw)^2)^2)})
    else if(kernel == "tricube")
      wts.all <- apply(dxs.all[mcs, ], 1, function(x){ifelse(x>bw,0, (1-(x/bw)^3)^3)})
    else if(kernel == "boxcar")
      wts.all <- apply(dxs.all[mcs, ], 1, function(x){ifelse(x>bw,0, 1)})
    else stop("ivalid kernel type")
  }

  if( any(wts.all[which(wts.all < 0) | is.na(wts.all) | sum(wts.all)==0])){
    stop(print("Invalid weights"))
  }


  cov.wt2_run <- function(y){
    out <- tryCatch(error_diagnostics(as.matrix(x), y) , error =function(e){
      ans <-  c(rmse = NA, mae = NA, difference = NA, calc.pointn = NA, cor.test = NA, cor.pval = NA
                #,mse = mse, bias2 = bias2, variance = v
                )
      return(ans)}
    )}

  ##cov.wt2 includes the calc of rmse & difference
  res.error <- apply(wts.all, 1, function(y)(cov.wt2_run(y)))

  #colname of res.error
  rnm <- c( "rmse_" %+% cn[2] %+% "_" %+% cn[1],
            "mae_" %+% cn[2] %+% "_" %+% cn[1],
            "difference_" %+% cn[2] %+% "_" %+% cn[1],
            "calcpointn_" %+% cn[2] %+% "_" %+% cn[1],
            "cor_" %+% cn[2] %+% "_" %+% cn[1],
            "cor.pvalue_" %+% cn[2] %+% "_" %+% cn[1]
            )
  rownames(res.error) <- rnm

  SDF_fill <- SpatialPointsDataFrame(coords = fp.mask.coord,
                                     data = data.frame(t(res.error)),
                                     proj4string = CRS(p4s))
  if(dim(fp.na.coord)[1]>0){
    SDF_NA <- SpatialPointsDataFrame(coords = fp.na.coord,
                                     data = data.frame(matrix(NA,nrow=dim(fp.na.coord)[1],ncol=dim(t(res.error))[2] )),
                                     proj4string = CRS(p4s))
    names(SDF_NA) <- names(SDF_fill)
    SDF <- rbind(SDF_fill, SDF_NA)
  } else {
    SDF <- SDF_fill
  }

  if (gridded)
    gridded(SDF) <- TRUE
  else if (!is.null(Polys) && fp.missing) {
    df <- data.frame(SDF@data)
    rownames(df) <- sapply(slot(Polys, "polygons"), function(i) slot(i, "ID"))
    SDF <- SpatialPolygonsDataFrame(Sr = Polys, data = df)
  }
  res <- list(SDF = SDF, bandwidth = bw, adapt = adapt)

  res
}
