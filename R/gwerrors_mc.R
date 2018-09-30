###function for gw.error
###edited from gwss.montecarlo
gwerrors_mc <- function(x, vars, fp, adapt = NULL, bw, longlat = NULL, distMatrix = NULL, nsim = 99)
{
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
  x <- as.matrix(data[, vars])
  if (any(is.na(x)))
    stop("x contains NAs")
  nc <- ncol(x)
  if (is.null(colnames(x)))
    colnames(x) <- paste("V", 1:nc, sep = "")
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
  fp.mask <-  which(!is.na(fp@data[,1]))
  #fp.mask <- !is.na(fp@data[,1])
  fp.coord <- coordinates(fp)
  fp.mask.coord <- fp.coord[fp.mask,]



  n2 <- nrow(fp.mask.coord)
  if (is.null(adapt)) {
    if (!missing(bw))
      bw <- rep(bw, nrow(fp.mask.coord))
    else stop("Bandwidth must be given for non-adaptive weights")
  }
  else {bw <- gw.adapt(dp = dp, fp = coordinates(fp), quant = adapt, longlat = longlat)

  ###mapping bw###
  ###grid <- SpatialPointsDataFrame(coords=fp, data=data.frame(bw))
  ###spplot(grid)

  bw <- bw[fp.mask]
  }

  if (!missing(distMatrix)) {
    dxs.all <- distMatrix
  } else {
    dxs.all <- spDists(dp, fp.mask.coord, longlat = longlat)
    dxs.all[!is.finite(dxs.all)] <- 0
  }


  dataI <- SpatialPointsDataFrame(dp, data.frame(x))
  SD <- SpatialPointsDataFrame(coords = fp.mask.coord,data.frame(val=rep(1, length(fp.mask.coord[,1]))), proj4string = CRS(p4s))

  resi <- (gwerrors(dataI, vars=vars, kernel = "bisquare", fp = SD, adapt = adapt, longlat = FALSE, distMatrix = dxs.all, random=FALSE))$SDF


  dp.n <- dim(dp)[1]
  nstats <- ncol(resi)
  lss.i <- array(0, dim = c(dim(resi)[1], dim(resi)[2], nsim + 1))
  lss.i[, , 1] <- as.matrix(resi@data)
  #    dMat1 <- dMat
  for (i in 1:nsim) {
    mcs <- sample(dp.n)
    dataI <- SpatialPointsDataFrame(dp[mcs, ], data.frame(x),  match.ID = F)
    # dMat1[mcs, ] <- dMat[1:dp.n, ]
    # dMat1[, mcs] <- dMat[, 1:dp.n]
    #   resi <- (gw.error(dataI, vars=vars, fp = SD, adapt = adapt, longlat = FALSE))$SDF

    resi <- (gwerrors(dataI, vars=vars,  kernel = "bisquare",fp = SD, adapt = adapt, longlat = FALSE,  random=TRUE))$SDF

    lss.i[, , i + 1] <- as.matrix(resi@data)
  }
  dimnames(lss.i)[[1]] <- seq(0, dim(resi)[1] - 1)
  dimnames(lss.i)[[2]] <- names(resi)
  dimnames(lss.i)[[3]] <- c("Original_Data", paste("Sample",
                                                   seq(1:nsim), sep = "_"))
  test <- apply(lss.i, c(1, 2), function(x, n) 1 - rank(x,
                                                        ties.method = "first")[1]/n, n = nsim + 1)
  test
}
