error_diagnostics <- function (x, wts.all, center = TRUE)
{
  if (is.data.frame(x))
    x <- as.matrix(x)
  else if (!is.matrix(x))
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x)))
    stop("'x' must contain finite values only")
  n <- nrow(x)
  if (with.wt <- !missing(wts.all)) {
    if (length(wts.all) != n)
      stop("length of 'wt' must equal the number of rows in 'x'")
    if (any(wts.all < 0) || (s <- sum(wts.all)) == 0)
      stop("weights must be non-negative and not all zero")
    wt <- wts.all/s      ###sum(wt). same as swts[i]
  }
  calc.pointn <- length(wt[wt!=0])
  if(calc.pointn <= 5){
    ###if kernel cannot find more than 5 points, results are NA.
    y <- c(rmse = NA, mae = NA, difference = NA, calc.pointn = calc.pointn, cor.test = NA, cor.pval = NA)
    return(y)

  } else {

    if (is.logical(center)) {
      center <- if (center)
        colSums(wt * x)
      else 0
    }
    else {
      if (length(center) != ncol(x))
        stop("length of 'center' must equal the number of columns in 'x'")
    }

    dif <- sum(wt  * x[ ,3], na.rm=TRUE)         #msd
    mae <- sum(wt * abs(x[ ,3]))                 #mae
    rmse <-  sqrt(sum( wt * x[ ,3]^2 ))          #rmse

    xx <- sqrt(wt) * sweep(x, 2, center, check.margin = FALSE)
    c.test <- cor.test(xx[ ,1], xx[ ,2], alternative = "two.sided", method="pearson" )

    y <- c(rmse = rmse, mae = mae, difference = dif,
           calc.pointn = calc.pointn, cor.test = c.test$estimate, cor.pval = c.test$p.value
    )

    return(y)
  }
}
