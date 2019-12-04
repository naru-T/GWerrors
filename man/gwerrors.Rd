\name{gwerrors}
\alias{gwerrors}
\title{
Geographically weighted error diagnostics
}
\description{
This function calculates the geographically weighted error diagnostic measures of mean signed deviation, mean absolute error, root mean squared error, and Pearson's correlation coefficient between two variables (predicted and reference values) in Spatial*DataFrame. The function is designed from gwss function from the GWmodel package.
}

\usage{
gwerrors(x, vara, fp, adapt = NULL, bw, kernel, longlat = NULL, distMatrix = NULL, random=FALSE)
}

\arguments{
  \item{x}{
a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package sp, which stores referenced validation values.
}

  \item{vars}{
a vector of two variable names which are both prediction and validation to be applied

}

  \item{fp}{
  fitted point defined by Spatial*DataFrame. If NULL, fp is the same locations as x.
  }

  \item{adapt}{
a proportion value (ranged over 0 until 1) to calculate an adaptive kernel. The value is a ratio of the number nearest neighbours (i.e. adaptive distance) against the number of total fitted point found in fp. Default is NULL. If NULL, bw should be specified as non-adaptive (fixed distance) kernel.

}

  \item{bw}{
bandwidth for a fixed distance kernel used in the weighting function. 

}
  \item{kernel}{
function chosen as follows:
gaussian: wgt = exp(-.5*(vdist/bw)^2);
exponential: wgt = exp(-vdist/bw);
bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
boxcar: wgt=1 if dist < bw, wgt=0 otherwise

}

  \item{longlat}{
if TRUE, great circle distances will be calculated. default = NULL.

}
  \item{distMatrix}{
a pre-specified distance matrix, it can be calculated by the function gw.dist. default=NULL.

}

  \item{random}{
This is used for monte carlo permutation test. If TRUE, data values are randomly allocated.

}


}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{SDF }{
a SpatialPointsDataFrame (may be gridded) or SpatialPolygonsDataFrame object (see package “sp”) with rmse, mae, msd (difference_*), number of point used for calculation (calcpointn_*), Pearson's correlation coefficient, p-value of Pearson's correlation.
}
  \item{bandwidth}{
  bandwidth distance at each point for the calculation.
}
   \item{adapt}{
   The size of adaptive kernel
}

}
\references{
Tsutsumida N., Rodríguez-Veiga P., Harris P., Balzter H., Comber A. Investigating Spatial Error Structures in Continuous Raster Data, accepted, International Journal of Applied Earth Observation and Geoinformation.
}
\author{
Tsutsumida N.
}
\note{
See the paper above for the details.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
#TBD
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }% use one of  RShowDoc("KEYWORDS")
\keyword{ multivariate }% __ONLY ONE__ keyword per line
