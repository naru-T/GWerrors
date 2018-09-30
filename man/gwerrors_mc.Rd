\name{gwerrors_mc}
\alias{gwerrors_mc}
\title{
Monte Carlo permutation test for Geographically weighted error diagnostics
}
\description{
This function applys Monte Carlo permutation test to the geographically weighted error diagnostic measures of mean signed deviation, mean absolute error, root mean squared error, Pearson's correlation coefficient between two variables (predicted and reference values) in Spatial*DataFrame. The function is designed from gwss function from the GWmodel package.
}
\usage{
gwerrors_mc(x, vars, fp, adapt = NULL, bw, longlat = NULL, distMatrix = NULL, nsim = 99)
}

\arguments{
  \item{x}{
a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package sp
}

  \item{vars}{
a vector of variable names to be applied

}

  \item{fp}{
  fitted point defined by Spatial*DataFrame. If NULL, fp is the same locations as x.
  }

  \item{adapt}{
if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is NULL, where a fixed kernel is found (bandwidth is a fixed distance)

}

  \item{bw}{
bandwidth used in the weighting function. if adaptive bandwidth is applied, bw should be proportional (eg, 0.2 when considering 20% nearest neighbours of the total data points). If fixed bandwidth, unit is 'km'.

}

  \item{longlat}{
if TRUE, great circle distances will be calculated. default = NULL.

}
  \item{distMatrix}{
a pre-specified distance matrix, it can be calculated by the function gw.dist. default=NULL.

}

  \item{nsim}{
The number of simulation time. The defaul is 99.
}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
P-values of each measure of mean signed deviation, mean absolute error, root mean squared error, Pearson's correlation coefficient are found at each column.

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
