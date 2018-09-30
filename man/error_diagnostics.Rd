\name{error_diagnostics}
\alias{error_diagnostics}
\title{
Core function for calculating Geographically weighted error diagnostics
}
\description{
This core function is used for calculation the geographically weighted error diagnostic measures of mean signed deviation, mean absolute error, root mean squared error, and Pearson's correlation coefficient between two variables (predicted and reference values) The function is designed from cov.wt function from the stats package.
}

\usage{
error_diagnostics(x, wts.all, center = TRUE)
}

\arguments{
  \item{x}{
a matrix or data frame. As usual, rows are observations and columns are variables.
}

  \item{wts.all}{
a non-negative and non-zero vector of weights for each observation. Its length must equal the number of rows of x.

}

  \item{center}{
either a logical or a numeric vector specifying the centers to be used when computing. If TRUE, the (weighted) mean of each variable is used, if FALSE, zero is used. If center is numeric, its length must equal the number of columns of x.
  }

}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{rmse }{
root mean squared error
}
  \item{mae}{
  mean absolute error
}
   \item{difference}{
   mean signed deviation
}
   \item{calc.pointn}{
  The number of point used for the calculation
}
   \item{calc.cor.test}{
  Pearson's correlation coefficient
}
   \item{cor.pval}{
  P-value of Pearson's correlation
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
