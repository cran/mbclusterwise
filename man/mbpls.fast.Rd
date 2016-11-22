\name{mbpls.fast}


\alias{mbpls.fast}


\title{Multiblock partial least squares}


\description{Function to perform a multiblock Partial Least Squares (PLS) of several explanatory blocks \eqn{(X_1, \dots, X_K)} defined as an object of class \code{ktab} (from \code{ade4}), to explain a dependent dataset \eqn{Y} defined as an object of class \code{dudi} (from \code{ade4}). This function is based on the same code and gives the same results as the \code{mbpls} function from the \code{ade4} package with additional ones developed for the clusterwise procedure.}


\usage{mbpls.fast(dudiY, ktabX, scale = FALSE, option = c("none", "uniform"), H)}


\arguments{
  \item{dudiY}{an object of class \code{dudi} (from \code{ade4}) containing the dependent variable(s)}
  \item{ktabX}{an object of class \code{ktab} (from \code{ade4}) containing the blocks of explanatory variables}
  \item{scale}{a logical value indicating whether the explanatory variables should be standardized}
  \item{option}{an option for the block weighting (by default, the first option is chosen): \cr
  	\option{none}{  the block weight is equal to the block inertia}\cr
  	\option{uniform}{  the block weight is equal to \eqn{1/K} for \eqn{(X_1, \dots, X_K)} and to 1 for \eqn{X} and \eqn{Y}}}
  \item{H}{an integer giving the number of dimensions}
}


\value{A list containing the following components is returned:
  \item{crit.reg}{the regression error}
  \item{lX}{a matrix of the global components associated with the whole explanatory dataset (scores of the individuals)}
  \item{XYcoef}{a list of matrices of the regression coefficients of the whole explanatory dataset onto the dependent dataset}
  \item{intercept}{a list of matrices of the regression intercepts of the whole explanatory dataset onto the dependent dataset}
  \item{fitted}{a list of matrices which contain the predicted dependent values}	
}   
 

\references{Bougeard, S., Qannari, E.M., Lupo, C. and Hanafi, M. (2011). From multiblock partial least squares to multiblock redundancy analysis. A continuum approach. Informatica, 22(1), 11-26}


\author{Stephanie Bougeard (\email{stephanie.bougeard@anses.fr})}


\seealso{\code{\link{cw.multiblock}}, \code{\link{cw.tenfold}}, \code{\link{cw.predict}}, \code{\link[ade4]{mbpls}}}



\examples{
  data(simdata.red) 
  Data.X <- simdata.red[c(1:15, 21:35), 1:10]
  Data.Y <- simdata.red[c(1:15, 21:35), 11:13]
  library(ade4)
  dudiy  <- dudi.pca(df = Data.Y, center = FALSE, scale = FALSE, scannf = FALSE)
  ktabx  <- ktab.data.frame(df = data.frame(Data.X), blocks = c(5,5), 
    tabnames = paste("Tab", c(1:2), sep = "."))
  res    <- mbpls.fast(dudiy, ktabx, scale = FALSE, option = "none", H = 2)
}



\keyword{multivariate}