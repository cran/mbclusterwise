\name{cw.predict}


\alias{cw.predict}


\title{Prediction procedure for clusterwise multiblock analyses}


\description{Function to perform the prediction of new observations by means of clusterwise multiblock analysis}


\usage{cw.predict(Xnew, res.cw)}


\arguments{
  \item{Xnew}{a data frame containing new observation values for the explanatory variables}
  \item{res.cw}{a list of results created by the function \code{\link{cw.multiblock}}}
}


\value{A list containing the following components is returned:
  \item{clusternew}{a vector containing the new observation assignation to the \eqn{G} expected clusters (when \eqn{G}>1 only)}
  \item{Ypred.cr}{a matrix that contain the predicted dependent values associated with the centered and scaled data for each of the \eqn{G} clusters}
  \item{Ypred.raw}{a matrix that contain the predicted dependent values associated with the raw data for each of the \eqn{G} clusters}
}


\references{Bougeard, S., Abdi, H., Saporta, G., Niang, N., Submitted, Clusterwise analysis for multiblock component methods.}


\author{Stephanie Bougeard (\email{stephanie.bougeard@anses.fr})}


\seealso{\code{\link{cw.multiblock}}, \code{\link{cw.tenfold}}}


\examples{ 
  data(simdata.red) 
  Data.X      <- simdata.red[c(1:10, 21:30), 1:10]
  Data.Y      <- simdata.red[c(1:10, 21:30), 11:13]
  Data.X.test <- simdata.red[c(16:20, 36:40), 1:10]
  ## Note that the options (INIT=2) and (parallel.level = "low") are chosen to quickly
  ## illustrate the function. 
  ## For real data, instead choose (INIT=20) to avoid local optima and (parallel.level = "high")
  ## to improve the computing speed.
  res.cw      <- cw.multiblock(Y = Data.Y, X = Data.X, blo = c(5, 5), option = "none", G = 2,
                  H = 1, INIT = 2, method = "mbpls", Gamma = NULL, parallel.level = "low")
  rescw.pred  <- cw.predict(Data.X.test, res.cw)
}


\keyword{multivariate}
\keyword{cluster}
