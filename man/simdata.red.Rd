\name{simdata.red}


\alias{simdata.red}


\docType{data}


\title{Simulated toy data with two groups to test the mbclusterwise package}


\description{This data frame is a toy example with a limited number of observations extracted from the data \code{simdata} given in the \code{plspm} package. These are simulated data organized in two clusters showing two different local regression models.}


\usage{data(simdata.red)}



\format{A data frame of simulated data with 40 observations on the following 14 variables.
  \describe{
  \item{\code{mv1}}{first variable of the block Price Fairness (\eqn{X_1})}
  \item{\code{mv2}}{second variable of the block Price Fairness (\eqn{X_1})}
  \item{\code{mv3}}{third variable of the block Price Fairness (\eqn{X_1})}
  \item{\code{mv4}}{fourth variable of the block Price Fairness (\eqn{X_1})}
  \item{\code{mv5}}{fifth variable of the block Price Fairness (\eqn{X_1})}
  \item{\code{mv6}}{first variable of the block Quality (\eqn{X_2})}
  \item{\code{mv7}}{second variable of the block Quality (\eqn{X_2})}
  \item{\code{mv8}}{third variable of the block Quality (\eqn{X_2})}
  \item{\code{mv9}}{fourth variable of the block Quality (\eqn{X_2})}
  \item{\code{mv10}}{fifth variable of the block Quality (\eqn{X_2})}
  \item{\code{mv11}}{first variable of the block Customer Satisfaction (\eqn{Y})}
  \item{\code{mv12}}{second variable of the block Customer Satisfaction (\eqn{Y})}
  \item{\code{mv13}}{third variable of the block Customer Satisfaction (\eqn{Y})}
}}




\references{Esposito Vinzi, V., Ringle, C., Squillacciotti, S. and Trinchera, L. (2007) Capturing and treating unobserved heterogeneity by response based segmentation in PLS path modeling. A comparison of alternative methods by computational experiments. Working paper, ESSEC Business School.}



\examples{
  data(simdata.red)
  simdata.red
  Data.X <- simdata.red[c(1:15, 21:35), 1:10]
  Data.Y <- simdata.red[c(1:15, 21:35), 11:13]  
}


\keyword{multivariate}
\keyword{cluster}
\keyword{datasets}

