% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AsyK.R
\name{mseLap}
\alias{mseLap}
\title{Calculate Mean Squared Error( MSE) when Laplace Kernel is used.}
\usage{
mseLap(y, k, h, type)
}
\arguments{
\item{y}{a numeric vector of positive values.}

\item{k}{gird points.}

\item{h}{the bandwidth}

\item{type}{mention distribution of vector.If exponential distribution then use "Exp".
if use gamma distribution then use "Gamma".If Weibull distribution then use "Weibull".}
}
\value{
Mean Squared Error
}
\description{
Calculate MSE by using Laplace Kernel.
}
\examples{
y<-rexp(100,1)
 h<-0.79 * IQR(y) * length(y) ^ (-1/5)
 mseLap(y, 200, h, "Exp")
}
\references{
Khan, J. A.; Akbar, A. Density Estimation by Laplace Kernel. \emph{Working paper,  Department of Statistics, Bahauddin Zakariya University, Multan, Pakistan.}
}
\seealso{
For further MSE by using RIG kernel see \code{\link{mseRIG}}. For density estimation by using Laplace Kernel \code{\link{plot.Laplace}} and for estimated values
of density \code{\link{Laplace}}.
}
