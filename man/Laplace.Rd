% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AsyK.R
\name{Laplace}
\alias{Laplace}
\title{Estimated Density Values by Laplace kernel}
\usage{
Laplace(y, k, h)
}
\arguments{
\item{y}{a numeric vector of positive values.}

\item{k}{gird points.}

\item{h}{the bandwidth}
}
\value{
\item{x}{grid points}
        \item{y}{estimated values of density}
}
\description{
Estimated Kernel density values by using Laplace Kernel.
}
\details{
Laplace kernel is developed by Khan and Akbar. Kernel is developed by using Chen's idea. Laplace kernel is;
\deqn{K_{Laplace\left(x,h^{\frac{1}{2}}\right)} (u)=\frac{1}{2\sqrt h}exp \left(-\frac{|{u-x}|}{\sqrt h}\right)}
}
\examples{
y <- rexp(100,1)
h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
Laplace(y,200,h)
}
\references{
Khan, J. A.; Akbar, A. Density Estimation by Laplace Kernel. \emph{Working paper,  Department of Statistics, Bahauddin Zakariya University, Multan, Pakistan.}
}
\seealso{
To examine laplace density plot see \code{\link{plot.Laplace}} and for Mean Squared Error \code{\link{mseLap}}. Similarly, for RIG kernel \code{\link{RIG}}.
}
\author{
Javaria Ahmad Khan, Atif Akbar.
}
