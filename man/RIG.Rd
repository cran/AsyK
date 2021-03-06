% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AsyK.R
\name{RIG}
\alias{RIG}
\title{Estimated Density Values by Resiprocal Inverse Gaussian kernel}
\usage{
RIG(y, k, h)
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
Estimated Kernel density values by using Resiprocal Inverse Gaussian Kernel.
}
\details{
Scaillet 2003. proposed Resiprocal Inverse Gaussian kerenl. He claimed that his proposed kernel share the same properties as those of gamma kernel estimator.
\deqn{K_{RIG \left( \ln{ax}4\ln {(\frac{1}{h})} \right)}(y)=\frac{1}{\sqrt {2\pi y}}  exp\left[-\frac{x-h}{2h} \left(\frac{y}{x-h}-2+\frac{x-h}{y}\right)\right]}
}
\examples{
y <- rexp(100,1)
h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
RIG(y,200,h)
}
\references{
Scaillet, O. 2004. Density estimation using inverse and reciprocal inverse Gaussian kernels. 	 \emph{Nonparametric Statistics}, \strong{16}, 217-226.
}
\seealso{
To examine RIG density plot see \code{\link{plot.RIG}} and for Mean Squared Error \code{\link{mseRIG}}. Similarly, for Laplace kernel \code{\link{Laplace}}.
}
\author{
Javaria Ahmad Khan, Atif Akbar.
}
