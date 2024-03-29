% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AsyK.R
\name{RIG}
\alias{RIG}
\title{Estimated Density Values by Reciprocal Inverse Gaussian kernel}
\usage{
RIG(x = NULL, y, k = NULL, h = NULL)
}
\arguments{
\item{x}{scheme for generating grid points}

\item{y}{a numeric vector of positive values.}

\item{k}{gird points.}

\item{h}{the bandwidth}
}
\value{
\item{x}{grid points}
        \item{y}{estimated values of density}
}
\description{
Estimated Kernel density values by using Reciprocal Inverse Gaussian Kernel.
}
\details{
Scaillet 2003. proposed Reciprocal Inverse Gaussian kerenl. He claimed that his proposed kernel share the same properties as those of gamma kernel estimator.
\deqn{K_{RIG \left( \ln{ax}4\ln {(\frac{1}{h})} \right)}(y)=\frac{1}{\sqrt {2\pi y}}  exp\left[-\frac{x-h}{2h} \left(\frac{y}{x-h}-2+\frac{x-h}{y}\right)\right]}
}
\examples{
#Data can be simulated or real data
## Number of grid points "k" should be at least equal to the data size.
### If user define the generating scheme of gridpoints than number of gridpoints should
####be equal or greater than "k"
###### otherwise NA will be produced.
y <- rexp(100, 1)
xx <- seq(min(y) + 0.05, max(y), length = 100)
h <- 2
den <- RIG(x = xx, y = y, k = 200, h = h)

##If scheme for generating gridpoints is unknown
y <- rexp(50, 1)
h <- 3
den <- RIG(y = y, k = 90, h = h)

\dontrun{
##If user do not mention the number of grid points
y <- rexp(23, 1)
xx <- seq(min(y) + 0.05, max(y), length = 90)
#any bandwidth can be used
require(KernSmooth)
h <- dpik(y)
den <- RIG(x = xx, y = y, h = h)
}
#if bandwidth is missing
y <- rexp(100, 1)
xx <- seq(min(y) + 0.05, max(y), length = 100)
den <- RIG(x = xx, y = y, k = 90)
}
\references{
Scaillet, O. 2004. Density estimation using inverse and reciprocal inverse Gaussian kernels. 	 \emph{Nonparametric Statistics}, \strong{16}, 217-226.
}
\seealso{
To examine RIG density plot see \code{\link{plot.RIG}} and for Mean Squared Error \code{\link{mse}}. Similarly, for Laplace kernel \code{\link{Laplace}}.
}
\author{
Javaria Ahmad Khan, Atif Akbar.
}
