#' AsyK
#'
#' Kernel Density Estimation
#'
#' @description  A collection of functions related to density estimation by using Chen's (2000) idea. For observing estimated values see \code{\link{Laplace}} and \code{\link{RIG}}. Plots by using these kernels can be drawn by \code{\link{plot.Laplace}} and \code{\link{plot.RIG}}. Mean squared errors (MSE) can be calculated by \code{\link{mse}}.
#' Here we also present a normal scale rule bandwidth which is given by Silverman (1986) for non-normal data.
#'@author Javaria Ahmad Khan, Atif Akbar.
"_PACKAGE"
#' Estimated Density Values by Reciprocal Inverse Gaussian kernel
#'
#' Estimated Kernel density values by using Reciprocal Inverse Gaussian Kernel.
#' @details Scaillet 2003. proposed Reciprocal Inverse Gaussian kerenl. He claimed that his proposed kernel share the same properties as those of gamma kernel estimator.
#' \deqn{K_{RIG \left( \ln{ax}4\ln {(\frac{1}{h})} \right)}(y)=\frac{1}{\sqrt {2\pi y}}  exp\left[-\frac{x-h}{2h} \left(\frac{y}{x-h}-2+\frac{x-h}{y}\right)\right]}
#' @param y a numeric vector of positive values.
#' @param x scheme for generating grid points
#' @param k gird points.
#' @param h the bandwidth
#' @import stats
#' @examples
#' #Data can be simulated or real data
#' ## Number of grid points "k" should be at least equal to the data size.
#' ### If user define the generating scheme of gridpoints than number of gridpoints should
#' ####be equal or greater than "k"
#' ###### otherwise NA will be produced.
#' y <- rexp(100, 1)
#' xx <- seq(min(y) + 0.05, max(y), length = 100)
#' h <- 2
#' den <- RIG(x = xx, y = y, k = 200, h = h)
#'
#' ##If scheme for generating gridpoints is unknown
#' y <- rexp(50, 1)
#' h <- 3
#' den <- RIG(y = y, k = 90, h = h)
#'
#'\dontrun{
#'##If user do not mention the number of grid points
#'y <- rexp(23, 1)
#'xx <- seq(min(y) + 0.05, max(y), length = 90)
#'#any bandwidth can be used
#'require(KernSmooth)
#'h <- dpik(y)
#'den <- RIG(x = xx, y = y, h = h)
#'}
#'#if bandwidth is missing
#'y <- rexp(100, 1)
#'xx <- seq(min(y) + 0.05, max(y), length = 100)
#'den <- RIG(x = xx, y = y, k = 90)
#' @return \item{x}{grid points}
#'         \item{y}{estimated values of density}
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Scaillet, O. 2004. Density estimation using inverse and reciprocal inverse Gaussian kernels. 	 \emph{Nonparametric Statistics}, \strong{16}, 217-226.
#' @seealso To examine RIG density plot see \code{\link{plot.RIG}} and for Mean Squared Error \code{\link{mse}}. Similarly, for Laplace kernel \code{\link{Laplace}}.
#' @export
RIG <- function(x = NULL, y, k = NULL, h = NULL){
   n <- length(y)
   if(is.null(x))
      x = seq(min(y) + 0.05, max(y), length = k)
   if(is.null(k))
      k = length(y)
   if(is.null(h))
      h = 0.79 * IQR(y) * length(y) ^ (-1/5)
   fhat <- rep(0, k)
   KRIG <- matrix(rep(0, k * n), ncol = k)

   for(j in 1:k) {
      for(i in 1:n) {

         KRIG[i, j] <- 1/(sqrt(2 * pi * h * y[i])) * exp(((-1 *(y[i] - h))/(2 * h) * (x[j]/(y[i] - h) - 2 +(y[i] - h)/x[j])))
      }
      fhat<- colMeans(KRIG)
   }
   results <- list(x = x, y = fhat)
   class ( results ) <-c('list', 'RIG')
   results
}

#' Density Plot by Reciprocal Inverse Gaussian kernel
#'
#' Plot density by using Reciprocal Inverse Gaussian Kernel.
#' @param x an object of class "RIG"
#' @param \dots Not presently used in this implementation
#' @import graphics
#' @import stats
#' @examples
#' y <- rexp(200, 1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' xx <- seq(min(y) + 0.05, max(y), length = 200)
#' den <- RIG(x = xx, y = y, k = 200, h = h)
#' plot(den, type = "l")
#'
#'
#' ##other details can also be added
#' y <- rexp(200, 1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' den <- RIG(x = xx, y = y, k = 200, h = h)
#' plot(den, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
#'
#' ## To add true density along with estimated
#' d1 <- density(y, bw = h)
#' lines(d1, type = "p", col = "red")
#' legend("topright", c("Real Density", "Density by RIG Kernel"),
#' col = c("red", "black"), lty = c(1, 2))
#' @return nothing
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Scaillet, O. 2004. Density estimation using inverse and reciprocal inverse Gaussian kernels. 	 \emph{Nonparametric Statistics}, \strong{16}, 217-226.
#' @seealso To examine RIG estimated values for density see \code{\link{RIG}} and for Mean Squared Error \code{\link{mse}}. Similarly, for plot of Laplace kernel \code{\link{plot.Laplace}}.
#' @export
plot.RIG <- function(x,...) {
   plot(x$x, x$y,...)
}
#' Estimate Density Values by Laplace kernel
#'
#' Estimated Kernel density values by using Laplace Kernel.
#' @details Laplace kernel is developed by Khan and Akbar. Kernel is developed by using Chen's idea. Laplace kernel is;
#' \deqn{K_{Laplace\left(x,h^{\frac{1}{2}}\right)} (u)=\frac{1}{2\sqrt h}exp \left(-\frac{|{u-x}|}{\sqrt h}\right)}

#' @param y a numeric vector of positive values.
#' @param x scheme for generating grid points
#' @param k gird points.
#' @param h the bandwidth
#' @import stats
#' @examples
#' #Data can be simulated or real data
#' ## Number of grid points "k" should be at least equal to the data size.
#' ### If user define the generating scheme of gridpoints than number of gridpoints should
#' ####be equal or greater than "k"
#' ###### otherwise NA will be produced.
#' y <- rexp(100, 1)
#' xx <- seq(min(y) + 0.05, max(y), length = 100)
#' h <- 2
#' den <- Laplace(x = xx, y = y, k = 200, h = h)
#'
#' ##If scheme for generating gridpoints is unknown
#' y <- rexp(50, 1)
#' h <- 3
#'den <- Laplace(y = y, k = 90, h = h)
#'
#'##If user do not mention the number of grid points
#'y <- rexp(23, 1)
#'xx <- seq(min(y) + 0.05, max(y), length = 90)
#'
#'\dontrun{
#'#any bandwidth can be used
#'require(KernSmooth)
#'h <- dpik(y)
#'den <- Laplace(x = xx, y = y, h = h)
#'}
#'
#'#if bandwidth is missing
#'y <- rexp(100, 1)
#'xx <- seq(min(y) + 0.05, max(y), length = 100)
#'den <- Laplace(x = xx, y = y, k = 90)
#' @return \item{x}{grid points}
#'         \item{y}{estimated values of density}
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Khan, J. A.; Akbar, A. Density Estimation by Laplace Kernel. \emph{Working paper,  Department of Statistics, Bahauddin Zakariya University, Multan, Pakistan.}
#' @seealso To examine Laplace density plot see \code{\link{plot.Laplace}} and for Mean Squared Error \code{\link{mse}}. Similarly, for RIG kernel \code{\link{RIG}}.
#' @export
Laplace <- function(x = NULL, y, k = NULL, h = NULL){
   n <- length(y)
   if(is.null(x))
      x = seq(min(y) + 0.05, max(y), length = k)
   if(is.null(k))
      k = length(y)
   if(is.null(h))
      h = 0.79 * IQR(y) * length(y) ^ (-1/5)
   KLaplace <- matrix(rep(0, k * n), ncol = k)
   ######Lap#####
   for(j in 1:k) {
      for(i in 1:n) {

         KLaplace[i, j] <-(1/(2 * sqrt(h))) * exp( - (abs(y[i] - (x[j]))) / sqrt(h))
      }
   }
      fhat <- fhat<- colMeans(KLaplace)
      results <- list(x = x,
                      y = fhat)
      class ( results ) <-c('list', 'Laplace')
      results
   }

#' Density Plot by Laplace kernel
#'
#' Plot density by using Laplace Kernel.
#' @param x an object of class "Laplace"
#' @param \dots Not presently used in this implementation
#' @import graphics
#' @import stats
#' @examples
#' y <- rexp(100, 1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' xx <- seq(min(y) + 0.05, max(y), length = 100)
#' den <- Laplace(x = xx, y = y, k = 100, h = h)
#' plot(den, type = "l")
#'
#'
#' ##other details can also be added
#' y <- rexp(100, 1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' den <- Laplace(x = xx, y = y, k = 100, h = h)
#' plot(den, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
#'
#' ## To add true density along with estimated
#' d1 <- density(y, bw = h)
#' lines(d1, type = "p", col = "red")
#' legend("topright", c("Real Density", "Density by RIG Kernel"),
#' col = c("red", "black"), lty = c(1, 2))

#' @return nothing
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Khan, J. A.; Akbar, A. Density Estimation by Laplace Kernel. \emph{Working paper,  Department of Statistics, Bahauddin Zakariya University, Multan, Pakistan.}
#' @seealso To examine Laplace estimated values for density see \code{\link{Laplace}} and for Mean Squared Error \code{\link{mse}}. Similarly, for plot of Laplace kernel \code{\link{plot.RIG}}.

#' @export
plot.Laplace <- function(x,...) {
   plot(x$x, x$y,...)
}

#' Calculate Mean Squared Error( MSE) by using different Kernels
#'
#' This function calculates the mean squared error (MSE) by using user specified kernel.This function is same as provided in package "DELTD". For details see \url{https://CRAN.R-project.org/package=DELTD}.
#' @param kernel type of kernel which is to be used
#' @param type mention distribution of vector.If exponential distribution then use \code{"Exp"}.
#'     If use gamma distribution then use \code{"Gamma"}.If Weibull distribution then use \code{"Weibull"}.
#' @import stats
#' @seealso This is also available in \pkg{DELTD}
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references  \url{https://CRAN.R-project.org/package=DELTD}
#'
#' @examples
#' y <- rexp(100, 1)
#' xx <- seq(min(y) + 0.05, max(y), length = 500)
#' h <- 2
#' gr <- Laplace(x = xx, y = y, k = 200, h = h)
#' mse(kernel = gr, type = "Exp")
#' ## if distribution is other than mentioned \code{type} is used then NaN will be produced.
#' \dontrun{
#' mse(kernel = gr, type ="Beta")
#' }
#' @return Mean Squared Error (MSE)
#' @export
mse<-function(kernel,type){
   ftrue<-switch(type,
                 Exp = dexp(kernel$x, (1 / mean(kernel$x))),
                 Gamma = dgamma(kernel$x, (mean(kernel$x) / (var(kernel$x) / mean(kernel$x))), (var(kernel$x) / mean(kernel$x))),
                 Weibull = dweibull(kernel$x, ((sd(kernel$x) / mean(kernel$x)) ^ - 1.806), scale = 1, log = FALSE)
   )
   return(mean((ftrue-kernel$y)^2))#mse of fhat w.r.t. the true density
}
 #'
 #' Calculate Bandwidth proposed by Silverman for non-normal data.
 #' @param y a numeric vector of positive values.
 #' @import graphics
 #' @import stats
 #' @author Javaria Ahmad Khan, Atif Akbar.
 #' @references Silverman, B. W. 1986. \emph{Density Estimation}. Chapman & Hall/ CRC, London.
 #' @examples
 #' y <- rexp(10, 1)
 #'  NSR(y)
 #' @return h
 #' @export
 NSR <- function(y){
   return(0.79 * IQR(y) * length(y) ^ (-1/5))
 }
