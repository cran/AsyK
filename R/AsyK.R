#' Plot Density by RIG kernel.
#'
#' Plot Kernel density by using Resiprocal Inverse Gaussian Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @import graphics
#' @import stats
#' @examples
#' y<-rexp(23,1)
#' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
#' graphrig(y,80,h)
#' @export

graphrig<-function(y,k,h){
  n <- length(y)
  x <- seq(min(y) + 0.05, max(y), length = k)

  fhat <- rep(0, k)
  KRIG <- matrix(rep(0, k * n), ncol = k)

  for(j in 1:k) {
    for(i in 1:n) {

      KRIG[i, j] <- 1/(sqrt(2 * pi * h * y[i])) * exp(((-1 *(y[i] - h))/(2 * h) * (x[j]/(y[i] - h) - 2 +(y[i] - h)/x[j])))
    }
    fhat[j] <- 1/n * (sum(KRIG[, j]))
    d1<-density(y,bw=h)
    plot(x, fhat, type = "l",ylab = "Density Function", lty = 1, xlab = "Time",ylim=c(min(fhat),max(fhat)))
    lines(d1,type="p",col="red")
    legend("topright", c("Real Density", "Density by  RIG  Kernel"),
           col=c("red", "black"), lty=c(1,2))
  }}
#' Plot Density by Laplace kernel.
#'
#' Plot Kernel density by using Laplace Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @import graphics
#' @import stats
#' @examples
#' y<-rexp(23,1)
#' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
#' graphLap(y,80,h)
#' @export
#'
 graphLap<-function(y,k,h){
n <- length(y)
x <- seq(min(y) + 0.05, max(y), length = k)

fhat <- rep(0, k)
KLaplace <- matrix(rep(0, k * n), ncol = k)

for(j in 1:k) {
  for(i in 1:n) {

    KLaplace[i, j] <-(1/(2*sqrt(h)))*exp(-(abs(y[i]-(x[j])))/sqrt(h))
  }
  fhat[j] <- 1/n * (sum(KLaplace[, j]))
  d1<-density(y,bw=h)
  plot(x,fhat, type="l", ylab = "Density Function",ylim=c(min(fhat),(max(fhat)+0.2))  ,lty = 2, xlab = "Time")
  lines(d1,type="l",col="green")
  legend("topright", c("Real Density", "Density by Laplace Kernel"),
         col=c("green", "black"), lty=c(1,2))
}}

 #' Calculate Mean Square Error( MSE) when RIG kernel is used.
 #'
 #' Calculate MSE by using Resiprocal Inverse Gaussian Kernel.
 #' @param y a numeric vector of positive values.
 #' @param k gird points.
 #' @param h the bandwidth
 #' @param type mention distribution of vector.If exponential distribution then use "Exp".
 #'     if use gamma distribution then use "Gamma".
 #' @import graphics
 #' @import stats
 #' @examples
 #' y<-rexp(100,1)
 #' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
 #' mserig(y,200,h,"exp")
 #' @return MSE
 #' @export
 #'
 mserig<-function(y,k,h,type){
   n <- length(y)
   x <- seq(min(y) + 0.05, max(y), length = k)
   ftrue<-switch(type,

                 Exp = dexp(x,(1/mean(x))),

                 Gamma = dgamma(x,(mean(x)/(var(x)/mean(x))),(var(x)/mean(x)))

   )

   fhat <- rep(0, k)
   KRIG <- matrix(rep(0, k * n), ncol = k)

   for(j in 1:k) {
     for(i in 1:n) {

       KRIG[i, j] <- 1/(sqrt(2 * pi * h * y[i])) * exp(((-1 *(y[i] - h))/(2 * h) * (x[j]/(y[i] - h) - 2 +(y[i] - h)/x[j])))
     }# 1st loop end
     fhat[j] <- 1/n * (sum(KRIG[, j]))

   }#2nd loop end
   return(mean((ftrue-fhat)^2))#mse of fhat w.r.t. the true density
 }#function end


 #' Calculate Mean Square Error( MSE) when Laplace Kernel is used.
 #'
 #' Calculate MSE by using Laplace Kernel.
 #' @param y a numeric vector of positive values.
 #' @param k gird points.
 #' @param h the bandwidth
 #' @param type mention distribution of vector.If exponential distribution then use "Exp".
 #'     if use gamma distribution then use "Gamma".
 #' @import graphics
 #' @import stats
 #' @examples
 #' y<-rexp(100,1)
 #'  h<-0.79 * IQR(y) * length(y) ^ (-1/5)
 #'  mselap(y,200,h,"exp")
 #' @return MSE
 #' @export
 #'
 mselap<-function(y,k,h,type){
   n <- length(y)
   x <- seq(min(y) + 0.05, max(y), length = k)
   ftrue<-switch(type,

                 Exp = dexp(x,(1/mean(x))),

                 Gamma = dgamma(x,(mean(x)/(var(x)/mean(x))),(var(x)/mean(x)))
   )

   fhat <- rep(0, k)
   KLaplace <- matrix(rep(0, k * n), ncol = k)

   for(j in 1:k) {
     for(i in 1:n) {

       KLaplace[i, j] <-(1/(2*sqrt(h)))*exp(-(abs(y[i]-(x[j])))/sqrt(h))
     }# 1st loop end
     fhat[j] <- 1/n * (sum(KLaplace[, j]))

   }#2nd loop end
   return(mean((ftrue-fhat)^2))#mse of fhat w.r.t. the true density
 }#function end

 #' Calculate Bandwidth.
 #'
 #' Calculate Bandwidth proposed by Silverman for nonnormal data.
 #' @param y a numeric vector of positive values.
 #' @import graphics
 #' @import stats
 #' @examples
 #' y<-rexp(10,1)
 #'  NSR(y)
 #' @return h
 #' @export
 NSR<-function(y){
   return(0.79 * IQR(y) * length(y) ^ (-1/5))
 }

